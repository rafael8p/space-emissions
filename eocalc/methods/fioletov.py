# -*- coding: utf-8 -*-
"""Fioletov functions"""
import os.path
import shutil
import gzip
import threading
from datetime import date, timedelta
from urllib.request import urlretrieve

import numpy
from pandas import Series
from shapely.geometry import MultiPolygon, shape
from geopandas import GeoDataFrame, overlay

from eocalc.context import Pollutant
from eocalc.methods.base import EOEmissionCalculator, DateRange, Status
import eocalc.methods.tools as tools

# Local directory we use to store downloaded and decompressed data
LOCAL_DATA_FOLDER = "/media/uba_emis/space_emissions/enrico/"
# Local directory we use to store downloaded and decompressed data
LOCAL_ERA5_FOLDER = "/media/uba_emis/space_emissions/enrico/ERA5/"
# Local directory we use to store downloaded and decompressed data
LOCAL_S5P_FOLDER = "/codede/Sentinel-5P/"
# Online resource used to download TEMIS data on demand
# TODO! implement download tool from Janot
# Scitools.sh/python version
resolution_lon, resolution_lat = 0.2,0.2
class MultiSourceCalculator(EOEmissionCalculator):

    @staticmethod
    def coverage() -> MultiPolygon:
        return shape({'type': 'MultiPolygon',
                      'coordinates': [[[[-180., -60.], [180., -60.], [180., 60.], [-180., 60.], [-180., -60.]]]]})

    @staticmethod
    def minimum_period_length() -> int:
        return 1

    @staticmethod
    def earliest_start_date() -> date:
        return date.fromisoformat('2019-06-01') # hardcapped for initial tests

    @staticmethod
    def latest_end_date() -> date:
        return date.fromisoformat('2019-07-01') # hardcapped for initial tests
        #return (date.today().replace(day=1) - timedelta(days=1)).replace(day=1) - timedelta(days=1)

    @staticmethod
    def supports(pollutant: Pollutant) -> bool:
        return pollutant == Pollutant.NO2



    def run(self, region= MultiPolygon, period= DateRange, pollutant= Pollutant, resolution=(resolution_lon,resolution_lat)) -> dict:
        self._validate(region, period, pollutant)
        self._state = Status.RUNNING
        self._progress = 0  # TODO Update progress below!
        
        # 1. Create a field of sources
        # TODO allow for domain of higher resolution within a coarser set. For example germany at 0.1 with surroundings at 0.2
        grid = self._create_grid(region, resolution[0], resolution[1], snap=True, include_center_cols=True)

        # 2. Read TEMIS data into the grid, use cache to avoid re-reading the file for each day individually
        cache = {}
        for day in period:
            month_cache_key = f"{day:%Y-%m}"
            if month_cache_key not in cache.keys():
                concentrations = tools._read_subset_data(region, tools._assure_data_availability(region,period))
                # value [1/cm²] * TEMIS scale [1] / Avogadro constant [1] * NO2 molecule weight [g] / to [kg] * to [km²]
                cache[month_cache_key] = [x * 10**13 / (6.022 * 10**23) * 46.01 / 1000 * 10**10 for x in concentrations]
                # TODO Correct for pollutant atmosphere lifetime and diurnal variation: pollutant.atmo_lifetime(day, latitude) * pollutant.diurnal_variation(day, instrument)

            # Here, values are actually [kg/km²], but the area [km²] cancels out below
            grid.insert(0, f"{day} {pollutant.name} emissions [kg]", cache[month_cache_key])

        # 3. Clip to actual region and add a data frame column with each cell's size
        grid = overlay(grid, GeoDataFrame({'geometry': [region]}, crs="EPSG:4326"), how='intersection')
        grid.insert(0, "Area [km²]", grid.to_crs(epsg=8857).area / 10 ** 6)  # Equal earth projection

        # 4. Update emission columns by multiplying with the area value and sum it all up
        grid.iloc[:, -(len(period)+3):-3] = grid.iloc[:, -(len(period)+3):-3].mul(grid["Area [km²]"], axis=0)
        grid.insert(1, f"Total {pollutant.name} emissions [kg]", grid.iloc[:, -(len(period)+3):-3].sum(axis=1))
        grid.insert(2, "Umin [%]", numpy.NaN)
        grid.insert(3, "Umax [%]", numpy.NaN)
        grid.insert(4, "Number of values [1]", len(period))
        grid.insert(5, "Missing values [1]", grid.iloc[:, -(len(period)+3):-3].isna().sum(axis=1))
        self._calculate_row_uncertainties(grid, period)  # Replace NaNs in Umin/Umax cols with actual values

        # 5. Add GNFR table incl. uncertainties
        table = self._create_gnfr_table(pollutant)
        total_uncertainty = self._combine_uncertainties(grid.iloc[:, 1], grid.iloc[:, 2])
        table.iloc[-1] = [grid.iloc[:, 1].sum() / 10**6, total_uncertainty, total_uncertainty]

        self._state = Status.READY
        return {self.TOTAL_EMISSIONS_KEY: table, self.GRIDDED_EMISSIONS_KEY: grid}


    # operations for class / or functions to fit emissions using fioletov approach

    #__some__basic__operation__


    #__get_all_plume_functions_


    #add_bias


    #fit_for_domain


    #errors_uncertainties_estimator


    #posteriori_adjustments