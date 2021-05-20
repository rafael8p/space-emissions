# -*- coding: utf-8 -*-
# %%
"""Fioletov Multi-source emission calculator"""
import os
import json
import time
from datetime import datetime

import pandas
import geopandas
from shapely.geometry import shape

from eocalc.context import Pollutant
from eocalc.methods.base import DateRange
from eocalc.methods.tools import *
from eocalc.methods.fioletov import *

# %%
# time test
s_t = time.time()
# select region for test
directory = "data/regions"
regions = {}

for filename in os.listdir(directory):
    if (filename.endswith(".geo.json") or filename.endswith(".geojson") and ('germany.geo.json' in filename)):
        with open(f"{directory}/{filename}", 'r') as geojson_file:
            regions[filename] = shape(json.load(geojson_file)["geometry"])

# del regions["europe.geo.json"] # Uncomment if this takes too long...

# Remove regions not covered
regions = {filename: region for filename, region in regions.items() if TropomiMonthlyMeanAggregator.covers(region)}

# %%


# %%