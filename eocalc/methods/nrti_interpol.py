import numpy
import copy
from osgeo import gdal
import h5py
from scipy import interpolate
from scipy.signal import triang
import sentinelsat


class GetData:  # only demo purposes do not use this class
    def __init__(self):
        pass

    @staticmethod
    def getnrti(geojsonfile, output: str):
        api = sentinelsat.SentinelAPI(user='', password='',
                                      api_url='https://s5phub.copernicus.eu/dhus')
        footprint = sentinelsat.geojson_to_wkt(sentinelsat.read_geojson(geojsonfile))
        products = api.query(area=footprint, date="[NOW-2DAYS TO NOW]", platformname='Sentinel-5',
                             processingmode="Near real time", producttype='L2__NO2___')
        api.download_all(products, output)
        products = api.query(area=footprint, date="[NOW-2DAYS TO NOW]", platformname='Sentinel-5',
                             processingmode="Near real time", producttype='L2__SO2___')
        api.download_all(products, output)
        products = api.query(area=footprint, date="[NOW-2DAYS TO NOW]", platformname='Sentinel-5',
                             processingmode="Near real time", producttype='L2__HCHO__')
        api.download_all(products, output)
        products = api.query(area=footprint, date="[NOW-2DAYS TO NOW]", platformname='Sentinel-5',
                             processingmode="Near real time", producttype='L2__CO____')
        api.download_all(products, output)
        products = api.query(area=footprint, date="[NOW-2DAYS TO NOW]", platformname='Sentinel-5',
                             processingmode="Near real time",
                             producttype='L2__CH4__')  # funzt nicht bzw. daten gab es nicht so schnell
        api.download_all(products, output)
        products = api.query(area=footprint, date="[NOW-2DAYS TO NOW]", platformname='Sentinel-5',
                             processingmode="Near real time", producttype='L2__O3___')


class L2NRT:
    def __init__(self):
        pass

    @classmethod
    def multall(cls, dat, gains):
        g = 0
        for i in gains:
            if i == numpy.nan:
                i = 1
            if i == numpy.NINF:
                i = 1
            if i == numpy.Inf:
                i = 1
            dat[:, g] = dat[:, g] * i
            g += 1
        return dat

    @classmethod
    def gained(cls, cubi, gains):
        cubi_gained = cls.multall(cubi, gains)
        return cubi_gained

    @classmethod
    def triang_normed(cls, N):  # Dreicksfilter mit der Flaeche 1 und edge != 0
        if N % 2 == 0:
            print("you need to input an odd number")
            exit(-1)
        filt = triang(N)
        filt = filt / ((N + 1) / 2)
        return filt

    @classmethod
    def average_columns(cls, data):
        dim = data.shape
        data = data / 1.0
        # daten = numpy.ndarray((dim[2]), dtype='float')
        dat = numpy.nanmedian(data, axis=0)
        return dat

    @classmethod
    def Im_cheating(cls, data, filter_normed):
        dat = cls.average_columns(data)
        dat = numpy.abs(numpy.convolve(dat, filter_normed, 'same'))
        return dat

    @classmethod
    def stanze_s5p(cls, data):
        data[data > 390000000] = 0
        fido = cls.triang_normed(11)  # 55 statt 11
        avg_dat = cls.average_columns(data)
        avg_dat_filt = cls.Im_cheating(data, fido)
        gains = abs(avg_dat_filt / avg_dat)
        numpy.place(gains, gains == numpy.nan, 1)
        destriped = cls.gained(data, gains)
        return destriped, gains

    @classmethod
    def read_l2_nrt(cls, file: str) -> (numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray):
        lss = h5py.File(file, 'r')
        lat = lss['PRODUCT']['latitude'][0, :, :]
        long = lss['PRODUCT']['longitude'][0, :, :]
        if "L2__CO" in file:
            data = lss['PRODUCT']['carbonmonoxide_total_column'][0, :, :]
            dataprec = lss['PRODUCT']['carbonmonoxide_total_column_precision'][0, :, :]
        elif "L2__HCHO" in file:
            data = lss['PRODUCT']['formaldehyde_tropospheric_vertical_column'][0, :, :]
            dataprec = lss['PRODUCT']['formaldehyde_tropospheric_vertical_column_precision'][0, :, :]
        elif "L2__NO2" in file:
            data = lss['PRODUCT']['nitrogendioxide_tropospheric_column'][0, :, :]
            dataprec = lss['PRODUCT']['nitrogendioxide_tropospheric_column_precision'][0, :, :]
        elif "L2__SO2" in file:
            data = lss['PRODUCT']['sulfurdioxide_total_vertical_column'][0, :, :]
            dataprec = lss['PRODUCT']['sulfurdioxide_total_vertical_column_precision'][0, :, :]
        elif "L2__O3" in file:
            data = lss['PRODUCT']['ozone_total_vertical_column'][0, :, :]
            dataprec = lss['PRODUCT']['ozone_total_vertical_column_precision'][0, :, :]
        elif "L2__CH4" in file:
            data = None
            dataprec = None
            lat = None
            long = None
        else:
            data = None
            dataprec = None
        quality = lss['PRODUCT']['qa_value'][0, :, :]
        return lat, long, data, quality, dataprec

    @classmethod
    def filter_vals(cls, data: numpy.ndarray) -> (numpy.ndarray, numpy.ndarray):
        data_indices = numpy.where(data > 1e30)
        data[data < 0] = 0
        data[data > 1e30] = 0
        return data_indices, data

    @classmethod
    def filter_vals_qa(cls, data: numpy.ndarray, qaval: numpy.ndarray) -> (numpy.ndarray, numpy.ndarray):
        data_indices = numpy.where(qaval > 50)
        data[data < 0] = 0
        data[data > 1e30] = 0
        return data_indices, data

    @classmethod
    def param_interpol(cls, lat: numpy.ndarray, long: numpy.ndarray) -> (float, float, numpy.ndarray):
        b = numpy.abs(long - numpy.roll(long, 1, axis=0))
        a = numpy.abs(long - numpy.roll(long, 1))
        c = numpy.abs(lat - numpy.roll(lat, 1, axis=0))
        d = numpy.abs(lat - numpy.roll(lat, 1))
        b50 = numpy.median(b[1:, :])
        a50 = numpy.median(a[:, 1:])
        c50 = numpy.median(c[1:, :])
        d50 = numpy.median(d[:, 1:])
        longstep = (a50 + b50) / 2.0
        latstep = (c50 + d50) / 2.0
        long_new = numpy.arange(numpy.min(long), numpy.max(long), longstep)  # times factor 2?
        lat_new = numpy.flip(
            numpy.arange(numpy.min(lat), numpy.max(lat), latstep))  # times factor 2? to reduce the gridsize?
        grid_long_new, grid_lat_new = numpy.meshgrid(long_new, lat_new)
        interpolgrid = numpy.asarray([grid_long_new, grid_lat_new])
        return longstep, latstep, interpolgrid

    @classmethod
    def schreibebsqsingle(cls, numpyarray: numpy.ndarray, out: str):
        if len(numpyarray.shape) != 2:
            print("no valid 2d-dataset!")
            return -1
        form = "ENVI"
        driver = gdal.GetDriverByName(form)
        dimat = numpy.shape(numpyarray)
        dataset_out = driver.Create(out, dimat[1], dimat[0], 1, gdal.GDT_Float32)
        dataset_out.GetRasterBand(1).WriteArray(numpyarray)
        dataset_out = None

    # interpolate the nulled data
    @classmethod
    def cubic_2_grid(cls, interplogrid: numpy.ndarray, longstep: float, latstep: float, long: numpy.ndarray,
                     lat: numpy.ndarray,
                     data: numpy.ndarray, outfile: str, subname: str):
        mapifo = "map info = {Geographic Lat/Lon,1, 1," + ' ' + str(numpy.min(interplogrid[0, :, :])) + ', ' + str(
            numpy.max(interplogrid[1, :, :])) + ', ' + str(longstep) + ', ' + str(latstep) + ",WGS-84}\n"
        print(mapifo)
        d = numpy.swapaxes(numpy.asarray([long.flatten(), lat.flatten()]), 0, 1)
        gridcd = interpolate.griddata(d, data.flatten(), (interplogrid[0, :], interplogrid[1, :]), method='cubic',
                                      fill_value=0)
        gridcd = numpy.nan_to_num(gridcd)
        namm = outfile + '_' + subname
        cls.schreibebsqsingle(gridcd, namm)
        open(namm + '.hdr', 'a').write(mapifo)
        return None

    @classmethod
    def cubic_2_grid_qa(cls, interplogrid: numpy.ndarray, longstep: float, latstep: float, long: numpy.ndarray,
                        lat: numpy.ndarray,
                        data: numpy.ndarray, qaindex: numpy.ndarray, outfile: str, subname: str):
        mapifo = "map info = {Geographic Lat/Lon,1, 1," + ' ' + str(numpy.min(interplogrid[0, :, :])) + ', ' + str(
            numpy.max(interplogrid[1, :, :])) + ', ' + str(longstep) + ', ' + str(latstep) + ",WGS-84}\n"
        print(mapifo)
        d = numpy.swapaxes(numpy.asarray([long[qaindex].flatten(), lat[qaindex].flatten()]), 0, 1)
        gridcd = interpolate.griddata(d, data[qaindex].flatten(), (interplogrid[0, :], interplogrid[1, :]),
                                      method='cubic',
                                      fill_value=0)
        gridcd = numpy.nan_to_num(gridcd)
        namm = outfile + '_' + subname
        cls.schreibebsqsingle(gridcd, namm)
        open(namm + '.hdr', 'a').write(mapifo)
        return None

    @classmethod
    def run(cls, infile: str, subname: str, flag: bool):
        bname = infile.split('.')[0] + '_gisproc'
        lat, long, data, quality, dataprec = cls.read_l2_nrt(infile)
        qaarr = quality[quality > 0.5]
        if (len(qaarr) / len(quality.flat)) > 0.5 and flag == True:
            data_dstr = cls.stanze_s5p(data)
            print('Staenz Destriping')  # please dont use staenz destriping on a broken L2 dset
        else:
            data_dstr = [data, data]
        longstep, latstep, interpolgrid = cls.param_interpol(lat, long)
        aux, filleddat = cls.filter_vals(data_dstr[0])
        filleddat = numpy.nan_to_num(filleddat)
        cls.cubic_2_grid(interpolgrid, longstep, latstep, long, lat, filleddat, bname, subname)
        return None

    @classmethod
    def run_qa(cls, infile: str, subname: str, flag: bool):
        bname = infile.split('.')[0] + '_gisproc'
        lat, long, data, quality, dataprec = cls.read_l2_nrt(infile)
        longstep, latstep, interpolgrid = cls.param_interpol(lat, long)
        qaarr = quality[quality > 0.5]
        if (len(qaarr) / len(quality.flat)) > 0.5 and flag == True:
            data_dstr = cls.stanze_s5p(data)
            print('Staenz Destriping')  # please dont use staenz destriping on a broken L2 dset
        else:
            data_dstr = [data, data]
        indices, filleddat = cls.filter_vals_qa(data_dstr[0], quality)
        filleddat = numpy.nan_to_num(filleddat)
        cls.cubic_2_grid_qa(interpolgrid, longstep, latstep, long, lat, filleddat, indices, bname, subname)
        return None
