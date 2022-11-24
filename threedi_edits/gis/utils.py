# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:52:59 2021
GDAL utils from raster-tools.

"""
from osgeo import ogr
from osgeo import gdal
from osgeo import gdal_array


def create(
    array,
    band=False,
    geotransform=None,
    spatial_reference=None,
    nodata_value=None,
):
    """
    Create and return a gdal dataset.

    :param array: A numpy array.
    :param geo_transform: 6-tuple of floats
    :param projection: wkt projection string
    :param no_data_value: integer or float

    This is the fastest way to get a gdal dataset from a numpy array, but
    keep a reference to the array around, or a segfault will occur. Also,
    don't forget to call FlushCache() on the dataset after any operation
    that affects the array.
    """
    # prepare dataset name pointing to array
    datapointer = array.ctypes.data
    if len(array.shape) == 2:
        lines, pixels = array.shape
        bands = 1
    else:
        bands, lines, pixels = array.shape
    datatypecode = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype.type)
    datatype = gdal.GetDataTypeName(datatypecode)
    if len(array.strides) == 2:
        lineoffset, pixeloffset = array.strides
        bandoffset = 0
    else:
        bandoffset, lineoffset, pixeloffset = array.strides

    dataset_name_template = (
        "MEM:::"
        "DATAPOINTER={datapointer},"
        "PIXELS={pixels},"
        "LINES={lines},"
        "BANDS={bands},"
        "DATATYPE={datatype},"
        "PIXELOFFSET={pixeloffset},"
        "LINEOFFSET={lineoffset},"
        "BANDOFFSET={bandoffset}"
    )
    dataset_name = dataset_name_template.format(
        datapointer=datapointer,
        pixels=pixels,
        lines=lines,
        bands=bands,
        datatype=datatype,
        pixeloffset=pixeloffset,
        lineoffset=lineoffset,
        bandoffset=bandoffset,
    )

    # access the array memory as gdal dataset
    dataset = gdal.Open(dataset_name, gdal.GA_Update)

    # set additional properties from kwargs
    if geotransform is not None:
        dataset.SetGeoTransform(geotransform)
    if spatial_reference is not None:
        dataset.SetProjection(spatial_reference)
    if nodata_value is not None:
        dataset.GetRasterBand(1).SetNoDataValue(nodata_value)
        # for i in range(len(array)):
        # dataset.GetRasterBand(i + 1).SetNoDataValue(no_data_value)

    if band:
        return dataset, dataset.GetRasterBand(1)
    else:
        return dataset


class Dataset(object):
    """
    Usage:
        >>> with Dataset(array) as dataset:
        ...     # do gdal things.
    """

    def __init__(self, array, **kwargs):
        self.array = array
        self.dataset = create(array, False, **kwargs)

    def __enter__(self):
        return self.dataset

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.dataset.FlushCache()


class Layer:
    """
    Usage:
        >>> with Layer(geometry) as layer:
        ...     # do ogr things.
    """

    def __init__(self, geometry):
        driver = ogr.GetDriverByName(str("Memory"))
        self.data_source = driver.CreateDataSource("")
        sr = geometry.GetSpatialReference()
        self.layer = self.data_source.CreateLayer(str(""), sr)
        layer_defn = self.layer.GetLayerDefn()
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(geometry)
        self.layer.CreateFeature(feature)

    def __enter__(self):
        return self.layer

    def __exit__(self, exc_type, exc_value, traceback):
        pass


DRIVERS = {}
for i in range(ogr.GetDriverCount()):
    drv = ogr.GetDriver(i)
    if (
        drv.GetMetadataItem("DCAP_VECTOR") == "YES"
        and drv.GetMetadataItem("DCAP_CREATE") == "YES"
    ):
        driver = drv.GetMetadataItem("DMD_EXTENSION")
        if driver:
            DRIVERS["." + driver] = drv.GetName()
