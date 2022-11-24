# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 18:20:00 2020

@author: chris.kerklaan

Raster wrapper used for gis programming
This raster wrapper has as peculiarity that for every function the raster is 
copied/loaded into memory. Hence it always has to be written. I.e.,
it cannot be edited in place.

TODO:
    1. Buffer toevoegen aan een raster --> Misschien iets voor rastergroup ThreeDi
    4. Optimalizeren van blocks
    5. Testen voor dit hele script, testen moeten ook kijken naar geheugen issuses
    
Notes:
    Memory handling:
        1. Numpy calculations are done in tiles as much as possible - speed/memory efficiency
        2. This is also why rasters are edited in place and not returned as new rasters.
        
    Bottleneck:
        1. Rasters are loaded into memory if needed, before the use of functions 
        2. Due to the fact that the rasters are loaded into memory, performance depends on your system
        
    Functions:
        1. Each function handles the object in tiles.
        2. After each function the old file is released and the file is referenced.
        3. Functions outside of the class or either extremely complicated or
           also use multiple external functions
    
    Iterating:
        1. One can iterate over this raster using for tile in raster:
           Blocksizes are optimzed and arrays must be closed after using
           Recommeded code:
               output = raster.copy(shell=True)
               for tile in raster:
                   array = tile
                   # do something with array
                   output.array = array, *tile.location
                   

"""
# First-party import
import logging
import os
import pathlib
import importlib
import tempfile
import numpy as np

# Third-party imports
from osgeo import ogr, gdal, gdal_array

# Third-party imports optional

# osgeo - base package
HAS_SCIPY = importlib.util.find_spec("scipy") is not None
if HAS_SCIPY:
    import scipy

# Local imports
from .vector import Vector
from .geotransform import GeoTransform
from .utils import Layer, Dataset
from .geometry import POLYGON
from .spatial_reference import SpatialReference
from threedi_edits.utils.progress import Progress
from threedi_edits.utils.project import Classes, Functions

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)


# GLOBALS
# Drivers
DRIVER_OGR_SHP = ogr.GetDriverByName("Esri Shapefile")
DRIVER_GDAL_GTIFF = gdal.GetDriverByName("GTiff")
DRIVER_GDAL_MEM = gdal.GetDriverByName("MEM")

# Options
WRITE_OPTIONS = ["COMPRESS=DEFLATE", "NUM_THREADS=4", "ZLEVEL=9", "TILED=YES"]

# Gdal/ogr band/field mapping
TYPE_MAPPING = {
    gdal.GDT_Byte: ogr.OFTInteger,
    gdal.GDT_UInt16: ogr.OFTInteger,
    gdal.GDT_Int16: ogr.OFTInteger,
    gdal.GDT_UInt32: ogr.OFTInteger,
    gdal.GDT_Int32: ogr.OFTInteger,
    gdal.GDT_Float32: ogr.OFTReal,
    gdal.GDT_Float64: ogr.OFTReal,
    gdal.GDT_CInt16: ogr.OFTInteger,
    gdal.GDT_CInt32: ogr.OFTInteger,
    gdal.GDT_CFloat32: ogr.OFTReal,
    gdal.GDT_CFloat64: ogr.OFTReal,
}
TYPE_MAPPING_INPUT = {
    float: gdal.GDT_Float32,
    int: gdal.GDT_Int32,
}
# memory counter
_loc_ras = 0

# logging
logger = logging.getLogger(__name__)


class Raster:

    instance = "gis.raster.Raster"

    def __init__(
        self,
        raster_input,
        name=None,
        optimize_blocksize=True,
        write=0,
        **kwargs,
    ):
        """

        Raster -- info

        Raster object for general raster GIS processing and reading
        Raster has functions such as:
            - read(geometry)
            - clip
            - merge
            - smooth
            - resample
            - polygonizing
            - push_vector
            - align
            - fill_nodata
            - reproject

        It also uses 'raster looping':
            Iterating:
                1. One can iterate over this raster using for tile in raster:
                    Blocksizes are optimzed and arrays must be closed after using
                    Recommended code:
                        output = raster.copy(shell=True)
                        for tile in raster:
                            array = tile.array
                            # do something with array
                            output.array = array, *tile.location
                        ouput.write(path)

            Inflated iterating with a buffer of 1:
                for tile, original_extent in raster.inflated_tiles(1, True):
                    tile.array
                    # do something with array
                    tile.clip_bbox(*original_extent)

        Use 'help()' for more info on the functions

        params:
            raster_input: str, pathlib.path or gdal.ds
            name: if None, extracted from path
            blocksize: blocksize for raster looping
            optimize_blocksize: if False, uses a standard blocksize
            **kwargs: geotransform, spatial_reference and nodata

        options:
            self.return_raster: Boolean Returns a Raster object at every function
                                instead of inplace changes
            self.nan_array: Boolean Sets all nodata to nan, usefull for raster calculations


        """
        if isinstance(raster_input, pathlib.PurePath):
            ds, name = self.from_path(raster_input, name, write)
        elif type(raster_input) is str:
            ds, name = self.from_path(raster_input, name, write)
        else:
            ds = raster_input

        if not name:
            name = "No name"

        self.ds = ds
        self.name = name

        # options
        self.return_raster = False
        self.nan_array = True
        self.optimize_blocksize = optimize_blocksize

        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def from_array(cls, array, name=None, **kwargs):
        """import from an array raster array"""
        cls.type = "array"
        if len(array.shape) != 2:
            raise ValueError("limited to 2d arrays, check array.shape")
        return cls(gdal_array.OpenArray(array), name, **kwargs)

    @classmethod
    def from_scratch(
        cls,
        name,
        geotransform,
        columns,
        rows,
        nodata_value=-9999,
        data_type=gdal.GDT_Float32,
        epsg=28992,
        memory=True,
    ):
        out, _ = create_source(columns, rows, data_type, memory)

        return cls(
            out,
            name=name,
            geotransform=geotransform,
            spatial_reference=epsg,
            nodata_value=nodata_value,
        )

    def from_path(self, path, name, write=0):
        if not os.path.exists(path):
            raise ValueError(f"Path does not exist: {path}")

        self.type = "path"
        if not name:
            name = pathlib.Path(path).stem
        return gdal.Open(str(path), write), name

    def __iter__(self):
        generator = tile_generator(
            self.ds, self.optimize_blocksize, self.blocksize[0], self.blocksize[1]
        )
        for block, location in generator:
            raster = Raster(block, name=f"{location} - {self.name}")
            raster.location = location
            yield raster

    def tiles_inflated(self, buffer_size=1, original_extent=False):
        """Generator for inflated blocks projection"""
        for or_ext, tile in zip(self.tile_extents, self):
            buffered = tile.extent_vector.buffer(buffer_size)

            clipped = self.clip_bbox(*buffered.extent, name=tile.name)
            if original_extent:
                yield clipped, or_ext
            else:
                yield clipped

    def __len__(self):
        return len(self.tiles)

    def __setitem__(self, key, value):
        x, y = key
        return self.band.WriteArray(value, x, y)

    def __getitem__(self, key):
        x, y = key
        return self.band.ReadAsArray(x, y, self.w, self.h)

    def __str__(self):
        return self.name

    def __repr__(self):
        settings = self.settings
        settings["size"] = self.size
        settings["resolution"] = self.resolution
        return f"({self.instance}) {self.name}: {self.settings}"

    @property
    def array(self):
        array = self.band.ReadAsArray()
        if self.nan_array:
            return self.nan(array)
        return array

    @array.setter
    def array(self, value: (np.ndarray, 0, 0)):
        """Takes '(numpy array, x, y)' to set the array"""
        if type(value) != tuple:
            x = 0
            y = 0
            array = value
        else:
            array, x, y = value

        self.band.WriteArray(array, x, y)
        self.band.FlushCache()

    @property
    def name(self):
        if not hasattr(self, "_name"):
            return "no-name"
        else:
            return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def geotransform(self):
        return GeoTransform(self.ds.GetGeoTransform())

    @geotransform.setter
    def geotransform(self, value):
        self.ds.SetGeoTransform(value)

    @property
    def nodata_value(self):
        return self.band.GetNoDataValue()

    @nodata_value.setter
    def nodata_value(self, value):
        self.band.SetNoDataValue(value)

    @property
    def spatial_reference(self):
        return SpatialReference(wkt=self.ds.GetProjection())

    @spatial_reference.setter
    def spatial_reference(self, sr_input):
        """input is either epsg or sr object"""
        if type(sr_input) == int:
            self.ds.SetProjection(SpatialReference.from_epsg(sr_input).wkt)
        else:
            self.ds.SetProjection(sr_input.wkt)

    @property
    def metadata(self):
        return self.ds.GetMetadata()

    @metadata.setter
    def metadata(self, metadata: dict):
        self.ds.SetMetadata(metadata)

    @property
    def band(self):
        return self.ds.GetRasterBand(1)

    @property
    def epsg(self):
        return self.spatial_reference.epsg

    @property
    def mask(self):
        return ~np.isnan(self.array)

    @property
    def data_type(self):
        return self.band.DataType

    @property
    def data_type_name(self):
        return gdal.GetDataTypeName(self.data_type)

    @property
    def np_data_type(self):
        return gdal_array.GDALTypeCodeToNumericTypeCode(self.data_type)

    @property
    def columns(self):
        return self.band.XSize

    @property
    def rows(self):
        return self.band.YSize

    @property
    def pixels(self):
        return self.columns * self.rows

    @property
    def resolution(self):
        return {"width": self.geotransform[1], "height": self.geotransform[5]}

    @property
    def extent(self):
        return self.get_extent()

    @property
    def extent_geometry(self):
        return ogr.CreateGeometryFromWkt(
            POLYGON.format(x1=self.minx, x2=self.maxx, y1=self.miny, y2=self.maxy)
        )

    @property
    def extent_vector(self):
        v = Vector.from_scratch("tile", ogr.wkbPolygon, self.epsg)
        v.add(geometry=self.extent_geometry)
        return v

    @property
    def minx(self):
        return self.extent[0]

    @property
    def maxx(self):
        return self.extent[1]

    @property
    def miny(self):
        return self.extent[2]

    @property
    def maxy(self):
        return self.extent[3]

    @property
    def width(self):
        return self.maxx - self.minx

    @property
    def height(self):
        return self.maxy - self.miny

    @property
    def size(self):
        return self.rows * self.columns

    @property
    def shape(self):
        return (self.rows, self.columns)

    @property
    def blocksize(self):
        if not hasattr(self, "_blocksize"):
            self._blocksize = self.band.GetBlockSize()

        return self._blocksize

    @blocksize.setter
    def blocksize(self, value=(256, 256)):
        self._blocksize = list(value)
        self.optimize_blocksize = False

    @property
    def optimized_blocksize(self):
        return optimize_blocksize(self.ds)

    @property
    def tiles(self):
        if not hasattr(self, "_tiles"):
            self._tileopt = self.optimize_blocksize
            self._tiles = tile_list(self.ds, self._tileopt, *self.blocksize)

        if self.optimize_blocksize != self._tileopt:
            self._tileopt = self.optimize_blocksize
            self._tiles = tile_list(self.ds, self._tileopt, *self.blocksize)

        return self._tiles

    @property
    def tile_extents(self):
        return [t.extent_vector.extent for t in self]

    @property
    def settings(self):
        return {
            "spatial_reference": self.spatial_reference,
            "geotransform": self.geotransform,
            "nodata_value": self.nodata_value,
        }

    def load_to_memory(self):
        """loads this datasource to memory, usefull for boosting read array speed"""
        self.change_source(DRIVER_GDAL_MEM.CreateCopy("", self.ds), name=self.name)

    def close(self):
        self.ds.FlushCache()
        del self.ds

    def change_source(self, in_ds, name=None):
        """closes old source and references this object to a new souce"""
        self.ds = None
        self.ds = in_ds
        if name == None:
            self.name = "No name"
        else:
            self.name = name

    def copy(self, shell=False, data_type=None, memory=True):
        if not shell:
            path = loc_dir(mem=memory, add=True)
            if memory:
                driver = DRIVER_GDAL_MEM
            else:
                driver = DRIVER_GDAL_GTIFF

            return Raster(
                driver.CreateCopy(path, self.ds), memory=memory, **self.settings
            )

        if not data_type:
            data_type = self.data_type

        return Raster.from_scratch(
            self.name,
            self.geotransform,
            self.columns,
            self.rows,
            self.nodata_value,
            data_type,
            self.epsg,
            memory,
        )

    def replace_nodata(self, nodata_value=-9999, quiet=False):
        """replaces the in-array nodata values in a tiled form"""
        target = self.copy(shell=True)
        for tile in Progress(self, f"Replace nodata {self.nodata_value}", quiet):
            array = tile.array
            array[np.isnan(array)] = nodata_value
            target.array = array, *tile.location
            tile.close()

        target.nodata_value = nodata_value
        self.change_source(target.ds, name=self.name)

    def change_data_type(self, data_type=float):
        """changes the data type of raster using float, int or gdal dt"""
        if data_type in [float, int]:
            data_type = TYPE_MAPPING_INPUT[data_type]
        out, path = create_source(self.columns, self.rows, data_type)
        target = Raster(out, name=self.name, **self.settings)
        target.array = self.array
        self.change_source(target.ds, name=self.name)

    def nan(self, array):
        """returns an np array which the nodata values are set to np.nan"""
        if np.dtype(array.dtype).char in np.typecodes["AllInteger"]:
            array = array.astype(float)
        array[array == self.nodata_value] = np.nan
        return array

    def resize_band(self, width, height):
        self.band.WriteArray(np.full([height, width], np.nan))

    def get_extent(self):
        """returns the extent of the geotransform"""

        transform = self.geotransform
        maxx = transform[0] + self.columns * transform[1] + self.rows * transform[2]

        miny = transform[3] + self.columns * transform[4] + self.rows * transform[5]

        minx = transform[0]
        maxy = transform[3]
        return minx, maxx, miny, maxy

    def clip(self, vector=None, geometry=None, quiet=False):
        """
        Clips the raster based on an geometry.
        params:
            vector: Clips the raster on a vector.
            geometry: Clips the raster on a geometry.
            quiet: Verbosity.
        return:
            A clipped raster.
        """
        if geometry:
            vector = Vector.from_scratch("", geometry.type, geometry.epsg)
            vector.add(geometry=geometry)

        with tempfile.TemporaryDirectory() as tempdir:

            logger.debug("Clip using the vector path")

            vector.write(tempdir + "/vector.gpkg")
            progress = Progress(total=100, message="Clip", quiet=quiet)

            options = {
                "cutlineDSName": tempdir + "/vector.gpkg",
                "cutlineLayer": vector.layer_name,
                "cropToCutline": True,
                "multithread": True,
                "creationOptions": WRITE_OPTIONS,
                "callback": progress.gdal,
            }

            ds = gdal.Warp(loc_dir(mem=True), self.ds, **options)

        return Raster(ds, name=self.name)

    def clip_bbox(self, minx, maxx, miny, maxy, name=None):
        """clips on a raster on a bbox"""
        return Raster(
            gdal.Translate(
                loc_dir(mem=True, add=True), self.ds, projWin=(minx, maxy, maxx, miny)
            ),
            name=name,
        )

    def reproject(self, epsg, quiet=True):
        progress = Progress(total=100, message="Reproject", quiet=quiet)

        options = {
            "dstSRS": f"EPSG:{epsg}",
            "creationOptions": WRITE_OPTIONS,
            "callback": progress.gdal,
        }

        ds = gdal.Warp(loc_dir(mem=True), self.ds, **options)

        if self.return_raster:
            return Raster(ds)
        else:
            self.change_source(ds, name="reprojected")

    def resample(self, x_res, y_res, algorithm=gdal.GRA_NearestNeighbour, quiet=True):

        progress = Progress(total=100, message="Resample", quiet=quiet)
        options = {
            "xRes": x_res,
            "yRes": y_res,
            "resampleAlg": algorithm,
            "creationOptions": WRITE_OPTIONS,
            "callback": progress.gdal,
        }

        ds = gdal.Translate(loc_dir(mem=True), self.ds, **options)

        if self.return_raster:
            return Raster(ds)
        else:
            self.change_source(ds, name="resampled")

    def smooth(
        self,
        method="gaussian",
        arguments={"sigma": 1, "mode": "nearest"},
        quiet=False,
    ):
        """
        Smooths a raster based on two methods, returns a raster object
        params:
            method: 'gaussian' or 'perona_malik'
            arguments:
                gaussian:
                    - sigma
                    - mode
                    see https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter.html
                perona_malik:
                    - iter
                    - K

        """

        if not HAS_SCIPY:
            raise ImportError("Please install scipy to use this option")

        target = self.copy(shell=True)
        if method == "gaussian":
            if not quiet:
                logger.debug(f"Applying gaussian filter with arguments: {arguments}")

            output_array = scipy.ndimage.gaussian_filter(
                self.array, output=self.np_data_type, **arguments
            )
        elif method == "perona_malik":
            output_array = perona_malik_filter(self.array, **arguments)

        target.array = output_array

        return target

    def push_vector(
        self,
        vector: Vector,
        value=None,
        field=None,
        set_to_nodata=False,
        quiet=False,
    ):
        """
        Pushes a vector this raster

        Params:
            Vector is a vector object
            Value (optional): pushes this value onto the raster
            Field (optional): pushes this field onto the raster
            set_to_nodata (optional): pushes nodata onto the raster

        """

        if self.epsg != vector.epsg:
            vector = vector.reproject(self.epsg)

        if value == self.nodata_value:
            set_to_nodata = True

        target = self.copy(shell=True)
        nodata = self.nodata_value

        # rasterize vector vector array
        rasterized_vector = Raster(
            vector.rasterize(
                rows=self.rows,
                columns=self.columns,
                geotransform=self.geotransform,
                nodata=nodata,
                field=field,
                return_ds=True,
            )
        )

        # make sure that it loops over the same blocks
        rasterized_vector.blocksize = self.blocksize
        self.optimize_blocksize = False

        for raster_tile, vector_tile in zip(self, rasterized_vector):
            array = raster_tile.array
            vector_array = vector_tile.array

            # if no value and field are filled, it is set to nodata?
            if set_to_nodata:
                array[vector_array == 1] = nodata
                target.array = (array, *raster_tile.location)

            elif value:
                array[~np.isnan(vector_array)] = value
                target.array = (array, *raster_tile.location)

            elif field:
                mask = ~np.isnan(vector_array)
                array[mask] = vector_array[mask]
                target.array = (array, *raster_tile.location)

            del array
            del vector_array

        return target

    def push_raster(self, raster, value=None):
        """pushes raster data to object, you can also use a value, rasters should be aligned"""

        if raster.shape != self.shape:
            raise ValueError(
                "Input is not aligned! Use input_raster = input_raster.align(raster)"
            )

        target = self.copy(shell=True)
        self.optimize_blocksize = False
        raster.optimize_blocksize = False

        logger.debug("raster - push_raster - start")
        for target_tile, raster_tile in zip(self, raster):
            target_array = target_tile.array
            raster_mask = raster_tile.mask

            if value:
                target_array[raster_mask] = value
            else:
                raster_array = raster_tile.array
                target_array[raster_mask] = raster_array[raster_mask]
            target.array = (target_array, *raster_tile.location)

        del target_array
        del raster_mask
        return target

    def polygonize(
        self,
        field_name="dn",
        layer_name="poly",
        mask_value=None,
        idx=0,
        no_mask=False,
        quiet=True,
    ):
        """polygonizes a raster
        params:
            field_name: field in vector
            layer_name: layer in vector
            mask_value: a value will be set to nodata
            no_mask: also polygonizes nodata
            quiet: if false it shows progress

        return:
            path: str

        """
        polygonize_type = TYPE_MAPPING[self.band.DataType]
        path = loc_dir(ext="shp", mem=True)
        ds_vector = DRIVER_OGR_SHP.CreateDataSource(path)
        layer_vector = ds_vector.CreateLayer(layer_name, srs=self.spatial_reference)
        raster_field = ogr.FieldDefn(field_name, polygonize_type)
        layer_vector.CreateField(raster_field)

        if mask_value is not None:
            logger.debug("loading raster in memory due to mask value edits")
            self.load_to_memory()
            array = self.array
            array[array == mask_value] = self.nodata_value
            self.array = array
            no_mask = False

        band = self.band
        array = self.array

        if no_mask:
            no_data_band = None
        else:
            no_data_band = band.GetMaskBand()

        band.FlushCache()

        if not quiet:

            progress = Progress(total=100, message="Polygonize", quiet=quiet)
            callback = progress.gdal
        else:
            callback = None

        if polygonize_type == ogr.OFTReal:
            gdal.FPolygonize(
                band,
                no_data_band,
                layer_vector,
                idx,
                [],
                callback=callback,
            )
        else:
            gdal.Polygonize(
                band,
                no_data_band,
                layer_vector,
                idx,
                [],
                callback=callback,
            )

        return path

    def align(self, template=None, fill_value=None, idw=False, nodata_align=True):
        return align(self, template, fill_value, idw, nodata_align)

    def fill_nodata(
        self,
        search_distance=10,
        smooth_iterations=0,
        max_iterations=10,
        array=None,
        quiet=False,
    ):

        """
        Inverse distance weighting - Fills a raster until all nodata is filled.
        For each iterations one pixel if filled
        Params:
            Search Distance: Maximum distance (in pixels) to search out for values to interpolate
                             This is also the amount of nodata pixels it will fill
            Smooth iterations: Number of smoothing iterations to run after the interpolation
            Max Iterations: Limits the amounts of IDW loops


        return a Raster object
        """

        if array is None:
            array = self.array  # nan
        else:
            array = self.nan(array)

        # back to  not np.nan values
        array[np.isnan(array)] = self.nodata_value

        kwargs = self.settings
        kwargs["spatial_reference"] = self.spatial_reference.wkt

        # combine both nans and nodata_values
        fill_mask = array == self.nodata_value
        nan_count = np.count_nonzero(array[fill_mask])

        target = self.copy(shell=True)
        idw_iterations = 0
        while nan_count > 0:
            idw_iterations += 1
            if not quiet:
                logger.debug(f"IDW iteration {idw_iterations}, nan count {nan_count}")

            with Dataset(array, **kwargs) as work:
                work_band = work.GetRasterBand(1)
                mask_band = work_band.GetMaskBand()

                gdal.FillNodata(
                    work_band,
                    mask_band,
                    search_distance,  # search distance
                    smooth_iterations,  # smoothing iterations
                )

            nan_count = np.nansum(array == self.nodata_value)

            if idw_iterations >= max_iterations:
                logger.debug(f"Reached max iterations {max_iterations}")
                break

        if idw_iterations < max_iterations:
            logger.debug("No more nodata left")

        target.array = array
        return target

    def read_point(self, point: ogr.Geometry, nodata_value=-9999):
        """
        Samples points, based on 'raster'.
        input is an ogr geometry
        return a dictionary with sampled values
        """
        # inverse
        p, a, b, q, c, d = self.geotransform
        D = 1 / (a * d - b * c)
        e, f, g, h = d * D, -b * D, -c * D, a * D

        x, y = point.GetX(), point.GetY()
        j = np.int64(e * (x - p) + f * (y - q))
        i = np.int64(g * (x - p) + h * (y - q))

        try:
            array = self.band.ReadAsArray(int(j), int(i), 1, 1).astype(np.float)
            return array[0][0]
        except:
            return nodata_value  # Point outside raster?

    def read_geometry(self, geometry: ogr.Geometry, count=False):
        """
        Return numpy array.

        :param inflate: inflate envelope to grid, to make sure that
            the entire geometry is contained in resulting indices.

        If the bounds fall outside the dataset, the result is padded
        with no data values.

        This ingenuity is from raster-tools
        """

        geotransfrom = self.geotransform
        spatial_reference = self.spatial_reference
        nodata_value = self.nodata_value

        x1, y1, x2, y2 = geotransfrom.get_indices(geometry, inflate=True)

        # overlapping bounds
        w, h = self.columns, self.rows
        p1 = min(w, max(0, x1))
        q1 = min(h, max(0, y1))
        p2 = min(w, max(0, x2))
        q2 = min(h, max(0, y2))

        # result array plus a view for what's actually inside datasets
        array = np.full((y2 - y1, x2 - x1), nodata_value, self.np_data_type)
        view = array[q1 - y1 : q2 - y1, p1 - x1 : p2 - x1]

        kwargs = {
            "xoff": p1,
            "yoff": q1,
            "win_xsize": p2 - p1,
            "win_ysize": q2 - q1,
        }
        data = self.band.ReadAsArray(**kwargs)
        index = data != nodata_value
        view[index] = data[index]

        if 0 in view.shape:
            if count:
                return view, 1
            return view

        # not using self.settings due to speed
        kwargs = {
            "spatial_reference": spatial_reference.wkt,
            "geotransform": geotransfrom.shifted(geometry, inflate=True),
            "nodata_value": nodata_value,
        }

        select_2d = np.zeros(view.shape, dtype="u1")
        with Layer(geometry) as layer:
            with Dataset(select_2d[np.newaxis], **kwargs) as dataset:
                gdal.RasterizeLayer(dataset, [1], layer, burn_values=[1])

        if array.shape != select_2d.shape:
            if count:
                return np.array([0]), 1
            return np.array([0])

        else:
            array_1d = array[select_2d.astype("b1")]

        # determine data or no data
        if array_1d.dtype.kind == "f":
            select_1d = ~np.isclose(array_1d, nodata_value)
        else:
            select_1d = ~np.equal(array_1d, nodata_value)

        if count:
            return array_1d[select_1d], len(array_1d)

        return array_1d[select_1d]

    def read(self, geometry: ogr.Geometry, count=False):

        """
        Return numpy array or a value
        Functions combines an read_geometry and read_point
        params:
            count: If True, returns the amount of selected pixels including nodata.

        """
        # find indices
        if geometry.GetGeometryType() == ogr.wkbPoint:
            return self.read_point(geometry)
        else:
            return self.read_geometry(geometry, count)

    def classify(self, new_classes, quiet=True):
        """
        classifies the raster based on a dictionary.
        e.g.,
            new_classes = {raster_value: [new_value]}
        """
        logger.info("start classifying raster")

        target = self.copy(shell=True, data_type=gdal.GDT_Float32)
        for tile in Progress(self, "Classifying", quiet):
            array = tile.array

            if np.isnan(array).all():
                target.array = array, *tile.location
                continue

            codes = np.unique(array[~np.isnan(array)])

            output_array = np.copy(array)
            for code in codes:
                output_array[array == code] = new_classes[code]

            target.array = output_array, *tile.location

            # # Delete array
            del array
            del output_array

        return target

    def write(self, filename, overwrite=True, write_options=[]):
        filename = str(filename)  # pathlib
        if overwrite:
            if os.path.exists(filename):
                DRIVER_GDAL_GTIFF.Delete(filename)

        WRITE_OPTIONS.extend(write_options)
        ds = DRIVER_GDAL_GTIFF.CreateCopy(filename, self.ds, options=WRITE_OPTIONS)
        ds.GetRasterBand(1).SetNoDataValue(self.nodata_value)
        ds.SetProjection(self.spatial_reference.wkt)
        ds.SetGeoTransform(self.geotransform)
        ds.SetMetadata(self.metadata)
        ds = None

        return


def loc_dir(_dir=None, _type="ras", ext="tif", mem=False, add=True):
    """Creates a string to be used in a datasource"""
    global _loc_ras
    if add:
        _loc_ras = _loc_ras + 1

    if mem:
        _dir = "/vsimem"
    else:
        _dir = tempfile.mkdtemp()

    return _dir + f"/{_type}_{_loc_ras}.{ext}"


def create_source(columns, rows, data_type, memory=True):
    """Takes columns, rows and a data_type to create an in memory gdal datasource"""
    path = loc_dir(mem=memory, add=True)
    if memory:
        driver = DRIVER_GDAL_MEM
    else:
        driver = DRIVER_GDAL_GTIFF

    out = driver.Create(
        path,
        columns,
        rows,
        1,
        data_type,
        options=WRITE_OPTIONS,
    )
    return out, path


def moving_window_generator(
    array, row_min=1, row_max=1, column_min=1, column_max=1, quiet=True
):
    """generates a moving window on top of the class"""
    rows, columns = array.shape
    for row in range(rows):
        for column in range(columns):

            # level 2 create window
            bottom = row - row_min
            top = row + row_max
            left = column - column_min
            right = column + column_max

            # window check
            if bottom < 0:
                bottom = 0
            if top > rows:
                top = rows
            if left < 0:
                left = 0
            if right > columns:
                right = columns
            # extract views
            yield array[bottom : top + 1, left : right + 1], (
                bottom,
                top,
                left,
                right,
            ), (row, column)


def optimize_blocksize(dataset, target_pixels=1000000):
    """Takes a band and creates an optimzed blocksize
    1. Blocksizes are at least orginal band block sizes
    2. New blocksizes are a multitude of original sizes
    3. Block size cannot be larger than the raster size
    4. If target pixels are approximated.

    """

    band = dataset.GetRasterBand(1)
    x_blocksize, y_blocksize = band.GetBlockSize()
    total_block_pixels = x_blocksize * y_blocksize
    total_raster_pixels = dataset.RasterXSize * dataset.RasterYSize

    x_optimized = x_blocksize
    y_optimized = y_blocksize

    while total_block_pixels < target_pixels:
        pre_x_optimized = x_optimized
        pre_y_optimized = y_optimized

        # block gets larger with multitude of blocks
        if x_optimized < y_optimized:
            x_optimized += x_blocksize
        else:
            y_optimized += y_blocksize

        # block cannot be larger than total_raster
        total_block_pixels = x_optimized * y_optimized
        if total_block_pixels > total_raster_pixels:
            x_optimized = pre_x_optimized
            y_optimized = pre_y_optimized
            break

    return x_optimized, y_optimized


def tile_list(dataset: gdal.Dataset, optimize=True, blocksize_x=2048, blocksize_y=2048):
    """Takes a gdal ds and returns a list of coords based on blocks"""
    if optimize:
        x_blocksize, y_blocksize = optimize_blocksize(dataset)
    else:
        x_blocksize = blocksize_x
        y_blocksize = blocksize_y

    x_size = dataset.RasterXSize
    y_size = dataset.RasterYSize

    tiles = []
    for x in range(0, x_size, x_blocksize):

        if x + x_blocksize < x_size:
            cols = x_blocksize
        else:
            cols = x_size - x

        for y in range(0, y_size, y_blocksize):

            if y + y_blocksize < y_size:
                rows = y_blocksize
            else:
                rows = y_size - y
            tiles.append((x, y, cols, rows))
    return tiles


def tile_generator(
    dataset: gdal.Dataset, optimize: bool = True, blocksize_x=2048, blocksize_y=2048
):
    """returns a list of tiles in memory"""
    band = dataset.GetRasterBand(1)
    transform = dataset.GetGeoTransform()
    nodata_value = band.GetNoDataValue()
    projection = dataset.GetProjection()

    for x, y, cols, rows in tile_list(dataset, optimize, blocksize_x, blocksize_y):
        tile_array = band.ReadAsArray(x, y, cols, rows)
        minx = transform[0] + transform[1] * x
        maxy = transform[3] + transform[5] * y
        tile_gt = (minx, transform[1], 0, maxy, 0, transform[5])
        tile_ds, path = create_source(cols, rows, gdal.GDT_Float32)
        tile_ds.SetGeoTransform(tile_gt)
        tile_ds.SetProjection(projection)
        tile_ds.GetRasterBand(1).SetNoDataValue(nodata_value)
        tile_ds.GetRasterBand(1).WriteArray(tile_array)
        tile_array = None
        yield tile_ds, (x, y)


def align(raster, template=None, fill_value=None, idw=False, nodata_align=True):
    """
    Takes a Raster object and returns a aligned raster.
    Assumes the current location of the rasters is correct
    and a similar pixel size
    params:
        template: raster to which the raster is aligned
        fill_value: if nodata needs to be filled, it will use this fill value
        idw: Uses idw to fill nodata needs
        nodata_align: If checked, nodata is aligned

    """

    # template info
    gtt = template.geotransform
    columns = template.columns
    rows = template.rows

    # coordinate difference
    gti = raster.geotransform
    gtt_x, gtt_y = (gtt[0], gtt[3])
    gti_x, gti_y = (gti[0], gti[3])
    x_dif = (gtt_x - gti_x) / gtt[1]
    y_dif = (gtt_y - gti_y) / abs(gtt[5])

    logger.debug("Padding and deleting raster pixels where needed")

    raster_array = raster.array
    raster_array[np.isnan(raster_array)] = raster.nodata_value

    if x_dif > 0:
        raster_array = raster_array[:, abs(int(x_dif)) :]
        # raster is to much to the left, delete values to the left
    elif x_dif < 0:
        # raster is too much to the right, add values on the left side
        raster_array = np.pad(
            raster_array,
            (
                (0, 0),  # pad upside  # pad downside
                (abs(int(x_dif)), 0),  # pad left side  # pad on the right side
            ),
            constant_values=raster.nodata_value,
        )

    col_dif = columns - raster_array.shape[1]
    if col_dif > 0:
        # add columns, width is too small
        raster_array = np.pad(
            raster_array,
            [
                (0, 0),  # pad upside  # pad downside
                (
                    0,
                    abs(int(col_dif)),
                ),  # pad left side  # pad on the right side
            ],
            constant_values=raster.nodata_value,
        )

    elif col_dif < 0:
        # delete columns
        raster_array = raster_array[:, 0 : raster_array.shape[1] - abs(int(col_dif))]

    # works
    if y_dif > 0:
        # raster is to low
        # add rows on the upside if raster is to low
        raster_array = np.pad(
            raster_array,
            [
                (abs(int(y_dif)), 0),  # pad upside  # pad downside
                (0, 0),  # pad left side  # pad on the right side
            ],
            constant_values=raster.nodata_value,
        )
    elif y_dif < 0:
        # raster is too high
        # delete rows from upside if raster is too hgh
        raster_array = raster_array[abs(int(y_dif)) :, :]

    row_dif = rows - raster_array.shape[0]
    if row_dif > 0:
        raster_array = np.pad(
            raster_array,
            [
                (0, abs(int(row_dif))),  # pad upside  # pad downside
                (0, 0),  # pad left side  # pad on the right side
            ],
            constant_values=raster.nodata_value,
        )

    elif row_dif < 0:
        # raster is too large in height
        raster_array = raster_array[0 : raster_array.shape[0] - abs(int(row_dif)) :, :]

    if nodata_align == True:
        logger.debug("Padding nodata values onto the to be aligned array")

        nodata_index = np.isnan(template.array)
        raster_array[nodata_index] = raster.nodata_value

        if np.sum(raster_array == raster.nodata_value) != np.sum(nodata_index):
            if idw == True:
                pass  # done later in the process

            elif not fill_value == None:
                (values, counts) = np.unique(raster_array, return_counts=True)
                fill_value = values[np.argmax(counts)]
                if fill_value == raster.nodata_value:
                    fill_value = values[np.argmax(counts) - 1]

                raster_array[
                    ~nodata_index & (raster_array == raster.nodata_value)
                ] = fill_value

                logger.debug(
                    f"Padding fill value {fill_value} onto the to be aligned array"
                )

    raster_array[raster_array == raster.nodata_value] = np.nan
    output = template.copy(shell=True)
    output.array = raster_array

    if idw == True:
        logger.debug("Using inverse distance weighting to fill up nodata")
        output.fill_nodata()

    output.name = raster.name

    return output


if HAS_SCIPY:

    def perona_malik_filter(data, iter=5, K=5, delta=0.2, method=2, **kwargs):
        """returns a filtered array"""

        dxW = np.array([[0, 0, 0], [1, -1, 0], [0, 0, 0]])
        dxE = np.array([[0, 0, 0], [0, -1, 1], [0, 0, 0]])
        dyN = np.array([[0, 1, 0], [0, -1, 0], [0, 0, 0]])
        dyS = np.array([[0, 0, 0], [0, -1, 0], [0, 1, 0]])

        pm_lambda = np.ones(data.shape) * delta

        k1 = np.ones(data.shape) * (-1.0 / K)
        k2 = np.ones(data.shape) * K * K

        for i in range(iter):

            dI_W = scipy.ndimage.filters.convolve(data, dxW)
            dI_E = scipy.ndimage.filters.convolve(data, dxE)
            dI_N = scipy.ndimage.filters.convolve(data, dyN)
            dI_S = scipy.ndimage.filters.convolve(data, dyS)

            if method == 1:

                cW = np.exp(dI_W * dI_W * k1)
                cE = np.exp(dI_E * dI_E * k1)
                cN = np.exp(dI_N * dI_N * k1)
                cS = np.exp(dI_S * dI_S * k1)

            elif method == 2:

                cW = np.ones(data.shape) / (np.ones(data.shape) + (dI_W * dI_W / k2))
                cE = np.ones(data.shape) / (np.ones(data.shape) + (dI_E * dI_E / k2))
                cN = np.ones(data.shape) / (np.ones(data.shape) + (dI_N * dI_N / k2))
                cS = np.ones(data.shape) / (np.ones(data.shape) + (dI_S * dI_S / k2))

            data = data + (pm_lambda * (cN * dI_N + cS * dI_S + cE * dI_E + cW * dI_W))

        return data
