# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 17:02:45 2021

@author: chris.kerklaan



"""
# First-party imports
import logging
import pathlib

# Third-party imports
import numpy as np
from osgeo import gdal
from osgeo import ogr

# local imports
from .raster import Raster, loc_dir
from .vector import Vector
from threedi_edits.utils.project import Progress
from threedi_edits.utils.project import Classes, Functions

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)


# GLOBALS
logger = logging.getLogger(__name__)
EPSG = 28992
DATA_TYPE = gdal.GDT_Float32
NODATA = -9999
NP_DATA_TYPE = "f4"
PLACEMENT = "stacked"
BLOCKSIZE = None

# DRIVERS
DRIVER_GDAL_GTIFF = gdal.GetDriverByName("GTiff")

# Write options
WRITE_OPTIONS = ["COMPRESS=DEFLATE", "NUM_THREADS=4", "ZLEVEL=9", "TILED=YES"]


class RasterGroup:
    """
    RasterGroup -- info
    RasterGroup is a class which makes handling multiple rasters easier.
    Multiple placement types of rastersgroups are available:
        - stacked (default)
        - adjacent

    Functions include:
            - Aligning all rasters
            - Finding the maximum extent
            - Clipping all rasters by  a vector file

            Rastergroup loops:
                for tiles in group.tiles:
                    tiles[0]# raster 0
                    tiles[1] # raster 1
                    etc.

            Inflated rastrgroup loops:
                inflated_tile_generator = group.inflated_tiles(buffer)
                for tiles, extent in zip(inflated_tile_generator, group.tile_extents):
                    # do something with tiles
                    tiles[0].array = tiles[0].array * tiles[1].array
                    clipped = tiles[0].clip_bbox(*extent)
                    clipped.write("something.tif")

    Params:
        rasters: list of Raster objects
        epsg: epsg format
        nodata: standard nodata for all rasters
        data_type: data_type for all rasters
        np_data_type: numpy variant of the data type

    """

    instance = "gis.rastergroup.RasterGroup"

    def __init__(
        self,
        rasters: list,
        placement=PLACEMENT,
        epsg=EPSG,
        nodata=NODATA,
        data_type=DATA_TYPE,
        np_data_type=NP_DATA_TYPE,
        blocksize=BLOCKSIZE,
    ):
        self.rasters = rasters
        self.names = [raster.name for raster in self.rasters]
        self.nodata = nodata
        self.epsg = epsg
        self.data_type = data_type
        self.np_data_type = np_data_type
        self.placement = placement
        self._blocksize = blocksize

        if len(self.rasters) > 0:
            self.blocksize = rasters[0].blocksize

        for name, raster in zip(self.names, self.rasters):
            setattr(self, name, raster)

    @classmethod
    def from_paths(
        cls,
        paths,
        epsg=EPSG,
        nodata=NODATA,
        data_type=DATA_TYPE,
        np_data_type=NP_DATA_TYPE,
    ):
        return cls([Raster(p) for p in paths])

    def __iter__(self):
        for raster in self.rasters:
            yield raster

    def __len__(self):
        return sum([len(raster) for raster in self])

    def __getitem__(self, name):
        if type(name) == "str":
            return self.rasters[self.names.index(name)]
        elif type(name) == int:
            return self.rasters[name]

    def __setitem__(self, name, value):
        """without alignemnt"""
        if name in self.names:
            self.rasters[self.names.index(name)] = value
        else:
            self.rasters.append(value)
            self.names.append(name)

        setattr(self, name, value)

    def __repr__(self):
        data = {"aligned": self.aligned}
        return f"({self.instance}) {self.names}: {data}"

    @classmethod
    def from_folder(cls, path):
        rasters = [r for r in pathlib.Path(path).rglob("*.tif")]
        return cls([Raster(str(r)) for r in rasters])

    @property
    def epsg(self):
        if not hasattr(self, "_epsg"):
            self._epsg = self[0].epsg
        return self._epsg

    @epsg.setter
    def epsg(self, value):
        self._epsg = value

    @property
    def count(self):
        return len(self.rasters)

    @property
    def aligned(self):
        self.check_alignment(False)
        if len(self._alignment["to_be_aligned"]) == 0:
            return True
        else:
            return False

    @property
    def aligned_nodata(self):
        self.check_alignment(True)
        if len(self._alignment["to_be_aligned"]) == 0:
            return True
        else:
            return False

    @property
    def tile_extents(self):
        return self.rasters[0].extent

    @property
    def extents(self):
        if not hasattr(self, "_extent"):
            _extent = Vector.from_scratch("extent", ogr.wkbPolygon, self.epsg)
            _extent.add_field("name", str)
            _extent.add_field("index", int)
            for i, raster in enumerate(self):
                _extent.add(
                    geometry=raster.extent_geometry,
                    items={"name": raster.name, "index": i},
                )
            self._extent = _extent
        return self._extent

    @property
    def stacked(self):
        return self.placement == "stacked"

    @stacked.setter
    def stacked(self, value):
        assert value in [True, False]
        if value == True:
            self.placement = "stacked"
        else:
            self.adjacent = True

    @property
    def adjacent(self):
        return self.placement == "adjacent"

    @adjacent.setter
    def adjacent(self, value):
        assert value in [True, False]
        if value == True:
            self.placement = "adjacent"
        else:
            self.stacked = True

    @property
    def blocksize(self):
        return self._blocksize

    @blocksize.setter
    def blocksize(self, blocksize):
        self._blocksize = blocksize
        for r in self:
            r.blocksize = blocksize
            r.optimize_blocksize = False

    def read(self, geometry):
        rasters = self.extents.spatial_filter(geometry)
        if len(rasters) == 0:
            return [np.nan]
        else:
            data = []
            for raster in rasters:
                data.append(self.rasters[raster["index"]].read(geometry))
            return data

    def tiles_inflated(self, buffer_size, original_extent):
        """
        Returns a generator based on the placement of the rastergroup.
        Returns multiple raster tiles if stacked.
        Returns a single raster tile if adjacent.
        """

        if self.stacked:
            if not self.aligned:
                raise ValueError("Not aligned, please align first")

            generators = [r.tiles_inflated(buffer_size) for r in self]

            for inflated_tiles in zip(*generators):
                yield inflated_tiles

        # TODO tiles cannot be inflated in this way.
        if self.adjacent:
            for raster in self:
                for tile in raster.tiles_inflated(buffer_size, original_extent):
                    yield tile

    def tiles(self):
        """
        Returns a generator based on the placement of the rastergroup.
        If stacked, we return mulitple raster tiles.
        If adjacent, we return a single raster tile.
        """

        if self.stacked:
            if not self.aligned:
                raise ValueError("Not aligned, please align first")

            for raster_tiles in zip(*self.rasters):
                yield raster_tiles

        if self.adjacent:
            for raster in self.rasters:
                for tile in raster:
                    yield tile

    def check_alignment(self, count_nodata=True):
        logger.info("Checking alignment of the rasters")
        self._alignment = check_alignment(self.rasters, count_nodata)
        return self._alignment

    def get_maximum_extent(self):
        if not hasattr(self, "maximum_extent"):
            self.maximum_extent = max_template(self.rasters)

    def align(self, index=0, use_maximum_extent=False, nodata=False):
        """
        This function aligns all rasters in the group
        params:
            index: This index is used for the raster on which all other
            rasters are aligned.
            use_maximum_extent: Calculates a maximum extent of all rasters
            and aligns all rasters to that maximum extent

        """
        logger.info("Aligning all rasters")

        if use_maximum_extent:
            self.get_maximum_extent()
            align_raster = self.maximum_extent
            rasters = self.rasters
        else:
            align_raster = self.rasters[index]
            rasters = [r for i, r in enumerate(self.rasters) if i != index]

        aligned_rasters = [r.align(align_raster, nodata_align=nodata) for r in rasters]

        if not use_maximum_extent:
            aligned_rasters.insert(index, align_raster)

        self.rasters = aligned_rasters

    def merge(self):
        """
        Merges the rasters based on the maximum extent
        loops over every tile and adds data on top of the nodata of the other
        """
        logger.info("start merging")

        # they all must be align, have the same extent an the same nodata value
        self.get_maximum_extent()
        self.align(use_maximum_extent=True)
        for raster in self.rasters:
            raster.optimize_blocksize = False

        pbar = Progress(len(raster), "Merging")

        target = self.rasters[0].copy()
        for tiles in zip(*tuple(self.rasters)):
            tile_array = tiles[0].array

            for tile in tiles[1:]:
                mask = np.isnan(tile_array)
                tile_array[mask] = tile.array[mask]

            pbar.update()
            target.array = tile_array, *tile.location

        return target

    def vrt(self):
        """
        Creates one raster out of all adjacent rasters.
        Returns a vrt of all rasters.


        """
        raster_list = [r.ds for r in self.rasters]
        path = loc_dir(mem=True).replace(".tif", ".vrt")
        vrt_options = gdal.BuildVRTOptions(resampleAlg="cubic")
        ds = gdal.BuildVRT(path, raster_list, options=vrt_options)
        return Raster(ds)

    def insert(self, raster: Raster):
        logger.info(f"Inserting a {raster.name} into the group")
        align_raster = self.rasters[0]
        raster.align(align_raster, nodata_align=False)
        self.rasters.append(raster)
        self.names.append(raster.name)

    def clip(self, vector):
        """clips all available rasters to the vector file"""
        for raster in self:
            logger.debug(f"Clipping {raster.name}")
            clipped = raster.clip(vector)
            self[raster.name] = clipped

    def write(self, filename):
        """writes it to one dataset, with different layers"""

        ds = DRIVER_GDAL_GTIFF.Create(filename, options=WRITE_OPTIONS)
        ds.GetRasterBand(1).SetNoDataValue(self.nodata_value)
        ds.SetProjection(self.spatial_reference.wkt)
        ds.SetGeoTransform(self.geotransform)
        ds.SetMetadata(self.metadata)


def check_alignment(raster_list, count_nodata=True):
    """checks data/nodata alignemnt and extent alignment"""
    dem = raster_list[0]

    output = {
        "extent": {},
        "counts": {},
        "location": {},
        "errors": [],
        "to_be_aligned": [],
    }
    for raster in raster_list:
        # extent
        if (raster.columns, raster.rows) != (dem.columns, dem.rows):
            msg = f"{raster.name} has unqual columns and rows "
            logger.debug(msg)
            output["errors"].append(("extent", msg))
            if not raster.name in output["to_be_aligned"]:
                output["to_be_aligned"].append(raster.name)

        output["extent"][raster.name] = {
            "rows": raster.rows,
            "columns": raster.columns,
        }
        if count_nodata:
            # data/nodata
            output["counts"][raster.name] = count_data_nodata(raster)

        # location
        if raster.geotransform != dem.geotransform:
            msg = f"{raster.name} has not a similar geotransform as the first raster"
            logger.debug(msg)
            output["errors"].append(("location", msg))
            if not raster.name in output["to_be_aligned"]:
                output["to_be_aligned"].append(raster.name)

        output["location"][raster.name] = raster.geotransform

    for key, values in output["counts"].items():
        if values != output["counts"][dem.name]:
            msg = f"{key} pixel data/nodata count not equal"
            logger.debug(msg)
            output["errors"].append(("counts", msg))
            if not key in output["to_be_aligned"]:
                output["to_be_aligned"].append(key)

    if len(output["errors"]) == 0:
        logger.debug("RasterGroup - Check alignment found no problems")

    return output


def count_data_nodata(raster):
    """input array has np.nan as nodata"""
    count_data = 0
    count_nodata = 0
    raster.optimize_blocksize = False
    for data in raster:
        arr = data.array
        total_size = arr.size
        add_cnt_nodata = np.count_nonzero(np.isnan(arr))
        add_cnt_data = total_size - add_cnt_nodata
        count_nodata += add_cnt_nodata
        count_data += add_cnt_data
        del arr
    return count_data, count_nodata


def max_template(rasters=[Raster, Raster], resolution=None):
    """Takes a list of Raster objects and returns the largest template
    Note that the resolution of the first raster is taken, if not overwritten
    by the resolution param.
    """

    x_min = 1e31
    x_max = -1e31
    y_min = 1e31
    y_max = -1e31

    if resolution == None:
        resolution = rasters[0].resolution["width"]

    for raster in rasters:

        gtt = raster.geotransform
        # xmin
        if gtt[0] < x_min:
            x_min = gtt[0]

        # xmax
        if gtt[0] + (raster.rows * gtt[1]) > x_max:
            x_max = gtt[0] + (raster.rows * gtt[1])

        # ymax
        if gtt[3] > y_max:
            y_max = gtt[3]

        # ymin
        if gtt[3] + (raster.columns * gtt[5]) < y_min:
            y_min = gtt[3] + (raster.columns * gtt[5])

    rows = (y_max - y_min) / resolution
    columns = (x_max - x_min) / resolution
    template = Raster.from_array(np.zeros((int(columns), int(rows))))
    template.geotransform = (x_min, resolution, 0.0, y_max, 0.0, -resolution)
    template.nodata_value = -9999
    template.spatial_reference = raster.epsg

    return template
