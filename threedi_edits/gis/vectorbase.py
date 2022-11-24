# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:44:42 2022

@author: chris.kerklaan
"""

# First-party imports
import os
import logging
import operator
import pathlib
import tempfile

# Third party imports
import shapely
from cached_property import cached_property
from osgeo import osr, ogr, _ogr, gdal

# Local imports
from threedi_edits.gis.geometry import (
    TRANSLATION as GEOMETRY_TRANSLATE,
    POLYGON,
    to_ogr as to_ogr_geometry,
    shapely_geometry,
    to_shapely as to_shapely_geometry,
)

from threedi_edits.gis.feature import Feature, to_ogr as to_ogr_feature
from threedi_edits.gis.fields import Fields, GET_FIELD
from threedi_edits.gis.polygon import Polygon
from threedi_edits.gis.spatial_reference import SpatialReference
from threedi_edits.gis.utils import DRIVERS
from threedi_edits.utils.project import Classes, Functions

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)


# Global DRIVERS
DRIVER_GDAL_MEM = gdal.GetDriverByName("MEM")
DRIVER_OGR_MEM = ogr.GetDriverByName("Memory")
DRIVER_OGR_GPKG = ogr.GetDriverByName("GPKG")
DRIVER_OGR_SHAPEFILE = ogr.GetDriverByName("ESRI Shapefile")
DRIVERS_SQL = ["PostgreSQL", "GPKG"]

# Memory counter
_mem_num = 0

# Logger
logger = logging.getLogger(__name__)


class VectorBase:
    """
    Basewrapper of an ogr layer. This provides the most optimalized
    functions to acces an ogr layer and its data.
    """

    instance = "gis.vector.VectorBase"

    def __init__(
        self,
        path: str = None,
        layer_name: str = None,
        write: int = 1,
        view: bool = False,
        name: str = None,
        epsg: int = 28992,
        datasource: ogr.DataSource = None,
    ):
        path = str(path)  # pathlib input
        self.layer_name = layer_name
        self.view = view
        self.standard_epsg = epsg
        self.memory = "vsimem" in path

        if path is None and datasource is None:
            raise FileNotFoundError("Bath path and datasource are not given.")

        if path is not None and datasource is None:
            self.path = pathlib.Path(path)
            if not self.path.exists() and not self.memory:
                raise FileNotFoundError("Vector not found")
            self.ds = ogr.Open(str(path), write)

        if datasource is not None:
            self.ds = datasource
            self.path = pathlib.Path("")

        if self.path.suffix == ".shp":
            layer_name = self.path.stem

        if layer_name:
            self.layer_name = layer_name
            self.layer = self.ds.GetLayer(layer_name)
            if self.layer is None:
                raise FileNotFoundError(f"{layer_name} does not exist.")
        else:
            try:
                self.layer = self.ds.GetLayer()
                self.layer_name = self.layer.GetName()
            except TypeError:
                raise FileNotFoundError("No layer name, please define")

        if not name:
            name = self.layer_name

        self.name = name
        self.table  # initialize

    @classmethod
    def from_pg(cls, host, port, user, password, dbname, layer_name, write=0, **kwargs):
        """
        returns a class from a postgres database
        takes a {"host":"x", "port":"x", "password":"x",
                 "user":"x", "host":"x"}
        """
        ogr_pg = ("PG:host={} port={} user='{}'" "password='{}' dbname='{}'").format(
            host, port, user, password, dbname
        )

        return cls(
            path=None,
            datasource=ogr.Open(ogr_pg, write),
            layer_name=layer_name,
            **kwargs,
        )

    @classmethod
    def from_ds(cls, datasource: ogr.DataSource = None, **kwargs):
        """returns a class from datasource input"""
        return cls(path=None, datasource=datasource, **kwargs)

    @classmethod
    def from_mem_path(cls, path: str = None, **kwargs):
        """returns a class from a memory path"""
        return cls(path=path, **kwargs)

    @classmethod
    def from_scratch(cls, layer_name: str, geometry_type: int, epsg: int, **kwargs):
        """Creates a layer from scratch
        Params:
            layer_name: name of the layer (str)
            geometry_type: 1 = point, 2 = linestring, 3 = polygon (int)
            epsg: spatial reference code (int)
        """
        return cls(
            path=None,
            datasource=mem_layer(layer_name, geometry_type=geometry_type, epsg=epsg)[0],
        )

    def get_feature(self, fid=None, feature=None):
        """helper for retrieving features"""
        if fid is not None:
            return Feature(
                _ogr.Layer_GetFeature(self.layer, fid), self.layer, self.fields
            )

        if feature:
            return Feature(feature, self.layer, self.fields)

    def __iter__(self):
        self.reset()
        while True:
            feature = _ogr.Layer_GetNextFeature(self.layer)
            if not feature:
                break
            yield self.get_feature(feature=feature)

    def __getitem__(self, i):
        return self.get_feature(i)

    def __len__(self):
        return self.count

    def __setitem__(self, key, value):
        if key == "feature":
            self.set_feature(value)

    def __repr__(self):
        settings = self.settings
        return f"({self.instance}) {self.name} ({len(self)}): {settings}, Type: {self.geometry_name}"

    def __next__(self):
        return self.get_feature(feature=self.layer.GetNextFeature())

    def __in__(self, fid):
        if hasattr(fid, "fid"):
            fid = fid.fid
        return fid in self.fids

    @property
    def count(self):
        return _ogr.Layer_GetFeatureCount(self.layer)

    @property
    def filepath(self):
        return str(self.path)

    @property
    def filename(self):
        return self.path.with_suffix("")

    @property
    def layers(self):
        return [layer.GetName() for layer in self.ds]

    @property
    def extent(self):
        return self.layer.GetExtent()

    @cached_property
    def fields(self):
        return Fields(self.layer)

    @cached_property
    def keys(self):
        return self.fields.keys

    @cached_property
    def field_indices(self):
        return self.fields.indices

    @property
    def fids(self):
        return self.table["fid"]

    @property
    def layer_defn(self):
        return self.layer.GetLayerDefn()

    @property
    def geometry(self):
        return self[0].geometry

    @property
    def geometry_type(self):
        return self.layer.GetGeomType()

    @property
    def geometry_name(self):
        return ogr.GeometryTypeToName(self.geometry_type)

    @property
    def fid_column(self):
        column = self.layer.GetFIDColumn()
        if column == "":
            return "id"
        return self.layer.GetFIDColumn()

    @property
    def driver(self):
        return self.ds.GetDriver().name

    @property
    def spatial_reference(self):
        sr = self.layer.GetSpatialRef()
        if not sr:
            return SpatialReference.from_epsg(self.standard_epsg)
        else:
            return SpatialReference.from_sr(sr)

    @property
    def epsg(self):
        return self.spatial_reference.epsg

    @property
    def extent_geometry(self):
        x1, x2, y1, y2 = self.layer.GetExtent()
        return Polygon.from_wkt(POLYGON.format(x1=x1, x2=x2, y1=y1, y2=y2))

    @property
    def extent_area(self):
        return self.extent_geometry.GetArea()

    @property
    def settings(self):
        return {
            "view": self.view,
            "layer_name": self.layer_name,
            "epsg": self.epsg,
            "name": self.name,
        }

    @property
    def as_gdal_dataset(self):
        """gdal and ogr not always compatible, we have to use openex"""
        mem_ds = mem_path() + ".shp"
        self.write(mem_ds, index=False)
        return gdal.OpenEx(mem_ds)

    @cached_property
    def table(self):
        return load_table(
            self.layer, self.fields.names, self.fields.types, self.has_geometry
        )

    @property
    def has_ogr_index(self):
        return self.layer.TestCapability("FastSpatialFilter")

    @property
    def has_geometry(self):
        return self.geometry_type != 100

    @property
    def geometries(self):
        return self.table["geometry"]

    def reset(self):
        _ogr.Layer_ResetReading(self.layer)

    def first(self):
        return self.get_feature(self.fids[0])

    def add_field(self, name, ogr_type):
        self.fields.add(name, ogr_type)
        self.add_table_field(name)
        self.delete_cached_fields()

    def delete_field(self, *field_names):
        for field_name in field_names:
            self.fields.delete(field_name)
            self.delete_table_field(field_name)

        self.delete_cached_fields()

    def _add(
        self,
        feature=None,
        geometry=None,
        fid=-1,
        fid_latest=False,
        **items,
    ):

        # unpack feature
        if isinstance(feature, Feature):
            geometry = feature.geometry
            items = feature.items

        # pick up items from kwargs
        attr = {}
        if "items" in items:
            attr = items["items"]
            del items["items"]
        attr.update(items)
        items = attr

        # account for none attributes
        if len(items) != 0:
            irrelevant = [k for k in items.keys() if k not in self.keys]
            for key in irrelevant:
                del items[key]

        # convert geometry if needed
        if geometry is not None:
            geometry = to_ogr_geometry(geometry)

        if fid_latest:
            fid = self._get_last_feature_id() + 1

        self.layer, fid = add_feature(
            self.layer,
            self.layer_defn,
            geometry,
            items,
            fid,
        )

        if self.has_geometry:
            items["geometry"] = shapely_geometry(geometry)

        items["fid"] = fid
        self.add_table_feature(items)

        return fid

    def delete(self, feature):
        """deletes by using feature ids
        Input can be a single feature id or of type Feature class
        """

        if type(feature) == int:
            feature = self[feature]

        feature = to_ogr_feature(feature)

        if self.driver == "PostgreSQL":
            self.layer.ResetReading()
            self.ds.StartTransaction()

        self.delete_table_feature(feature)
        self.layer.DeleteFeature(feature.GetFID())

        if self.driver == "PostgreSQL":
            self.ds.CommitTransaction()

    def delete_all(self):
        self.reset()
        for feature in self:
            self.delete(feature)

        self.layer.SyncToDisk()
        self.ds.FlushCache()

        self.layer = None
        self.layer = self.ds.GetLayer(self.layer_name)

    def set_feature(self, feature):
        self.layer.SetFeature(to_ogr_feature(feature))

    def _copy(
        self,
        geometry_type=None,
        shell=False,
        fields=True,
        driver="GPKG",
        spatial_reference=None,
    ) -> (ogr.DataSource, str):
        """In memory ogr object, with features
        params:
            geometry_type: change the geometry type
            shell: Does not add the features
            fields: Set to True is fields have to be included
        """

        if type(geometry_type) == str:
            geometry_type = GEOMETRY_TRANSLATE[geometry_type]

        if not geometry_type:
            geometry_type = self.geometry_type

        ds, layer, path = mem_layer(
            self.layer_name,
            self.layer,
            geometry_type=geometry_type,
            shell=shell,
            fields=fields,
            sr=self.spatial_reference,
            return_path=True,
            driver=driver,
        )
        self.path = path
        return ds, path

    def write(self, path, overwrite=False, driver_name=None):
        """
        Write is done by creating a new file and copying features
        params:
            path: output path of the vector.
            overwrite: Deletes the path is exists.
            driver_name: Short of the ogr driver.
        """
        path = str(path)  # pathlib

        if driver_name is not None:
            driver = ogr.GetDriverByName(driver_name)
        else:
            driver = self.identify_driver(path)
        driver_name = driver.name

        if overwrite:
            if os.path.exists(path):
                driver.DeleteDataSource(path)

        # create output file
        out_ds = driver.CreateDataSource(path)
        out_layer = out_ds.CreateLayer(
            self.layer_name, self.spatial_reference, self.geometry_type
        )

        # add fields
        for i in range(self.layer_defn.GetFieldCount()):
            out_layer.CreateField(self.layer_defn.GetFieldDefn(i))

        # use out layer definition
        layer_defn = out_layer.GetLayerDefn()

        # if uploading to sql based drivers, we have to commit once in a while
        sql = driver_name in DRIVERS_SQL
        if sql:
            out_layer.StartTransaction()

        for i, feature in enumerate(self):
            add_feature(
                out_layer,
                layer_defn,
                to_ogr_geometry(feature.geometry),
                feature.items,
            )
            if sql and (i % 1500 == 0):  # commit once every 1500 features
                out_layer.CommitTransaction()
                out_layer.StartTransaction()

        if sql:
            out_layer.CommitTransaction()

        out_layer = None
        out_ds.Destroy()

    def close(self):
        self.layer = None
        self.ds.Destroy()

    """tables and index"""

    def add_table_feature(self, items):
        for key in self.table:
            if key not in items:
                value = None
            else:
                value = items[key]
            self.__dict__["table"][key].append(value)

    def add_table_field(self, field):
        self.__dict__["table"][field] = [None] * len(self.table["fid"])

    def delete_table_field(self, field):
        del self.__dict__["table"][field]

    def delete_table_feature(self, feature: ogr.Feature):
        """delete a feature in the table"""
        idx = self.table["fid"].index(feature.GetFID())
        for key in self.table:
            del self.__dict__["table"][key][idx]

    def delete_cached_fields(self):
        if "fields" in self.__dict__:
            del self.__dict__["fields"]
        if "keys" in self.__dict__:
            del self.__dict__["keys"]

    """spatial and table filters"""

    def set_spatial_filter(self, geometry):
        self.layer.SetSpatialFilter(geometry)
        self.info(self.layer)

    def _sfilter(self, geometry, method="intersects") -> list:
        """
        Spatial filter of all that is present within the predicates of pygeos.
        https://pygeos.readthedocs.io/en/latest/predicates.html

        Parameters
        ----------
        geometry : shapely.BaseGeometry
            Shapely base geometry.
        method : string
            Intersection method based on the link above.
            The default is "intersects".

        Returns
        -------
        list of ids

        Speed-test with 450.000 features:
            ~32 seconds using this method.
            ~120 with the old method.
            New method is about 4 times faster.

        """
        if not hasattr(shapely.predicates, method):
            raise ValueError(
                f"""
                             Method not in shapely.predicates,
                             please choose from {dir(shapely.predicates)}
                             """
            )
        geometry = to_shapely_geometry(geometry)
        predicate = getattr(shapely.predicates, method)
        data = predicate(self.geometries, geometry).nonzero()
        if len(data[0]) > 0:
            output = operator.itemgetter(*data[0])(self.fids)
            if type(output) not in [list, tuple]:
                return [output]
            return output

        return []

    def _subset(self, *fids) -> str:
        """
        Returns a subset of the vector based on the fids
        fastest way to retrieve fids is via attribute filter and some sql.

        Also when there are more fids that +- 4000 it'll sometimes return an
        error.
        Maybe we should look a settings:
            https://courses.spatialthoughts.com/gdal-tools.html#tips-for-improving-performance
        """

        if len(fids) == 0:
            return self._copy(shell=True)[1]

        if type(fids[0]) == list:
            fids = tuple(fids[0])

        if len(fids) == 1 and type(fids[0]) == int:
            self.layer.SetAttributeFilter(f"FID IN ({fids[0]})")
            return self._copy()[1]

        max_size = 4000
        sets = []
        for i in range(0, len(fids), max_size):
            self.layer.SetAttributeFilter(f"FID IN {fids[i:i + max_size]}")
            sets.append(self._copy()[0])
        output = merge(sets, self.name)

        self.layer.SetAttributeFilter(None)
        return output

    def _filter(self, **filtering) -> str:
        """Filters based  on a dict, returns a Vector"""
        loop_fids = self.table["fid"]
        for field_key, field_value in filtering.items():
            loop_fids = lookup(loop_fids, self.table[field_key], field_value)
        return self._subset(loop_fids)

    def _get_fids(self):
        SQL = f"SELECT fid FROM {self.layer_name}"
        result = self.ds.ExecuteSQL(SQL)
        return [f.GetFID() for f in result]

    def identify_driver(self, path):
        name = DRIVERS[os.path.splitext(path)[1]]
        return ogr.GetDriverByName(name)


def load_table(
    layer: ogr.Layer, field_names: list, field_types: list, has_geometry: bool
) -> dict:
    """
    Retrieves all data from a layer in the quickest way possible.
    We circle around ogr and go directly to the _ogr module.
    This way we are directly accesing the C++ functions.

    Speed-test with 450.000 features:
        ~23 secondds using this method.
        ~163 with the old method.
        New method is about 6 times faster.
        Loading speed increases with more features.

    """

    names = field_names
    if has_geometry:
        extra_fields = ["fid", "geometry"]
    else:
        extra_fields = ["fid"]

    data = {field: [] for field in names + extra_fields}
    indices_types = list(zip(list(range(len(names))), field_types))

    _ogr.Layer_ResetReading(layer)
    feature = _ogr.Layer_GetNextFeature(layer)
    while feature:
        for i, itype in indices_types:
            data[names[i]].append(GET_FIELD[itype](feature, i))

        data["fid"].append(_ogr.Feature_GetFID(feature))
        if has_geometry:
            geometry = _ogr.Feature_GetGeometryRef(feature)

            if geometry is not None:
                data["geometry"].append(shapely_geometry(geometry))
        feature = _ogr.Layer_GetNextFeature(layer)

    return data


def lookup(lookup_table: list, table: list, value: int) -> list:
    """
    Fastest way to lookup a value by using an index in a unordered list
    Looping over a list only once.
    returns a list of matching values.

    A faster way should be with np.searchsorted, however,
    we also have string data.

    """

    values = []

    if value not in table:
        return values

    index = 0

    # last index of value in table
    end = len(table) - operator.indexOf(reversed(table), value) - 1

    # when the first value is the only one
    if end == 0:
        return [lookup_table[0]]

    while index <= end:
        index = table.index(value, index)
        values.append(lookup_table[index])
        index += 1

    return values


def mem_path(extension: str = None) -> str:
    """Generates a new vsimem path with or without an extension"""
    global _mem_num
    location = f"/vsimem/mem{_mem_num}"
    _mem_num += 1
    if extension is not None:
        location = location + f".{extension}"
    else:
        location = location + ".vrt"

    return location


def mem_ds(driver_name: str = "Memory") -> (ogr.DataSource, str):
    """creates a new ogr datasource base on a driver name"""
    driver = ogr.GetDriverByName(driver_name)
    meta = driver.GetMetadata()
    extension = meta.get("DMD_EXTENSION")
    location = mem_path(extension)
    ds = driver.CreateDataSource(location)
    return ds, location


def mem_layer(
    layer_name: str,
    layer: ogr.Layer = None,
    geometry_type: int = None,
    epsg: int = None,
    sr: osr.SpatialReference = None,
    shell: bool = True,
    fields: bool = True,
    options: list = ["SPATIAL_INDEX=YES"],
    return_path: bool = False,
    driver="GPKG",
) -> (ogr.DataSource, ogr.Layer):
    """
    Standard creates a 'shell' of another layer,
    To create a new memory layer,
    or to create a memory copy of a layer
    params:
        epsg/sr: one can choose between epsg/sr
        shell: if True, returns a shell with no  feature
        field: if False, returns a shell without the fields
    returns:
        ds : ogr.Datasource
    """
    ds, path = mem_ds(driver)

    if not shell and layer is not None:
        new_layer = ds.CopyLayer(layer, layer_name)
        if return_path:
            return ds, new_layer, path
        else:
            return ds, new_layer

    if epsg:
        sr = SpatialReference.from_epsg(epsg)
    elif sr == None:
        sr = layer.GetSpatialRef()

    if not geometry_type and layer:
        geometry_type = layer.GetGeomType()
    new_layer = ds.CreateLayer(layer_name, sr, geometry_type, options)

    if fields and layer:
        layer_defn = layer.GetLayerDefn()
        for i in range(layer_defn.GetFieldCount()):
            new_layer.CreateField(layer_defn.GetFieldDefn(i))

    if return_path:
        return ds, new_layer, path
    else:
        return ds, new_layer


def add_feature(
    layer: ogr.Layer,
    feature_defn: ogr.FeatureDefn,
    geometry: ogr.Geometry,
    attributes: dict,
    fid: int = -1,
) -> (ogr.Layer, int):
    """
    Append geometry and attributes as new feature, starting at 0

    for fid=-1 there are three options:
        1. fid is similar to the current count: exisiting fids = [1,4,10], new fid  = 3
        2. fid is 'next number in line': exisiting fids = [1,10], new fid  = 2
        3. fid continues at the latest number: exisiting fids = [1,10], new fid = 11

    """
    feature = ogr.Feature(feature_defn)
    _ogr.Feature_SetFID(feature, fid)

    if geometry is not None:
        _ogr.Feature_SetGeometry(feature, geometry)

    for key, value in attributes.items():
        if value == None:
            _ogr.Feature_SetFieldNull(feature, key)
        else:
            try:
                _ogr.Feature_SetField(feature, key, value)
            except Exception as e:
                raise ValueError("error:", e, "key", key, "value", value)

    try:
        _ogr.Layer_CreateFeature(layer, feature)
    except RuntimeError as e:
        raise RuntimeError("error:", e, "geom:", geometry, "attributes:", attributes)
    finally:
        fid = _ogr.Feature_GetFID(feature)
        _ogr.delete_Feature(feature)

    return layer, fid


def merge(
    datasources: [ogr.DataSource, ogr.DataSource], name="merged"
) -> ogr.DataSource:
    """
    Returns a path for the merged datasources.
    This seems to be the fastest way to merge different vector formats.
    Based on ogrmerge.py and https://gdal.org/drivers/vector/vrt.html
    it uses a vrt combined with a temporary directory to merge the vectors.

    TODO:
        We should lose the temp directory. Maybe we can try this with vsimem.

    """
    with tempfile.TemporaryDirectory() as folder:
        vrt_filename = folder + "/temp.vrt"
        f = gdal.VSIFOpenL(vrt_filename, "wb")
        writer = XMLWriter(f)
        writer.open_element("OGRVRTDataSource")
        gpkg_driver = ogr.GetDriverByName("GPKG")
        writer.open_element("OGRVRTUnionLayer", attrs={"name": name})

        for i, src_ds in enumerate(datasources):
            # first we copy it to the tempfolder
            gpkg_driver.CopyDataSource(src_ds, folder + f"/{i}.gpkg")
            layer = src_ds.GetLayer()
            writer.open_element("OGRVRTLayer", attrs={"name": f"{i}"})
            writer.write_element_value("SrcDataSource", folder + f"/{i}.gpkg")
            writer.write_element_value("SrcLayer", layer.GetName())
            writer.close_element("OGRVRTLayer")

        writer.write_element_value("PreserveSrcFID", "ON")
        writer.close_element("OGRVRTUnionLayer")
        writer.close_element("OGRVRTDataSource")
        gdal.VSIFCloseL(f)

        target = mem_path("gpkg")
        gdal.VectorTranslate(target, vrt_filename, options=["-preserve_fid"])
        gdal.Unlink(vrt_filename)

    return target


# from osgeo utils
def _VSIFPrintfL(f, s):
    gdal.VSIFWriteL(s, 1, len(s), f)


def _Esc(x):
    if gdal.__version__ == "3.5.2":
        return gdal.EscapeString(x, gdal.CPLES_XML).decode("UTF-8")
    else:
        return gdal.EscapeString(x, gdal.CPLES_XML)


class XMLWriter:
    def __init__(self, f):
        self.f = f
        self.inc = 0
        self.elements = []

    def _indent(self):
        return "  " * self.inc

    def open_element(self, name, attrs=None):
        xml_attrs = ""
        if attrs is not None:
            for key in attrs:
                xml_attrs = xml_attrs + ' %s="%s"' % (
                    key,
                    _Esc(attrs[key].encode("utf-8")),
                )
        x = "%s<%s%s>\n" % (self._indent(), name, xml_attrs)
        x = x.encode("utf-8")
        _VSIFPrintfL(self.f, x)
        self.inc = self.inc + 1
        self.elements.append(name)

    def write_element_value(self, name, value, attrs=None):
        xml_attrs = ""
        if attrs is not None:
            for key in attrs:
                xml_attrs = xml_attrs + ' %s="%s"' % (
                    key,
                    _Esc(attrs[key].encode("utf-8")),
                )
        x = "%s<%s%s>%s</%s>\n" % (
            self._indent(),
            name,
            xml_attrs,
            _Esc(value.encode("utf-8")),
            name,
        )
        x = x.encode("utf-8")
        _VSIFPrintfL(self.f, x)

    def close_element(self, closing_name=None):
        self.inc = self.inc - 1
        name = self.elements[-1]
        if closing_name is not None:
            assert name == closing_name
        self.elements = self.elements[0:-1]
        _VSIFPrintfL(self.f, "%s</%s>\n" % (self._indent(), name))
