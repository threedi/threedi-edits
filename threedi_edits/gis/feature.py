# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:52:22 2021

@author: chris.kerklaan
"""

# First-party imports
import logging

# Third-party imports
import shapely
from osgeo import ogr, _ogr

# Local imports
from .fields import Fields, GeometryFields
from . import Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon

# GLOBALS
logger = logging.getLogger(__name__)
GEOMETRIES = [
    Point,
    MultiPoint,
    LineString,
    MultiLineString,
    Polygon,
    MultiPolygon,
]
OGR_GEOMETRY_CONVERSION = {}
for geometry in GEOMETRIES:
    OGR_GEOMETRY_CONVERSION.update(
        {coverage: geometry for coverage in geometry.ogr_coverage}
    )


class Feature:
    """wrapper of ogr feature

    Edits to geometry and fields are never set, just do object['field']
    to set the field
    """

    def __init__(self, feature: ogr.Feature, layer: ogr.Layer, fields: Fields):
        self.feature = feature

        if feature is not None:
            self.exists = True
        else:
            layer.ResetReading()
            first = None
            if layer.GetFeatureCount() > 0:
                first = layer.GetNextFeature().GetFID()
            raise IndexError(
                f"Feature does not exist, starts at: {first}, use vector.fids for all."
            )

        self.layer = layer
        self.fields = fields

    def __iter__(self):
        for key, value in self.items.items():
            yield key, value

    def __setitem__(self, field_name, value):
        _ogr.Feature_SetField(self.feature, field_name, value)
        _ogr.Layer_SetFeature(self.layer, self.feature)

    def __getitem__(self, field_name):
        return self.get_item(field_name)

    def __str__(self):
        return str(self.all_items)

    def __repr__(self):
        return str(self.all_items)

    def copy(self):
        return Feature(_ogr.Feature_Clone(self.feature))

    @property
    def all_items(self):
        return {"items": self.items, "geometry": self.geometry}

    @property
    def fid(self):
        return _ogr.Feature_GetFID(self.feature)

    @fid.setter
    def fid(self, value):
        _ogr.Feature_SetFID(value)
        _ogr.Layer_SetFeature(self.layer, self.feature)

    def get_item(self, field_name):
        i = self.fields.indices[field_name]
        if _ogr.Feature_IsFieldNull(self.feature, i):
            return None
        else:
            func = self.fields.get_field_func(field_name)
            return func(self.feature, i)

    def get_geometry(self):
        geometry = _ogr.Feature_GetGeometryRef(self.feature)
        if not geometry:
            geometry = _ogr.Feature_GetGeomFieldRef(self.feature, 1)
        return geometry

    @property
    def geometry(self):
        return geometry_factory(self.get_geometry())

    @geometry.setter
    def geometry(self, value):
        if isinstance(value, shapely.geometry.base.BaseGeometry):
            value = ogr.CreateGeometryFromWkb(value.wkb)

        _ogr.Feature_SetGeometry(self.feature, value)
        _ogr.Layer_SetFeature(self.layer, self.feature)

    @property
    def items(self):
        items = {}
        for field_name in self.fields.names:
            items[field_name] = self[field_name]
        return items

    @property
    def id(self):
        return self.fid

    @property
    def spatial_reference(self):
        return self.geometry.spatial_reference

    @property
    def area(self):
        return self.geometry.area

    @property
    def geometry_fields(self):
        return GeometryFields(self.feature)

    def simplify(self, wkb_threshold, start=10000000, quiet=False):
        """simplifying a geometry to a wkb size threshold"""
        geometry = self.geometry
        while geometry.WkbSize() > wkb_threshold:
            geometry = geometry.simply(1 / start)
            start = start / 10
        return geometry

    def reproject(self, out_epsg):
        return self.geometry.reproject(out_epsg)


def to_ogr(feature):
    """Returns an ogr feature"""
    if hasattr(feature, "feature"):
        return feature.feature
    return feature


def geometry_type_class(ogr_geometry_type):
    """returns the correct geometry type given an ogr_geometry"""
    try:
        geometry_class = OGR_GEOMETRY_CONVERSION[ogr_geometry_type]
    except KeyError:
        name = ogr.GeometryTypeToName(ogr_geometry_type)
        raise ValueError(f"Found unexpected geometry type {name}")
    else:
        return geometry_class


def geometry_factory(geometry: ogr.Geometry):
    """converts an ogr geometry to the custom ogr geometry classes"""
    if geometry:
        geometry_type = _ogr.Geometry_GetGeometryType(geometry)
        return geometry_type_class(geometry_type)(geometry)
