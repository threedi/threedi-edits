# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 13:49:48 2021

@author: chris.kerklaan

Wrapper for an ogr geometry

#TODO:
    1. Take over the most important shapely functions.
    1. Convert/fix geometry should have an input geom_type and a target_geom_type
    
    3. Take over shapely's affine functions.
    4. Take over triangulate.
    5. Take over vonoroi shapely.ops.voronoi_diagram.
    6. Snapping.
    
"""
# First-party imports
import logging
import numpy as np

# Third-party imports
import shapely
import shapely.lib as shapely_lib
from osgeo import ogr, _ogr
from cached_property import cached_property

# Local imports
from .spatial_reference import SpatialReference

# GLOBALS
logger = logging.getLogger(__name__)

POLYGON = "POLYGON (({x1} {y1},{x2} {y1},{x2} {y2},{x1} {y2},{x1} {y1}))"
POINT = "POINT ({x1} {y1})"
LINESTRING = "LINESTRING ({x1} {y1}, {x2} {y2})"

TRANSLATION = {
    "point": ogr.wkbPoint,
    "line": ogr.wkbLineString,
    "polygon": ogr.wkbPolygon,
    "multipoint": ogr.wkbMultiPoint,
    "multiline": ogr.wkbMultiLineString,
    "multiPolygon": ogr.wkbMultiPolygon,
}

SINGLE_TYPES = {
    ogr.wkbMultiSurface: ogr.wkbPolygon,
    ogr.wkbMultiPolygon: ogr.wkbPolygon,
    ogr.wkbMultiLineString: ogr.wkbLineString,
    ogr.wkbMultiPoint: ogr.wkbPoint,
    ogr.wkbPolygon: ogr.wkbPolygon,
    ogr.wkbLineString: ogr.wkbLineString,
    ogr.wkbPoint: ogr.wkbPoint,
}

MULTI_TYPES = {
    ogr.wkbPoint: ogr.wkbMultiPoint,
    ogr.wkbLineString: ogr.wkbMultiLineString,
    ogr.wkbPolygon: ogr.wkbMultiPolygon,
    ogr.wkbMultiPoint: ogr.wkbMultiPoint,
    ogr.wkbMultiLineString: ogr.wkbMultiLineString,
    ogr.wkbMultiPolygon: ogr.wkbMultiPolygon,
}

# TODO Can be written easier, also does not cover everything.
ESTIMATES = {
    # multipolygons
    ogr.wkbMultiPolygon: ogr.wkbMultiPolygon,
    ogr.wkbMultiSurface: ogr.wkbMultiPolygon,
    ogr.wkbMultiSurfaceM: ogr.wkbMultiPolygon,
    ogr.wkbMultiSurfaceZ: ogr.wkbMultiPolygon,
    ogr.wkbMultiSurfaceZM: ogr.wkbMultiPolygon,
    ogr.wkbMultiPolygon25D: ogr.wkbMultiPolygon,
    ogr.wkbMultiPolygonM: ogr.wkbMultiPolygon,
    ogr.wkbMultiPolygonZM: ogr.wkbMultiPolygon,
    # polygons
    ogr.wkbPolygon: ogr.wkbPolygon,
    ogr.wkbCurvePolygon: ogr.wkbPolygon,
    ogr.wkbCurvePolygonM: ogr.wkbPolygon,
    ogr.wkbCurvePolygonZ: ogr.wkbPolygon,
    ogr.wkbCurvePolygonZM: ogr.wkbPolygon,
    ogr.wkbSurface: ogr.wkbPolygon,
    ogr.wkbSurfaceM: ogr.wkbPolygon,
    ogr.wkbSurfaceZ: ogr.wkbPolygon,
    ogr.wkbSurfaceZM: ogr.wkbPolygon,
    ogr.wkbPolygon25D: ogr.wkbPolygon,
    ogr.wkbPolygonM: ogr.wkbPolygon,
    ogr.wkbPolygonZM: ogr.wkbPolygon,
    # multilinestring
    ogr.wkbMultiLineString: ogr.wkbMultiLineString,
    ogr.wkbMultiCurve: ogr.wkbMultiLineString,
    ogr.wkbMultiCurveM: ogr.wkbMultiLineString,
    ogr.wkbMultiCurveZ: ogr.wkbMultiLineString,
    ogr.wkbMultiCurveZM: ogr.wkbMultiLineString,
    ogr.wkbMultiLineString25D: ogr.wkbMultiLineString,
    ogr.wkbMultiLineStringM: ogr.wkbMultiLineString,
    ogr.wkbMultiLineStringZM: ogr.wkbMultiLineString,
    # linestring
    ogr.wkbLineString: ogr.wkbLineString,
    ogr.wkbCircularString: ogr.wkbLineString,
    ogr.wkbCircularStringM: ogr.wkbLineString,
    ogr.wkbCircularStringZ: ogr.wkbLineString,
    ogr.wkbCircularStringZM: ogr.wkbLineString,
    ogr.wkbCurve: ogr.wkbLineString,
    ogr.wkbCurveM: ogr.wkbLineString,
    ogr.wkbCurveZ: ogr.wkbLineString,
    ogr.wkbCurveZM: ogr.wkbLineString,
    ogr.wkbLineString25D: ogr.wkbLineString,
    ogr.wkbLineStringM: ogr.wkbLineString,
    ogr.wkbLineStringZM: ogr.wkbLineString,
    # multipoints
    ogr.wkbMultiPoint: ogr.wkbMultiPoint,
    ogr.wkbMultiPoint25D: ogr.wkbMultiPoint,
    ogr.wkbMultiPointM: ogr.wkbMultiPoint,
    ogr.wkbMultiPointZM: ogr.wkbMultiPoint,
    # point
    ogr.wkbPoint: ogr.wkbPoint,
    ogr.wkbPoint25D: ogr.wkbPoint,
    ogr.wkbPointM: ogr.wkbPoint,
    ogr.wkbPointZM: ogr.wkbPoint,
}

SINGLE_TYPES_INVERT = dict(map(reversed, SINGLE_TYPES.items()))
SINGLE_TO_MULTIPLE = {1: [1, 4], 2: [2, 5], 3: [3, 6], 4: [4], 5: [5], 6: [6]}
ALLOWED = {1: [1, 4], 2: [2, 5], 3: [3, 6], 4: [4, 1], 5: [5, 2], 6: [6, 3]}
SHAPELY_INVALID = np.uint8(2)


class Geometry(ogr.Geometry):
    """
    The Geometry class is an inhertied version of an ogr geometry.

    Memory:
        The input is copied, hence the pointer is lost.
        The geometry in a feature must be set with set_feature.
        Every function returns a geometry.

    """

    def __init__(self, geometry: ogr.Geometry):
        super().__init__(wkb=convert_geometry(geometry).ExportToWkb())
        self.AssignSpatialReference(geometry.GetSpatialReference())

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        wkt = self.wkt
        max_length = 60
        if len(wkt) > max_length:
            return "<{}...>".format(wkt[: max_length - 3])

        return f"<{wkt}>"

    @classmethod
    def from_shapely(cls, geometry):
        return cls(ogr.CreateGeometryFromWkt(geometry.wkt))

    def check_type(self, types):
        if not self.type in types:
            raise TypeError(f"Wrong ogr geom type {self.name} {self.type_name}")

    @cached_property
    def shape(self):
        """laods a shapely geometry from wkt"""
        geometry = shapely.from_wkt(self.wkt)
        shapely.prepare(geometry)
        return geometry

    @property
    def area(self):
        return self.GetArea()

    @property
    def length(self):
        return self.Length()

    @property
    def name(self):
        return self.GetGeometryName()

    @property
    def type(self):
        return self.GetGeometryType()

    @property
    def is_multi(self):
        return self.type > 3

    @property
    def type_name(self):
        return ogr.GeometryTypeToName(self.type)

    @property
    def valid(self):
        return self.IsValid()

    @property
    def fixed(self):
        return fix(self)

    @property
    def points(self):
        return self.GetPoints()

    @property
    def envelope(self):
        return self.GetEnvelope()

    @property
    def wkt(self):
        return self.ExportToWkt()

    @property
    def wkb(self):
        return self.ExportToWkb()

    @property
    def wkb_size(self):
        return self.WkbSize()

    @property
    def point_on_surface(self):
        return Geometry(self.PointOnSurface())

    @property
    def centroid(self):
        return Geometry(self.Centroid())

    @property
    def boundary(self):
        return Geometry(self.GetBoundary())

    @property
    def spatial_reference(self):
        return SpatialReference.from_sr(self.GetSpatialReference(), self.centroid)

    @spatial_reference.setter
    def spatial_reference(self, sr_input):
        """takes an epsg or a spatial_reference input"""
        if type(sr_input) == int:
            sr_input = SpatialReference.from_epsg(sr_input)
        self.AssignSpatialReference(sr_input)

    @property
    def epsg(self):
        return self.spatial_reference.epsg

    @property
    def ogr(self):
        """returns a pure ogr geometry"""
        return release(self)

    def copy(self):
        return Geometry(geometry=self)

    def buffer(self, size):
        return Geometry(self.Buffer(size))

    def simplify(self, simplification):
        return Geometry(self.Simplify(simplification))

    def difference(self, geometry, ignore_errors=False):
        """returns a differenced geometry, or none if invalid or not intersects"""
        return Geometry(difference(self, geometry, ignore_errors))

    def clip(self, mask):
        """returns a list of clipped geometry in the correct output type
        or none
        """
        clipped = clip(self, mask, self.type)
        if clipped:
            return [Geometry(i) for i in clipped]
        else:
            return

    def reproject(self, out_epsg, use_points=False):
        """reprojects the geometries"""
        geometry = self.ogr
        geometry.Transform(self.spatial_reference.transform(out_epsg))
        return Geometry(geometry)

    def dissolve(self):
        """dissolved geometry, returns a single geometry"""
        return Geometry(dissolve(self), SINGLE_TYPES[self.type])

    def to_single(self):
        """creates single parts, returns a list or none"""
        return [Geometry(i) for i in multipart_to_singlepart(self)]

    def fix(self, geometry_type=None):
        return fix(self, geometry_type)

    def union(self, geometry):
        return Geometry(self.Union(geometry))

    def intersection(self, geometry):
        return Geometry(self.Intersection(geometry))


def shapely_geometry(geometry: ogr.Geometry) -> shapely.geometry.base.BaseGeometry:
    if geometry:
        geometry = convert_geometry(geometry)
        wkb = bytes(_ogr.Geometry_ExportToWkb(geometry))
        array = np.asarray(wkb, dtype=object)
        geometry = shapely_lib.from_wkb(array, SHAPELY_INVALID)
        shapely_lib.prepare(geometry)
        return geometry
    return None


def to_ogr(geometry):
    if isinstance(geometry, shapely.geometry.base.BaseGeometry):
        return ogr.CreateGeometryFromWkt(geometry.wkt)
    if hasattr(geometry, "ogr"):
        return geometry.ogr
    if isinstance(geometry, ogr.Geometry):
        return geometry
    raise ValueError("Cannot convert ogr geometry")


def to_shapely(geometry):
    if isinstance(geometry, ogr.Geometry):
        return shapely_geometry(geometry)

    elif hasattr(geometry, "ogr"):
        return shapely_geometry(geometry.ogr)

    elif hasattr(geometry, "shape"):
        return geometry.shape

    else:
        if isinstance(geometry, shapely.geometry.base.BaseGeometry):
            return geometry

    raise ValueError("Cannot convert to shapely geometry")


def release(geometry: ogr.Geometry):
    """releases the c pointer from the geometry"""
    output = ogr.CreateGeometryFromWkb(geometry.ExportToWkb())
    output.AssignSpatialReference(geometry.GetSpatialReference())
    return output


def estimate_geometry_type(geometry_type):
    """
    Estimates the target geometry type
    """
    if geometry_type > ogr.wkbUnknown and geometry_type < ogr.wkbMultiPolygon:
        return geometry_type

    if not geometry_type in ESTIMATES:
        name = ogr.GeometryTypeToName(geometry_type)
        raise ValueError(f"Cannot estimate geometry type for: {name}")

    return ESTIMATES[geometry_type]


def convert_geometry(geometry, target_type=None) -> ogr.Geometry:
    """
    We only use planar and linear geometries.
    Therefore we convert the geometry to the target type.
    if the target type is not present we estimate it.

    We account for the following cases:
        - geometry type is 3D.
        - geometry type is curved.
        - geometry type is a collecton.

    returns a single, linear and planar geometry, or raises an error.


    """
    geometry_type = geometry.GetGeometryType()

    if target_type is None:
        target_type = estimate_geometry_type(geometry_type)

    if geometry_type == target_type:
        return geometry

    # 3D conversion
    if geometry.Is3D():
        geometry.FlattenTo2D()

    # Curved geometrya
    if geometry.HasCurveGeometry():
        geometry = geometry.GetLinearGeometry()

    # Geometry type, also looks for Z, M and ZM.
    if str(ogr.wkbGeometryCollection) in str(geometry_type):
        logger.debug("Found a geometry collection, retrieving only the first geometry")
        for part in geometry:
            if part.GetGeometryType() == target_type:
                geometry = part
                break

    # Finally, if all of the above does not work, we try to force.
    geometry = ogr.ForceTo(geometry, target_type)

    return geometry


def fix(geometry, target_geometry_type=None):
    """
    Fixes a geometry is it is invalid or has the wrong geometry type.
    Fixes of the wrong type using convert_geometry.
    Things that we fix:
        - Point count of 1 in a linestring returns an error.
        - Invalid polygons of multipolygons receive a buffer.

    Everything else is done with make valid.


    """
    if geometry is None:
        return None

    geometry_type = geometry.GetGeometryType()
    geometry_name = ogr.GeometryTypeToName(geometry_type)

    if not target_geometry_type:
        target_geometry_type = estimate_geometry_type(geometry_type)

    geometry = convert_geometry(geometry, target_geometry_type)

    if geometry.IsValid():
        return geometry

    if geometry_type == ogr.wkbLineString:
        if geometry.GetPointCount() == 1:
            logger.debug("Geometry point count of linestring = 1")
            return geometry, False

    if geometry_type in [ogr.wkbPolygon, ogr.wkbMultiPolygon]:
        logger.debug("Buffered the geometry")
        geometry = geometry.Buffer(0)

    if not geometry.IsValid():
        logger.deug("Using MakeValid")
        geometry = geometry.MakeValid()  # multipolygon to collection.

    if not geometry.IsValid():
        logger.warning("Failed to fix the geometry, returning None.")
        return None

    if not geometry.GetGeometryType() == target_geometry_type:
        target_name = ogr.GeometryTypeToName(target_geometry_type)
        logger.warning(f"Could not convert {geometry_name} to {target_name}.")

    return geometry


def multipart_to_singlepart(geometry: ogr.Geometry):
    """returns a list of geometries, if invalid none"""
    if "MULTI" in geometry.GetGeometryName():
        return [release(part) for part in geometry]
    else:
        return [geometry]


def clip(geometry: ogr.Geometry, mask: ogr.Geometry, output_type: int):
    """Clips a single geometry on a mask
    returns invalid geoms as none
    returns a geometry of an allowed output type
    returns a list
    """
    within = geometry.Within(mask)

    # check for the boundary
    if geometry.Intersects(mask) and not within:
        intersect = geometry.Intersection(mask)
        intersect_type = intersect.GetGeometryType()
        if intersect_type not in ALLOWED[output_type]:
            logger.debug("Output has incorrect geometry type")
            return

        if intersect_type > 3:
            return [intersect]
        else:
            return multipart_to_singlepart(intersect)

    elif within:
        return [geometry]

    else:
        logger.debug("Outside of geometry")
        return  # outside


def difference(geometry: ogr.Geometry, mask: ogr.Geometry, ignore_errors=False):
    """checks if intersects before doing difference"""

    try:
        if geometry.Intersects(mask):
            return geometry.Difference(mask)
    except RuntimeError:
        if not ignore_errors:
            raise RuntimeError("Difference")
        logger.debug("Geometry difference runtime error - ignore = True")

    # if not intersects or ignore errors return original geometry
    return geometry


def dissolve(multigeometry: ogr.Geometry):
    """dissolves a multigeometry into a single geometry"""
    geometry_type = multigeometry.GetGeometryType()

    if geometry_type == ogr.wkbMultiPolygon:
        union = multigeometry.UnionCascaded()
    else:
        union = ogr.Geometry(SINGLE_TYPES[geometry_type])
        for geometry in multigeometry:
            union = union.Union(geometry)

    return union
