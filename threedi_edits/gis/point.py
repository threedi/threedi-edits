# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:13:51 2021

@author: chris.kerklaan
"""
# Third-party imports
from osgeo import ogr
from shapely.ops import voronoi_diagram

# Local imports
from .geometry import Geometry

# Globals
POINT_COVERAGE = [
    ogr.wkbPoint,
    ogr.wkbPoint25D,
    ogr.wkbPointM,
    ogr.wkbPointZM,
]
MULTIPOINT_COVERAGE = [
    ogr.wkbMultiPoint,
    ogr.wkbMultiPoint25D,
    ogr.wkbMultiPointM,
    ogr.wkbMultiPointZM,
]


class Point(Geometry):
    ogr_coverage = POINT_COVERAGE

    def __init__(self, geometry: ogr.wkbPoint = None, epsg=None, flatten=True):
        if flatten:
            geometry.FlattenTo2D()

        super().__init__(geometry)

        self.check_type(Point.ogr_coverage)
        if epsg:
            self.spatial_reference = epsg

    @classmethod
    def from_point(cls, points: tuple = None, flatten=True, epsg=None):
        """takes tuple points and creates an ogr point"""
        output_geom = ogr.Geometry(ogr.wkbPoint)
        output_geom.AddPoint(*points)

        return cls(output_geom, epsg, flatten)

    @classmethod
    def from_wkt(cls, wkt, flatten=True, epsg=None):
        """takes tuple points and creates an ogr point"""
        return cls(ogr.CreateGeometryFromWkt(wkt), epsg, flatten)

    @property
    def point(self):
        return self.points[0]


class MultiPoint(Geometry):
    ogr_coverage = MULTIPOINT_COVERAGE

    def __init__(self, geometry: ogr.wkbMultiPoint, epsg=None):
        super().__init__(geometry)
        if epsg:
            self.spatial_reference = epsg

        if self.type == ogr.wkbPoint:
            raise TypeError(
                """geometry is multi while vector is single
                                use vector.to_single()
                            """
            )
        self.check_type(MultiPoint.ogr_coverage)

    @classmethod
    def from_points(cls, points, flatten=True, epsg=None):
        """takes list of tuple points and creates an ogr multipoint"""
        output_geom = ogr.Geometry(ogr.wkbMultiPoint)
        for point in points:
            point_geometry = ogr.Geometry(ogr.wkbPoint)
            point_geometry.AddPoint(*point)
            output_geom.AddGeometry(point_geometry)

        if flatten:
            output_geom.FlattenTo2D()

        return cls(output_geom, epsg)

    @classmethod
    def from_wkt(cls, wkt, flatten=True, epsg=None):
        """takes tuple points and creates an ogr point"""
        return cls(ogr.CreateGeometryFromWkt(wkt), epsg, flatten)

    @property
    def points_geometry(self):
        return [Point(point) for point in self]

    @property
    def points(self):
        return [point.GetPoint() for point in self]

    @property
    def point(self):
        return self.points[0]

    def voronoi(self, envelope=None, tolerance=0.0, edges=False):
        """returns a tuple of list of regions and shapely collection"""
        data = voronoi_diagram(self.shape, envelope, tolerance, edges)
        regions = [self.from_shapely(r) for r in data.geoms]
        return regions, data
