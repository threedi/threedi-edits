# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:14:38 2021

@author: chris.kerklaan
"""
# First-party imports
import logging

# Third-party imports
from osgeo import ogr

# optional imports
import importlib

HAS_SCIPY = importlib.util.find_spec("scipy") is not None
if HAS_SCIPY:
    from scipy import spatial

# Local imports
from .geometry import dissolve, Geometry
from .linestring import LineString, MultiLineString

# GLOBALS
logger = logging.getLogger(__name__)

POLYGON_COVERAGE = [
    ogr.wkbPolygon,
    ogr.wkbPolygon25D,
    ogr.wkbPolygonM,
    ogr.wkbPolygonZM,
    ogr.wkbCurvePolygon,
    ogr.wkbCurvePolygonM,
    ogr.wkbCurvePolygonZ,
    ogr.wkbCurvePolygonZM,
]

MULTIPOLYGON_COVERAGE = [
    ogr.wkbMultiPolygon,
    ogr.wkbMultiPolygon25D,
    ogr.wkbMultiPolygonM,
    ogr.wkbMultiPolygonZM,
    ogr.wkbMultiSurface,
    ogr.wkbMultiSurface,
    ogr.wkbMultiSurfaceM,
    ogr.wkbMultiSurfaceZ,
    ogr.wkbMultiSurfaceZM,
]


class Polygon(Geometry):
    ogr_coverage = POLYGON_COVERAGE

    def __init__(self, geometry: ogr.wkbPolygon = None, epsg=None):
        super().__init__(geometry)
        self.check_type(Polygon.ogr_coverage)
        if epsg:
            self.spatial_reference = epsg

    @classmethod
    def from_points(cls, points, flatten=True, close=True, epsg=None):
        """takes list of tuple points and creates an ogr polygon"""

        output_geom = ogr.Geometry(ogr.wkbPolygon)

        ring = ogr.Geometry(ogr.wkbLinearRing)
        for point in points:
            ring.AddPoint(*point)

        output_geom.AddGeometry(ring)

        if close:
            output_geom.CloseRings()

        if flatten:
            output_geom.FlattenTo2D()

        if not output_geom.IsValid():
            logger.warning(
                """
                            Is it a self-intersection polygon?
                            Are the points in the form a ring? E.g., 
                            left-upper, left-lower, right-lower, right-upper"""
            )

        return cls(output_geom, epsg)

    @classmethod
    def from_wkt(cls, wkt, epsg=None):
        """takes tuple points and creates an ogr point"""
        return cls(ogr.CreateGeometryFromWkt(wkt), epsg)

    @property
    def narrowity(self):
        """
        returns the area of the geometry divided by the length/perimeter.

        This is generally used as a pre-filter to derive long and narrow polygons.
        Interpretation depends on the use

        """

        return self.boundary.length / self.area

    @property
    def narrow(self):
        """
        returns True if self.narrowity is smaller than 0.35
        """
        return self.narrowity > 0.35

    def centerline(self, resolution=0.5, intersection=True):
        """
        returns the center of this polygon
        params:
            resolution: resolution of the voronoi boundary points
            intersection: cuts the multiline at intersections

        A centerline is returned by creating points on the boundary of the line.
        The points are used to create vonoroi polygons.
        Filtering on the ridges inside the polygon returns the centerline

        """
        return MultiLineString(centerline(self, resolution, intersection))

    def width(self, samples: int):
        """
        returns a list of possible widths, made along the centerline of the polygon
        params:
            samples: amount of samples taken along the centerline

        The direction of a polygon is hard to describe. We describe it as
        the centerline of the polygon. Perpendicular lines are taken
        and clipped to the geometry to sample some widths of the polygon.
        Which is returned in this function

        It is a very costly operation, so use with care.

        """
        return width(self, samples)


class MultiPolygon(Geometry):
    ogr_coverage = MULTIPOLYGON_COVERAGE

    def __init__(
        self, geometry: ogr.wkbMultiPolygon = None, points: list = None, epsg=None
    ):

        if points:
            geometry = self.create_multipolygon(points)

        super().__init__(geometry)

        if epsg:
            self.spatial_reference = epsg

        # check types
        if self.type == ogr.wkbPolygon:
            raise TypeError(
                """geometry is multi while vector is single
                                use vector.to_single()
                            """
            )
        self.check_type(MultiPolygon.ogr_coverage)

    @classmethod
    def from_points(cls, points, flatten=True, close=True, epsg=None):
        output_geom = ogr.Geometry(ogr.wkbMultiPolygon)

        for point in points:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for p in point:
                ring.AddPoint(p)

            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            output_geom.AddGeometry(poly)

        if close:
            output_geom.CloseRings()

        if flatten:
            output_geom.FlattenTo2D()

        if not output_geom.IsValid():
            logger.warning(
                """
                            Is it a self-intersection polygon?
                            Are the points in the form a ring? E.g., 
                            left-upper, left-lower, right-lower, right-upper"""
            )
        return cls(output_geom, epsg)

    @classmethod
    def from_wkt(cls, wkt, epsg=None):
        """takes tuple points and creates an ogr point"""
        return cls(ogr.CreateGeometryFromWkt(wkt), epsg)


def centerline(geometry: Polygon, resolution: int, intersection: bool, clean=False):
    if not HAS_SCIPY:
        raise ImportError("Please install scipy before using this function")

    boundary = geometry.boundary
    points = LineString(boundary).points_on_line(resolution)

    vor = spatial.Voronoi([p.point for p in points])
    vertices = vor.vertices
    ridges = vor.ridge_vertices

    multi = ogr.Geometry(5)
    for ridge in ridges:
        if "-1" not in ridge:  # is infinite
            start_point = (vertices[ridge[0]][0], vertices[ridge[0]][1])
            end_point = (vertices[ridge[1]][0], vertices[ridge[1]][1])
            line = LineString.from_points([start_point, end_point])

            if line.Within(geometry):
                multi.AddGeometry(line)

    dissolved = dissolve(multi)
    if clean:
        return MultiLineString(dissolved).clean(intersection)
    return MultiLineString(dissolved)


def width(geometry: Polygon, samples: int = 10):
    centerlines = geometry.centerline(0.5, intersection=True)

    # keep the longest line
    for i, centerline in enumerate(centerlines):
        if i == 0:
            line = centerline
        else:
            if centerline.Length() > line.Length():
                line = centerline

    line = LineString(line)
    perpendicular_lines = line.perpendicular_lines(
        line.length / (samples + 1), line.length, start=False, end=False
    )

    widths = []
    for line in perpendicular_lines:
        intersect = line.clip(geometry)
        if intersect:
            widths.append(intersect[0].length)
    return widths
