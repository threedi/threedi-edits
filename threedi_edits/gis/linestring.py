# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:14:17 2021

@author: chris.kerklaan
"""
# First-party import
import math
import logging
from copy import deepcopy
from operator import attrgetter
from collections import namedtuple

# Third-party imports
import numpy as np
from osgeo import ogr
from itertools import chain

# Local imports
from .geometry import Geometry
from .point import Point, MultiPoint

# Globals
LINESTRING_COVERAGE = [
    ogr.wkbLineString,
    ogr.wkbLineString25D,
    ogr.wkbLineStringM,
    ogr.wkbLineStringZM,
]
MULTILINESTRING_COVERAGE = [
    ogr.wkbMultiLineString,
    ogr.wkbMultiLineString25D,
    ogr.wkbMultiLineStringM,
    ogr.wkbMultiLineStringZM,
]


# Logger
logger = logging.getLogger(__name__)


class LineString(Geometry):

    ogr_coverage = LINESTRING_COVERAGE

    def __init__(self, geometry: ogr.wkbLineString = None, epsg=None):
        super().__init__(geometry)
        self.check_type(LineString.ogr_coverage)
        if epsg:
            self.spatial_reference = epsg

    @classmethod
    def from_points(cls, points, flatten=True, epsg=None):
        """takes list of tuple points and creates an ogr linestring"""
        output_geom = ogr.Geometry(ogr.wkbLineString)
        for point in points:
            output_geom.AddPoint(*point)
        if flatten:
            output_geom.FlattenTo2D()

        return cls(output_geom, epsg)

    @classmethod
    def from_wkt(cls, wkt, flatten=True, epsg=None):
        """takes tuple points and creates an ogr point"""
        return cls(ogr.CreateGeometryFromWkt(wkt), epsg, flatten)

    def __iter__(self):
        """iterates over the vertices"""
        points = self.points
        for pid in range(0, len(points) - 1):
            yield LineString.from_points([points[pid], points[pid + 1]])

    @property
    def points_geometry(self):
        return [Point.from_point(point) for point in self.points]

    @property
    def start_point(self):
        return self.points_geometry[0]

    @property
    def end_point(self):
        return self.points_geometry[-1]

    @property
    def middle_point(self, interval=1):
        points = self.points_on_line(interval, vertices=True)
        return points[int(len(points) / 2)]

    def reversed(self):
        """reverses the line"""
        points = self.points
        points.reverse()
        return LineString.from_points(points)

    def intersection(self, geometry):
        """intersection of a linestring, returns a point, multipoint or none"""
        if self.Intersects(geometry):
            intersection = self.Intersection(geometry)

            if intersection.GetGeometryType() == ogr.wkbPoint:
                return Point(intersection)
            elif intersection.GetGeometryType() == ogr.wkbLineString:
                return LineString(intersection)
            else:
                return MultiPoint(intersection)

        else:
            logger.debug("Found no intersection")
            return None

    def points_on_line(
        self,
        interval=1,
        custom_interval=None,
        start=True,
        end=True,
        vertices=False,
        geometry=True,
    ):
        """Return points on a linestring

        Params:
            interval: distance between poitns
            custom_interval: a list of distances
            start: include start point
            end: include end point
            vertices: include vertices
            geometry: returns as a point
        """

        return points_on_line(
            self.points,
            interval,
            custom_interval,
            start,
            end,
            vertices,
            geometry,
        )

    def perpendicular_lines(
        self,
        distance,
        perpendicular_length,
        start=True,
        end=True,
        vertices=False,
    ):

        """returns a perpendicular linestring on point1 with a certain dist
        params:
            distance: distance between perpendicular lines
            perpendicular_length: length of the lines
        """
        return perpendicular_lines(
            self,
            distance,
            perpendicular_length,
            start=start,
            end=end,
            vertices=vertices,
        )

    def add_vertice(self, geometry: Point):
        """
        Returns a linestring with an extra vertice on the point
        If the vertice is already present, we will not adjust it.

        most likely there are two types of linestrings
        1. With one coordinate
        2. With a connection between coordinates
        this add vertice is for assumption 1.

        """
        # First we find all intersects
        distances = [(i, v.Distance(geometry)) for i, v in enumerate(self)]

        distances.sort(key=lambda tup: tup[-1])
        closest = distances[0]

        points = []
        snap_count = 0
        for i, vertice in enumerate(self):
            vertice_point = vertice.points[0]
            if vertice_point not in points:
                points.append(vertice_point)

            if i == closest[0]:
                points.append(geometry.point)
                snap_count += 1

        # last one
        if vertice.points[-1] not in points:
            points.append(vertice.points[-1])

        if snap_count == 0:
            logger.debug("snap not found")

        return LineString.from_points(points)

    def split(self, interval=1, custom_interval=None):
        points = self.points_on_line(interval, custom_interval, geometry=False)
        return [
            LineString.from_points([points[pid], points[pid + 1]])
            for pid in range(len(points) - 1)
        ]

    def split_on_vertices(self):
        return [vertice for vertice in self]

    def split_on_point(self, geometry: Point):

        points = self.points
        assert geometry.point in points

        index = points.index(geometry.point)
        line1 = LineString.from_points(points[: index + 1])
        line2 = LineString.from_points(points[index:])
        return [line1, line2]

    def as_multi(self):
        return MultiLineString.from_points(self.points)

    def transform(self, epsg):
        return LineString(self.reproject(epsg))

    def angle(self):
        return angle(*self.points)


class MultiLineString(Geometry):
    """Can also use linestrings"""

    ogr_coverage = MULTILINESTRING_COVERAGE

    def __init__(self, geometry: ogr.wkbMultiLineString = None, epsg=None):
        super().__init__(geometry)
        if epsg:
            self.spatial_reference = epsg

        self.check_type(MultiLineString.ogr_coverage)

    @classmethod
    def from_points(cls, points, flatten=True, epsg=None):
        """creates a ogr multilinestring from points"""

        output_geom = ogr.Geometry(ogr.wkbMultiLineString)

        # check if input is a single linestring
        if type(points[0]) == tuple:
            points = [[point] for point in points]

        for point in points:
            line = ogr.Geometry(ogr.wkbLineString)
            for p in point:
                line.AddPoint(*p)
            output_geom.AddGeometry(line)

        if flatten:
            output_geom.FlattenTo2D()

        return cls(output_geom, epsg)

    @classmethod
    def from_wkt(cls, wkt, flatten=True, epsg=None):
        """takes tuple points and creates an ogr point"""
        return cls(ogr.CreateGeometryFromWkt(wkt), epsg, flatten)

    def dissolve(self):
        return self.clean(False).to_single()[0]

    @property
    def points(self):  # overwrits the property of geometry class
        return [l.GetPoints() for l in self]

    def clean(self, intersection):
        """
        Returns a rearranged multilinestrings

        params:
            intersection: If True, splits the linestring on large intersects

        A multiline can be a multiline even if there are double
        coordinates and they are not arranged.
        This function arranges the multilinestring and cleans out all
        double points.

        """
        return MultiLineString(clean_multilinestring(self, intersection))

    def points_on_line(
        self,
        interval=1,
        custom_interval=None,
        start=True,
        end=True,
        vertices=False,
        geometry=True,
        single=True,
    ):

        all_points = []
        for points in self.points:
            all_points.append(
                points_on_line(
                    points,
                    interval,
                    custom_interval,
                    start,
                    end,
                    vertices,
                    geometry,
                )
            )

        if single:
            return list(chain.from_iterable(all_points))
        else:
            return all_points

    def perpendicular_lines(
        self,
        distance,
        perpendicular_length,
        start=True,
        end=True,
        vertices=False,
        single=True,
    ):

        """returns a perpendicular linestring on point1 with a certain dist
        params:
            distance: distance between perpendicular lines
            perpendicular_length: length of the lines
        """
        lines = []
        for linestring in self:
            lines.extend(
                perpendicular_lines(
                    self,
                    distance,
                    perpendicular_length,
                    start=start,
                    end=end,
                    vertices=vertices,
                )
            )

        if single:
            return list(chain.from_iterable(lines))
        else:
            return lines


def calc_dist(x1, y1, x2, y2):
    """returns the distance beteen two points"""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def lower_bound(x, l):
    if l[0] > x and not l[0] == x:
        return
    for i, y in enumerate(l):
        if y > x:
            return l[i - 1]


def points_with_distance(point1, point2, distance):
    """
    Returns a point on a certain distance between two points
    Note that the distance is the distance from the first point in the coordinate unit
    """
    t = distance / calc_dist(*point1, *point2)
    p = (
        ((1 - t) * point1[0] + t * point2[0]),
        ((1 - t) * point1[1] + t * point2[1]),
    )

    if type(p[0]) is np.ndarray:
        return t, (p[0][0], p[1][0])
    else:
        return t, p


def points_on_line(
    points,
    interval=1,
    custom_interval: list = None,
    start=True,
    end=True,
    vertices=False,
    geometry=True,
):
    """
    returns a point on the line for every interval that is given
    use combine to also include begin, endpoint and vertices
    Params:
        interval: distance between poitns
        custom_interval: a list of distances
        start: include start point
        end: include end point
        vertices: include vertices
        geometry: returns as a point
    """

    sections = list(
        np.cumsum(
            [
                calc_dist(*points[pid], *points[pid + 1])
                for pid in range(0, len(points) - 1)
            ]
        )
    )
    total_dist = sections[-1]

    if custom_interval is not None:
        new_sections = custom_interval
    else:
        new_sections = np.linspace(interval, total_dist, int(total_dist / interval))

    if vertices:
        new_sections = sorted(list(new_sections) + sections)

    new_points = []
    for i in new_sections:
        bound = lower_bound(i, sections)
        if not bound:
            index = 0
            dist = i
        else:
            index = sections.index(bound) + 1
            dist = i - sections[sections.index(bound)]

        ratio, point = points_with_distance(
            points[index], points[index + 1], distance=dist
        )

        if 0 <= ratio <= 1:
            new_points.append(point)
        else:
            pass

    if start:
        new_points.insert(0, points[0])

    if end:
        new_points.append(points[-1])

    if geometry:
        return [Point.from_point(point) for point in new_points]

    return new_points


def angle(pt1, pt2):
    x_diff = pt2[0] - pt1[0]
    y_diff = pt2[1] - pt1[1]
    return math.degrees(math.atan2(y_diff, x_diff))


def perpendicular_points(pt, bearing, dist):
    """returns perpendicular points at point"""
    bearing_pt1 = math.radians(bearing + 90)
    bearing_pt2 = math.radians(bearing - 90)
    points = []
    for bearing in [bearing_pt1, bearing_pt2]:
        x = pt[0] + dist * math.cos(bearing)
        y = pt[1] + dist * math.sin(bearing)
        points.append((x, y))
    return points


def perpendicular_line(pt1, pt2, dist):
    return LineString.from_points(perpendicular_points(pt1, angle(pt1, pt2), dist))


def perpendicular_lines(
    linestring,
    line_dist=10,
    perp_dist=10,
    start=True,
    end=True,
    vertices=False,
):
    """returns perpendicular lines on a linestring"""
    points = linestring.points_on_line(
        interval=line_dist,
        custom_interval=None,
        start=start,
        end=end,
        vertices=vertices,
        geometry=False,
    )

    lines = [
        perpendicular_line(points[index], points[index + 1], perp_dist)
        for index in range(0, len(points) - 1)
    ]

    # add last point
    if len(points) > 1:
        lines.append(perpendicular_line(points[-1], points[-2], perp_dist))

    return lines


def arrange_multiline(points):
    """arranges the points in a multilinestring
    e.g.,
    [[2, 1] ,[4, 5], [0, 1], [5, 6], [9, 8]]

    to
    [[0, 1], [1, 2], [4, 5], [5, 6], [8, 9]]

    """

    def distance(coords1, coords2):
        return math.sqrt(
            math.pow(coords1[0] - coords2[0], 2) + math.pow(coords1[1] - coords2[1], 2)
        )

    MinDistance = namedtuple("MinDistance", "target distance offset reverse_target")
    lines = deepcopy(points)
    sorted_multistring = [lines.pop(0)]

    while lines:
        min_distances = []

        for line in lines:
            source_a = sorted_multistring[0][0]
            source_b = sorted_multistring[-1][-1]

            target_a = line[0]
            target_b = line[-1]

            distances = [
                MinDistance(
                    target=line,
                    distance=distance(source_b, target_a),
                    offset=1,
                    reverse_target=False,
                ),
                MinDistance(
                    target=line,
                    distance=distance(source_a, target_a),
                    offset=-1,
                    reverse_target=True,
                ),
                MinDistance(
                    target=line,
                    distance=distance(source_b, target_b),
                    offset=1,
                    reverse_target=True,
                ),
                MinDistance(
                    target=line,
                    distance=distance(source_a, target_b),
                    offset=-1,
                    reverse_target=False,
                ),
            ]

            min_distance = min(distances, key=attrgetter("distance"))
            min_distances.append(min_distance)

        min_distance = min(min_distances, key=attrgetter("distance"))
        target = min_distance.target

        if min_distance.reverse_target:
            target.reverse()

        if min_distance.offset == 1:
            sorted_multistring.append(target)
        else:
            sorted_multistring.insert(0, target)

        lines.remove(target)

    return sorted_multistring


def clean_multilinestring(geometry, intersections=True):
    """cleans a messy linestring order and turns them into logical
    seperations.

    returns a list of linestrings

    intersections:
        Splits the lines at the intersctions of multiple lines

    """

    points = MultiLineString(geometry).points
    if len(points) == 1:
        return geometry

    arranged = arrange_multiline(points)

    # combine the the points of the arranged points,
    lines = []
    line = [arranged[0][0], arranged[0][1]]
    for i in range(len(arranged) - 1):
        coord2 = arranged[i][1]
        coord3, coord4 = arranged[i + 1][0], arranged[i + 1][1]

        match = coord2 == coord3
        if match:  # we have a match
            line.append(coord4)

        final_lap = i == len(arranged) - 2
        if not match:
            lines.append(line)
            line = [coord3, coord4]
        elif match and final_lap:
            lines.append(line)

    lines = [LineString.from_points(l) for l in lines]

    if intersections:
        # get intersects    intersects = []
        line_geometries = [(i, l) for i, l in enumerate(lines)]
        final_geometries = []
        for i1, line1 in line_geometries:

            new_lines = [line1]
            for i2, line2 in line_geometries:
                if i1 == i2:
                    continue

                # find intersections of the whole line
                intersects = []
                for new_line in new_lines:

                    if new_line.GetPointCount() < 2:
                        continue

                    intersect = new_line.intersection(line2)

                    if not intersect:
                        continue

                    points = []
                    if intersect.is_multi:
                        points = [p for p in intersect.points_geometry]
                    else:
                        points = [intersect]

                    for point in points:
                        on_start_point = point == new_line.start_point.point
                        on_end_point = point == new_line.end_point.point

                        if on_start_point or on_end_point:
                            continue

                        intersects.append(point)

                # split the line per point and add to the main line list
                line_count = len(new_lines)
                while len(intersects) > 0:
                    point = intersects.pop(0)

                    for i in range(line_count):
                        new_line = new_lines[i]
                        if len(new_line.GetPoints()) <= 2:
                            continue

                        if point.Intersects(new_line):
                            for geometry in new_line.split_on_point(point):
                                new_lines.append(geometry)

            if len(new_lines) > 1:
                del new_lines[0]

            for line in new_lines:
                final_geometries.append(line)
        lines = final_geometries

    multi = ogr.Geometry(5)
    for line in lines:
        multi.AddGeometry(line)

    return multi
