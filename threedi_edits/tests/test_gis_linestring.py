# -*- coding: utf-8 -*-
"""
@author: chris.kerklaan

#TODO:

#Done:
    1. from_points
    2. __iter__
    3. get_points
    4. reversed
"""
# First-party imports
import pathlib

# Local imports
from threedi_edits.gis.linestring import LineString
from threedi_edits import Vector
from threedi_edits.gis.point import Point

# Globals
DATA = pathlib.Path(__file__).parent / "data" / "gis_linestring.gpkg"
line = Vector(DATA, layer_name="line")[0].geometry.ogr
split_point = Vector(DATA, layer_name="split_point")[0].geometry.ogr


def test_from_points():
    """tests if a geometry can be created from points"""
    points = [(1116651.4, 637392.6), (1188804.1, 652655.7)]
    line = LineString.from_points(points)
    assert type(line) == LineString


def test__iter__():
    """tests if ouput if iterations has count and is unequal"""
    line_geom = LineString(line)
    lines = [l for l in line_geom]
    assert len(lines) == 2
    assert lines[0].wkt != lines[1].wkt


def test_start_point():
    """tests if start point has original value"""
    line_geom = LineString(line)
    assert line_geom.start_point.point == (
        117565.43122004971,
        459721.3204805324,
    )


def test_end_point():
    """tests if start point has original value"""
    line_geom = LineString(line)
    assert line_geom.end_point.point == (
        117486.66500118251,
        459702.01357949915,
    )


def test_reversed():
    line_geom = LineString(line)
    reverse = line_geom.reversed()

    assert reverse.start_point.point == (
        117486.66500118251,
        459702.01357949915,
    )
    assert reverse.end_point.point == (117565.43122004971, 459721.3204805324)


def test_points_on_line_distance():
    line_geom = LineString(line)
    points = line_geom.points_on_line(10)
    assert int(points[3].Distance(points[4])) == 10


def test_points_on_line_wo_start_end():
    line_geom = LineString(line)
    points = line_geom.points_on_line(10, start=False, end=False)
    assert len(points) == 8


def test_points_on_line_w_vertices():
    line_geom = LineString(line)
    points = line_geom.points_on_line(10, vertices=True)
    assert len(points) == 11


def test_perendicular_lines():
    """tests if the lines have correct line length"""
    line_geom = LineString(line)
    perps = line_geom.perpendicular_lines(10, 1)

    assert int(perps[0].length) == 2
    assert len(perps) == 10


def test_add_vertice():
    """tests if extra vertice is added at point"""
    line_geom = LineString(line)
    extra = line_geom.add_vertice(Point(split_point))
    assert len(extra.points) == 4
