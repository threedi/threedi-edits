# -*- coding: utf-8 -*-
"""
@author: chris.kerklaan

#TODO functions:
    1. Convert geometry
    2. fix geometry

Currently supported functions
    1. buffer
    2. difference
    3. clip
    4. reproject
    5. to_single

"""

# First-party imports
import pathlib

# Local imports
from threedi_edits.gis.geometry import Geometry
from threedi_edits import Vector

# Globals
DATA = pathlib.Path(__file__).parent / "data" / "gis_geometry.gpkg"
polygon = Vector(DATA, layer_name="geometry")[0].geometry
multiline = Vector(DATA, layer_name="lines").dissolve()[1].geometry


def test_buffer():
    """tests if buffer is correctly performed"""
    geom = Geometry(polygon)
    original_area = geom.area
    buffered = geom.buffer(1)
    assert original_area < buffered.area


def test_difference():
    """tests if difference length is smaller than original length"""
    line = Geometry(multiline)
    differenced = line.difference(polygon)
    assert line.length > differenced.length


def test_clip():
    """tests if clip is correctly performed with correct geom_type"""
    lines = Geometry(multiline)
    clipped = lines.clip(polygon)[0]
    assert clipped.type == 5
    assert clipped.length < lines.length


def test_reproject():
    """tests if geometry is reprojected and correct epsg"""
    geom = Geometry(polygon)
    reprojected = geom.reproject(4326)
    assert reprojected.area < geom.area
    assert reprojected.epsg == 4326


# def test_dissolve():
#     """ tests if ouput of dissolve is 1"""
#     lines = Geometry(multiline)
#     lines.to_single()
#     dissolved = lines.dissolve()


def test_to_single():
    """tests if to single produced multiple geometries"""
    lines = Geometry(multiline)
    single_lines = lines.to_single()
    assert len(single_lines) > 1
