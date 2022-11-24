# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:54:38 2021

@author: chris.kerklaan

Functions to test:
#Done:
    # functions
    - rasterize
    - clip 
    - dissolve
    - difference
    - reproject
    - multipart2singlepart
    - centroid
    - simplify

"""
# First-party imports
import sys
import pathlib


# Third-party imports
import numpy as np
from osgeo import ogr

# Local imports
from threedi_edits.gis.vector import Vector
from threedi_edits import Raster

# Globals
# __file__ = "C:/Users/chris.kerklaan/Documents/Github/threedi-raster-edits/threedi_edits/tests/test_gis_vector.py"
DATA = pathlib.Path(__file__).parent / "data" / "gis_vector.gpkg"
RASTER_FOLDER = pathlib.Path(__file__).parent / "data" / "raster"


def test_fix(benchmark):
    """tests if multisurfaces are fixed"""
    bgt = Vector(DATA, layer_name="bgt_panden")
    bgt = bgt.copy()
    fixed = bgt.fix(geometry_type=ogr.wkbMultiPolygon)
    assert fixed.geometry_type == ogr.wkbMultiPolygon
    assert fixed[1].geometry.type == ogr.wkbMultiPolygon

    # benchmark
    # threshold = 0.3
    # benchmark(bgt.fix, geometry_type=ogr.wkbMultiPolygon)
    # assert benchmark.stats.stats.mean < threshold


def test_rasterize_ones():
    """tests rasterization based on counts of array using field"""
    lines = Vector(DATA, layer_name="lines")
    geotransform = (117366.7415, 1.0, 0.0, 459774.6888, 0.0, -1.0)
    rows = 151
    columns = 280
    ds = lines.rasterize(rows, columns, geotransform)
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    array[array == -9999] = np.nan
    assert np.nansum(array) == 1757.0
    lines.close()


def test_rasterize_fields():
    """tests rasterization based on counts of array using field"""

    lines = Vector(DATA, layer_name="lines")
    geotransform = (117366.7415, 1.0, 0.0, 459774.6888, 0.0, -1.0)
    rows = 151
    columns = 280
    ds = lines.rasterize(rows, columns, geotransform, field="test_field")
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    array[array == -9999] = np.nan
    assert np.nansum(array) == 1778.0
    lines.close()


def test_buffer(benchmark):
    """tests if buffer algorithm is performed correctly"""
    lines = Vector(DATA, layer_name="lines")
    buffered = lines.buffer(10, quiet=True)
    assert buffered[1].geometry.area > 0
    lines.close()
    buffered.close()

    # threshold = 0.05
    # benchmark.pedantic(lines.buffer, kwargs={"buffer_size":10}, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_simplify(benchmark):
    """tests if simplify algorithm is performed correctly"""
    shape = Vector(DATA, layer_name="strange_shape")
    simplified = shape.simplify(10)
    assert sys.getsizeof(simplified[1].geometry.wkb) < sys.getsizeof(
        shape[0].geometry.wkb
    )

    # threshold = 0.05
    # benchmark.pedantic(shape.simplify, 10, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_centroid(benchmark):
    """tests if centroid algorithm is performed correctly"""
    lines = Vector(DATA, layer_name="lines")
    centroids = lines.centroid()
    original = lines[0].geometry.Centroid().ExportToWkt()
    new = centroids[1].geometry.centroid.wkt
    assert original == new


def test_multipart_to_singlepart(benchmark):
    """tests if multipart to singlepart on geometry type and counts"""
    multi = Vector(DATA, layer_name="multi")
    singles = multi.to_single()
    assert singles.geometry_type == 3
    assert singles.count == 2
    assert singles[1].geometry.type == 3

    # threshold = 0.05
    # benchmark.pedantic(multi.to_single, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_reproject(benchmark):
    """tests if reprojects has a new bbox, epsg is changed and geometry has changed"""
    lines = Vector(DATA, layer_name="lines")
    new_lines = lines.reproject(4326)
    assert max(new_lines.extent) < 60
    assert new_lines[1].geometry.points[0][0] < 60

    # threshold = 0.05
    # benchmark.pedantic(lines.reproject, 4326, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_clip(benchmark):
    """tests if a clip works by count and feature geometry"""
    lines = Vector(DATA, layer_name="lines")
    mask = Vector(DATA, layer_name="geometry")
    clipped = lines.clip(mask)
    assert clipped.count == 21

    feat_14 = "LINESTRING (117556.454477141 459720.712475407,117576.403681566 459722.13699034)"
    assert clipped[14].geometry.wkt == feat_14

    # threshold = 0.05
    # benchmark.pedantic(lines.clip, mask, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_difference(benchmark):
    """tests if a difference works by count and feature geometry"""
    lines = Vector(DATA, layer_name="lines")
    mask = Vector(DATA, layer_name="geometry")
    differenced = lines.difference(mask)
    assert differenced.count == 75
    feat_51 = "LINESTRING (117575.529616571 459712.808025775,117577.996659787 459713.235393432)"
    assert differenced[52].geometry.wkt == feat_51

    # threshold = 0.05
    # benchmark.pedantic(lines.difference, mask, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_polygon_to_lines():
    """tests polygon to lines"""
    geometry = Vector(DATA, layer_name="geometry")
    lines = geometry.polygon_to_lines()
    assert lines.geometry_type == 2
    assert lines[1].geometry.type == 2


def test_dissolve(benchmark):
    """tests if a difference works by count and wkbsize"""
    lines = Vector(DATA, layer_name="lines")
    dissolved = lines.dissolve()
    assert dissolved.count == 1
    assert dissolved[1].geometry.wkb_size == 3617
    assert dissolved[1].geometry.type == ogr.wkbMultiLineString
    assert dissolved.geometry_type == ogr.wkbMultiLineString

    # threshold = 0.05
    # benchmark.pedantic(lines.dissolve, iterations=10)
    # assert benchmark.stats.stats.mean < threshold


def test_dissolve_field():
    """tests if a difference works by count and wkbsize"""
    lines = Vector(DATA, layer_name="lines")
    dissolved = lines.dissolve(field="HOID")
    assert dissolved.count == 9


def test_integrated_delete_dissolve():
    """tests if integrated works after polygonizing, looping, deleting and getting count
    problem:
        python crashes when this used to be done.

    """
    lines_raster = Raster(RASTER_FOLDER / "lines.tif")
    lines_vector_path = lines_raster.polygonize(quiet=True)
    lines_vector = Vector(lines_vector_path)
    lines_vector = lines_vector.copy()
    count = lines_vector.count

    for i, feature in enumerate(lines_vector):
        if i == 2:
            lines_vector.delete(feature.fid)

    assert lines_vector.count == count - 1

    dissolved = lines_vector.dissolve()
    single_dissolved = dissolved.to_single()
    assert single_dissolved.count == count - 1


def test_voronoi():
    """tests if voronoi can be done"""
    points = Vector(DATA, layer_name="points")
    voronoi = points.voronoi()
    assert voronoi.count == points.count


# def test_spatial_clusters():
#     """  tests if spatial clusters works"""
#     panden = Vector(DATA, layer_name="bgt_panden_fixed")
#     clusters = panden.spatial_clusters()
