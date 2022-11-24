# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 16:04:45 2022

@author: chris.kerklaan

#TODO:
    Create similar tests as vector.py however add some speed tests.
    
Files used:
    lines.shp
    single.shp
    
"""

# First-party imports
import tempfile
import pathlib

# Third-party imports
from osgeo import ogr

# Local imports
from threedi_edits.gis.vectorbase import VectorBase, load_table

# Globals
# __file__ = "C:/Users/chris.kerklaan/Documents/Github/threedi-raster-edits/threedi_edits/tests/test_gis_vectorbase.py"
DATA = pathlib.Path(__file__).parent / "data" / "gis_vector.gpkg"


def test_load_path():
    """
    Tests if vector is correctly loaded using a path
    A table is also loaded so, it should influence the speed.
    """

    lines = VectorBase(DATA, "lines")
    assert lines.count > 0
    lines.close()


# def test_load_pg():
#     """ tests if vector is correctly loaded using a postgres ds"""
#     pg_test = Vector.from_pg(
#         host="localhost",
#         port=5432,
#         password="postgres",
#         user="postgres",
#         dbname="postgres",
#         layer_name="buurten",
#     )
#     assert pg_test.count > 0
#     pg_test.close()


def test_load_ds():
    """tests if vector is correctly loaded using an ogr datasource"""
    lines_ds = ogr.Open(str(DATA))
    lines = VectorBase.from_ds(lines_ds)
    assert lines.count > 0
    lines.close()


def test_load_table(benchmark):
    lines = VectorBase(DATA, "bgt_panden_fixed")
    layer = lines.layer
    fields = lines.fields

    threshold = 0.05  # or 100 ms or 100.000 us for 3000 features.
    benchmark(load_table, layer, fields.names, fields.types, lines.has_geometry)
    assert benchmark.stats.stats.mean < threshold
    lines.close()


def test_load_scratch():
    """tests if vector can be created from scratch"""
    multi = VectorBase.from_scratch("temp", ogr.wkbMultiPolygon, 28992)
    assert isinstance(multi, VectorBase)


def test_spatial_filtering(benchmark):
    """tests spatial filtering on witin and intersect, with and without layer filter"""

    lines = VectorBase(DATA, "lines")
    lines = VectorBase(lines._copy()[1])

    geometry = VectorBase(DATA, "geometry")[0].geometry

    fids = lines._sfilter(geometry, method="intersects")
    assert len(fids) == 21

    fids = lines._sfilter(geometry, method="intersects")
    assert len(fids) == 21

    fids = lines._sfilter(geometry, method="within")
    assert len(fids) == 14

    fids = lines._sfilter(geometry, method="within")
    assert len(fids) == 14

    threshold = 1 / 1000  # seconds to ms
    benchmark(lines._sfilter, geometry=geometry)
    #     kwargs={"geometry": geometry}, rounds=10
    # )
    assert benchmark.stats.stats.mean < threshold
    lines.close()


def test_table_filter():
    """tests on table filtering
    1. Single filter
    2. double filter
    3. Feature output
    4. fid output
    """

    lines = VectorBase(DATA, "lines")
    lines = VectorBase(lines._copy()[1])
    assert len(VectorBase(lines._filter(HOID="44466"))) == 15
    assert len(VectorBase(lines._filter(HOID="44466", test_field=2))) == 1
    lines.close()


def test_add_feature():
    """tests if a feature is correctly added by using a feature,
    a dictionary and a geometry
    """

    lines = VectorBase(DATA, "lines")
    lines = VectorBase(lines._copy()[1])
    single = VectorBase(DATA, "single")

    lines._add(single[0])
    assert lines.count == 90
    assert lines[89].items["test_field"] == 1

    for value in lines.table.values():
        assert len(value) == 90
    single.close()


def test_delete_feature():
    """tests if feature is properly deleted"""
    lines = VectorBase(DATA, "lines")
    lines = VectorBase(lines._copy()[1])
    lines.delete(0)
    assert lines.count == 88
    assert len(lines.table["fid"]) == 88


def test_copy():
    """tests if features are copied"""
    lines = VectorBase(DATA, "lines")
    copy = VectorBase(lines._copy()[1])
    assert lines.settings == copy.settings
    assert lines.count == copy.count
    assert lines[10].items == copy[10].items
    lines.close()


def test_copy_shell():
    """tests if they have no count"""
    lines = VectorBase(DATA, "lines")
    copy = VectorBase(lines._copy(shell=True)[1])
    assert lines.settings == copy.settings
    assert copy.count == 0
    lines.close()


def test_add_field():
    """tests if field is added and if one can add something"""
    lines = VectorBase(DATA, "lines")
    lines = VectorBase(lines._copy()[1])
    lines.add_field("x", float)
    assert "x" in lines.fields.items

    lines._add(items={"x": 1.0})
    assert lines[89].items["x"] == 1
    lines.close()


def test_delete_field():
    lines = VectorBase(DATA, "lines")
    lines = VectorBase(lines._copy()[1])
    lines.delete_field("test_field")
    assert "test_field" not in lines.fields.items
    assert "test_field" not in lines[80].items
    lines.close()


def test_write_shape():
    """tests if a shape can be written"""

    with tempfile.TemporaryDirectory() as folder:
        lines = VectorBase(DATA, "lines")
        lines = VectorBase(lines._copy()[1])
        lines.write(folder + "/lines_write.shp")
        lines2 = VectorBase(folder + "/lines_write.shp")
        assert lines.count == lines2.count
        assert lines[0].items == lines2[0].items
        assert lines[0].geometry.wkt == lines2[0].geometry.wkt
        lines.close()
        lines2.close()


def test_write_gpkg():
    """tests if a geopackage can be written"""
    with tempfile.TemporaryDirectory() as folder:
        lines = VectorBase(DATA, "lines")
        lines = VectorBase(lines._copy()[1])
        lines.write(folder + "/lines_write.gpkg")
        lines2 = VectorBase(folder + "/lines_write.gpkg")
        assert lines.count == lines2.count
        assert lines[0].items == lines2[1].items
        assert lines[0].geometry.wkt == lines2[1].geometry.wkt
        lines.close()
        lines2.close()
