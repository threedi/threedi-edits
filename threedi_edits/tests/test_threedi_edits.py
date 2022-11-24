# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 16:02:59 2021

@author: chris

#TODO

#DONE

    # input
    - Create empty sqlite
    - Loading in read mode.
    - Copy a model.
    - Create a grid with a path.
    - Create a grid
    
    # adding 
    - Add a node with items.
    - Add a node with template.
    - add a node with feature
    
    # deleting
    - delete a node
    
    # changing
    - Setting field from a feature.
    - Setting field from a vector.
    - Setting a raster array.
    - Setting a whole new raster.
    - None correcting for memory.
    
    # output
    - Writing from memory.
    
"""
# First-party imports
import pathlib
import tempfile

# Local imports
from threedi_edits.gis.point import Point
from threedi_edits.threedi.edits import ThreediEdits
from threedi_edits.globals import SUPPORTED_THREEDI_VERSIONS
from threedi_edits.threedi.constants import get_version

V = get_version(SUPPORTED_THREEDI_VERSIONS[0])


# Globals
# __file__ = "C:/Users/chris.kerklaan/Documents/Github/threedi-raster-edits/threedi_edits/tests/test_threedi_edits.py"

Templates = V.Templates
DATA = pathlib.Path(__file__).parent / "data" / "threedi_edits"
DEM_PATH = DATA / "dem.tif"

BWN_PATH = DATA / "bwn_test.sqlite"
KLONDIKE_PATH = DATA / "empty_klondike.sqlite"
BWN_SCHEMA_MEMORY = ThreediEdits(BWN_PATH, "memory")


def test_create_empty_model():
    """tests if an empty model can be created"""
    with tempfile.TemporaryDirectory() as folder:

        model = ThreediEdits.from_scratch()
        model.write(folder + "/empty.sqlite")
        model = ThreediEdits(folder + "/empty.sqlite", "memory")

    assert type(model) == ThreediEdits


def test_load_read_mode():
    BWN_SCHEMA = ThreediEdits(BWN_PATH, "read")
    assert BWN_SCHEMA.nodes[15]["code"] == "15"


def test_copy():
    BWN_SCHEMA_COPY = BWN_SCHEMA_MEMORY.copy()
    assert BWN_SCHEMA_COPY.nodes[15]["code"] == "15"


def test_setting_raster():
    """we tests if rasters can be edited and put back in the object"""
    BWN_SCHEMA_MEMORY.rasters  # initialize rasters
    dem = BWN_SCHEMA_MEMORY.dem
    new_dem = dem.copy()
    array = new_dem.array
    array[0, 0] = 2
    new_dem.array = array
    BWN_SCHEMA_MEMORY.dem = new_dem

    assert BWN_SCHEMA_MEMORY.dem.array[0, 0] == 2


def test_create_grid_path():

    BWN_SCHEMA = ThreediEdits(BWN_PATH)
    grid = BWN_SCHEMA.grid()

    assert len(grid["cells"]) == 422


def test_create_grid_memory():
    grid = BWN_SCHEMA_MEMORY.grid()
    assert len(grid["cells"]) == 422


def test_add_delete_node_template():
    templates = Templates()
    node = templates.node
    node["initial_waterlevel"] = 1
    node["the_geom"] = Point.from_point((1, 1))
    new_fid = BWN_SCHEMA_MEMORY.nodes.add(node)

    assert new_fid == 484
    for key in BWN_SCHEMA_MEMORY.nodes.table:
        assert len(BWN_SCHEMA_MEMORY.nodes.table[key]) == 73

    BWN_SCHEMA_MEMORY.nodes.delete(new_fid)

    # deletion in table
    assert BWN_SCHEMA_MEMORY.nodes.fids[-1] == 483
    for key in BWN_SCHEMA_MEMORY.nodes.table:
        assert len(BWN_SCHEMA_MEMORY.nodes.table[key]) == 72


def test_add_items_geometry():
    items = {"initial_waterlevel": 0.1, "storage_are": 0.64, "code": "x"}
    geometry = Point.from_point((1, 1))

    new_fid = BWN_SCHEMA_MEMORY.nodes.add(items=items, geometry=geometry)

    assert new_fid == 484
    for key in BWN_SCHEMA_MEMORY.nodes.table:
        assert len(BWN_SCHEMA_MEMORY.nodes.table[key]) == 73

    BWN_SCHEMA_MEMORY.nodes.delete(new_fid)


def test_add_feature():

    feature = BWN_SCHEMA_MEMORY.nodes[BWN_SCHEMA_MEMORY.nodes.fids[0]]
    new_fid = BWN_SCHEMA_MEMORY.nodes.add(feature)
    assert new_fid == 484
    for key in BWN_SCHEMA_MEMORY.nodes.table:
        assert len(BWN_SCHEMA_MEMORY.nodes.table[key]) == 73

    BWN_SCHEMA_MEMORY.nodes.delete(new_fid)


def test_type_correcting_none():

    # None float
    assert BWN_SCHEMA_MEMORY.global_settings[2]["table_step_size_1d"] is None

    # None int
    assert (
        BWN_SCHEMA_MEMORY.global_settings[1]["simple_infiltration_settings_id"] is None
    )

    # None string
    assert BWN_SCHEMA_MEMORY.global_settings[1]["start_time"] is None


def test_setting_field_from_feature():

    nodes = BWN_SCHEMA_MEMORY.nodes
    node = nodes[15]
    node["code"] = 25

    # first an internal check
    assert BWN_SCHEMA_MEMORY.nodes[15]["code"] == str(25)

    BWN_SCHEMA_MEMORY.nodes[15]["code"] = str(15)


def test_setting_field_from_vector():
    """within this test we retrieve edit a feature and try to put it back"""

    BWN_SCHEMA_MEMORY.nodes[15]["code"] = str(25)

    nodes = BWN_SCHEMA_MEMORY.nodes
    new_nodes = nodes.copy()
    new_nodes[15]["code"] = 26

    # now they should be loose, nothing happens.
    assert BWN_SCHEMA_MEMORY.nodes[15]["code"] == str(25)

    # reset so somthing should happen.
    BWN_SCHEMA_MEMORY.nodes = new_nodes
    assert BWN_SCHEMA_MEMORY.nodes[15]["code"] == str(26)

    BWN_SCHEMA_MEMORY.nodes[15]["code"] = str(15)


def test_setting_raster_array():

    BWN_SCHEMA_MEMORY.rasters

    dem = BWN_SCHEMA_MEMORY.dem

    array = dem.array
    array[0, 0] = 1
    dem.array = array

    assert BWN_SCHEMA_MEMORY.dem.array[0, 0] == 1


def test_write_from_memory():
    """memory causes losses in type checking, so we have to readjust
    and check if everything is written correctly.
    """

    with tempfile.TemporaryDirectory() as tempdir:

        path = str(tempdir) + "/bwn.sqlite"
        BWN_SCHEMA_MEMORY.write(path)

        # check if all fields are correct
        BWN_READ = ThreediEdits(path, mode="read")
        for table in BWN_READ:
            mem_table = BWN_SCHEMA_MEMORY[table.name]
            for feature in table:
                if feature in mem_table:
                    mem_feature = mem_table[feature.id]
                    assert mem_feature.items == feature.items
        table = None
        feature = None
        BWN_READ = None


# def test_node_view():
#     """tests if node view is working"""
#     model = ThreediEdits(BWN_PATH)
#     view = model.node_view
#     assert view[483]["v2_weir"][0] == 10171


# def test_node_delete():
#     model = ThreediEdits(BWN_PATH, "memory")
#     model.delete_node(483)
#     assert 483 not in model.nodes.fids
