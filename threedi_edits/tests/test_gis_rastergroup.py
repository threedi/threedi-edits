# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:38:51 2021

@author: chris.kerklaan


Script for automated testing of the raster group object


Currently supported functions
    1. Check alignment
        1. Data/nodata value
        2. geotransform
        3. rows, columns
        
    

"""
# First-party imports
import pathlib

# Third-party imports
import numpy as np

# Local imports
from threedi_edits.gis.raster import Raster
from threedi_edits.gis.rastergroup import RasterGroup

# Globals
DATA = pathlib.Path(__file__).parent / "data" / "raster" / "groups"


RASTERS = []
for raster_path in DATA.glob("*.tif"):
    print(raster_path)
    raster = Raster(raster_path)
    RASTERS.append(raster.copy(shell=False, memory=False))


def test_check_alignment_geotransform():
    """Tests is alignment checking notices geotransform changes"""
    test_rasters = [raster.copy() for raster in RASTERS]
    test_rasters[1].geotransform = (96798.0, 0.5, 0.0, 417515.0, 0.0, -0.5)
    rastergroup = RasterGroup(test_rasters)

    output = rastergroup.check_alignment()
    assert len(output["errors"]) > 0


def test_check_alignment_rows_columns():
    """Tests is alignment checking notices rows/columns changes"""
    test_rasters = [raster.copy() for raster in RASTERS]
    new = Raster.from_scratch(
        "new",
        test_rasters[1].geotransform,
        1,
        1,
    )
    test_rasters.append(new)
    rastergroup = RasterGroup(test_rasters)

    output = rastergroup.check_alignment()
    assert len(output["errors"]) > 0


def test_check_alignment_data_nodata():
    """Tests is alignment checking notices data/nodata changes"""
    test_rasters = [raster.copy() for raster in RASTERS]
    array = test_rasters[2].array

    array[array > 0] = np.nan
    test_rasters[2].array = array
    test_rasters[2].name = "Test raster"
    rastergroup = RasterGroup(test_rasters)

    output = rastergroup.check_alignment()
    errors_present = len(output["errors"]) > 0
    assert errors_present

    if errors_present:
        assert "counts" in output["errors"][0]


def test_tiling():
    """tests if tiling works"""
    test_rasters = [raster.copy() for raster in RASTERS]
    group = RasterGroup(test_rasters)

    for tiles in group.tiles():
        assert tiles[0].shape == (1, 1032)
        break
