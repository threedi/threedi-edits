# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 18:53:32 2021

@author: chris.kerklaan
"""
# First-party imports
import pathlib

# Local imports
from threedi_edits.gis.vector import Vector

# Globals
DATA = pathlib.Path(__file__).parent / "data" / "gis_feature.gpkg"


def test_set_feature():
    """tests if feature has been set after field is added"""
    lines = Vector(DATA, layer_name="lines")
    lines = lines.copy()
    lines.add_field("test", float)

    lines[0]["test"] = 0.1

    assert lines[0]["test"] == 0.1
    lines.close()
