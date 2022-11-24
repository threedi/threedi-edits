# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 13:52:55 2019

@author: chris.kerklaan
"""
import importlib


def get_version(version):
    v = importlib.import_module(f"threedi_edits.threedi.constants.v{version}")
    return v
