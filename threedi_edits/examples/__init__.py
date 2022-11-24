# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 14:39:17 2019

@author: chris.kerklaan
"""
import pathlib
import shutil
from threedi_edits.utils.project import Files
from threedi_edits.utils.project import Modules

files = Files(__file__)
modules = Modules(__name__)


def get_examples(destination_dir):
    current_directory = str(pathlib.Path(__file__).parent)
    shutil.copytree(current_directory, destination_dir)
