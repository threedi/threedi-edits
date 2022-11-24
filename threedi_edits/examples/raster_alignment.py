# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 10:43:59 2022

@author: chris.kerklaan
"""


import threedi_raster_edits as tre


raster_path_a = r"/path/to/raster.tif"

raster_path_b = r"/path/to/rasterb.tif"

raster_a = tre.Raster(raster_path_a)
raster_b = tre.Raster(raster_path_b)

help(raster_a.align)
raster_a.align(raster_b)
