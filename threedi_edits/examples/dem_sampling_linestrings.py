# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:14:30 2021

@author: chris.kerklaan
"""
import os
import numpy as np
from threedi_edits import Vector, Raster, Progress

os.chdir(r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing")
levee_path = "levees_19gz2_types_2.shp"
dem_path = r"\\utr-3fs-01.nens.local\WorkDir\J_vLieshout\HHNK_AHN_basiskaarten\AHN4_Interpolated\ahn4_interpolated\19gz2.tif"
levees = Vector(levee_path)
dem = Raster(dem_path)
output = levees.copy(shell=True)
output.add_field("height", float)

output_merged = levees.copy(shell=True)
output_merged.add_field("height", float)
for levee in Progress(levees, "Sampling levees"):
    lines = levee.geometry.perpendicular_lines(distance=2, perpendicular_length=5)

    heights = []
    for line in lines:
        try:
            height = dem.read(line)
        except Exception:
            height = 0

        output.add(geometry=line, items=levee.items, height=height)
        heights.extend(list(height))
    average = np.nanmean(heights)
    output_merged.add(geometry=levee.geometry, items=levee.items, height=average)

output.write("levee_lines_example2.shp")
output_merged.write("levee_lines_example_merged2.shp")
