# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 15:08:00 2021

@author: chris.kerklaan
"""
# Import the Vector and Raster classes from threedi_edits
from threedi_raster_edtis import Vector, Raster
import numpy as np

# path to data
path_to_raster = "path/to/raster.tif"
path_to_vector = "path/to/vector.shp"
output_path = "path/to/output.shp"

# load raster and vector data
vector = Vector(path_to_vector)
raster = Raster(path_to_raster)

# copy the vector without the data
output = vector.copy(shell=True)

# add a field to the output
output.add_field("sum", float)

# loop over each feature in the vector
for feature in vector:

    # read the raster data
    raster_data = raster.read(feature.geometry)

    # do something with the data
    data_sum = np.nansum(raster_data)

    # write the output
    items = {"sum": data_sum}
    output.add(items=items, geometry=feature.geometry)


# Write to output_path
output.write(output_path)
