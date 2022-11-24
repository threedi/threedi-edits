# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:37:33 2022

@author: chris

A simple version of editing a threedimodel in python
"""

# variables
import threedi_edits as tre

sqlite_path = r"C:\Users\chris.kerklaan\Documents\Projecten\overijssels_kanaal\processing\model\winter.sqlite"
point_path = r"C:\Users\chris.kerklaan\Documents\Projecten\overijssels_kanaal\processing\point.shp"
output = r"C:\Users\chris.kerklaan\Documents\Projecten\overijssels_kanaal\processing"

# First show all logging

# Load in memory so it is editable.
model = tre.ThreediEdits(sqlite_path, mode="memory")

# resettings all friction values of the cross section
cross_section_locations = model.cross_section_locations
for cross_section_location in cross_section_locations:
    cross_section_location["friction_value"] = 0.1

# show the edited part values
print(cross_section_locations[cross_section_locations.fids[0]])

# add an vertice to the channel
# get the data
points = tre.Vector(point_path)
point_geometry = points[0].geometry

# get the right channel
channels = model.channels

# get the channel by filter
channel = channels.filter(code="18049")[0]

# get the channel by id
channel = channels[26064]

# see what you can do with the geometry
dir(channel.geometry)
help(channel.geometry.add_vertice)

# show differences in linestring:
# Keep snapping size as low as possible
print(channel.geometry.add_vertice(point_geometry))
print(channel.geometry)

channel.geometry = channel.geometry.add_vertice(point_geometry)

# Templates are used to add a row.
templates = tre.Templates()
cross_section = templates.cross_section

# You can look up some of the documentation by printing the template
print(cross_section)

# Edit the template
cross_section["height"] = 5

# print again to show
print(cross_section)

# add to the cross sections
cross_section_defintions = model.cross_section_definitions

# Note this gives an error, because code cannot be None
cross_section_defintions.add(cross_section)

# Edit all other things
cross_section["code"] = "x"
cross_section["shape"] = 1
cross_section["width"] = 1

# Add
cross_section_defintions.add(cross_section)

# You can write and check the whole model at once
model.write(output + "/model2.sqlite", rasters=True)

# or write only the cross sections for a faster look
cross_section_locations.write(output + "/csl.gpkg")
channels.write(output + "/channel.gpkg")

raster = model.rasters

raster.check_properties()

dem = raster.dem
infiltratie = raster.infiltration

output = dem.copy()
for tile in dem:
    array = tile.array * 1
    output.array = (tile.location, array)

raster.dem = output
