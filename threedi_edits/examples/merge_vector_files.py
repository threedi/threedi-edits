# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:20:45 2021

@author: chris.kerklaan

Example shows how to merge multiple vector files

"""
import os
from threedi_edits import Vector, Progress

os.chdir(r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing")

# define paths
dikes_path = "data/waterkering_werkbestand.shp"
roads_path = "data/wegen.gpkg"
rail_path = "data/spoorwegen.gpkg"

# load data
dikes = Vector(dikes_path)
roads = Vector(roads_path)
rails = Vector(rail_path)

# create output -- help
help(Vector.from_scratch)

# create a vector with the name levees, a linstring and projection of the Netherlands
output = Vector.from_scratch("levees", 2, 28992)

# add a field to track which type of data is added
output.add_field("type", str)

# add dike keep geometry and 'code', 'category' and 'name'
# first we add fields to the output
output.add_field("waterkering_code", str)
output.add_field("waterkering_naam", str)
output.add_field("waterkering_categorie", str)
output.add_field("weg_fysiekvoorkomen", str)
output.add_field("weg_typeweg", str)
output.add_field("weg_id", int)
output.add_field("weg_naam", str)

output.add_field("spoor_fuuid", str)
output.add_field("spoor_geocode_naam", str)


# show progress as well
for dike in Progress(dikes, "Adding dikes"):
    # extract data
    dike_code = dike["CODE"]
    dike_name = dike["NAAM"]
    dike_category = dike["CATEGORIE"]

    # transform to new names
    if dike_category == 1:
        type = "kering-primair"
    else:
        type = "kering-secundair"
    items = {
        "waterkering_code": dike_code,
        "waterkering_naam": dike_name,
        "waterkering_categorie": dike_category,
        "type": type,
    }
    # 3D geom
    dike_geometry = dike.geometry.fix()
    # add to output
    output.add(items=items, geometry=dike_geometry)


# add roads keep geometry and 'fysiekvoorkomen', 'typeweg', 'id', 'name'
# first we add fields to the output

# show progress as well
for road in Progress(roads, "Adding roads"):
    # extract data
    road_id = road["id"]
    road_name = road["Naam"]
    road_fysiekvoorkomen = road["fysiekvoorkomen"]
    road_type = road["typeweg"]

    # transform to new names
    items = {
        "weg_fysiekvoorkomen": road_fysiekvoorkomen,
        "weg_typeweg": road_type,
        "weg_id": road_id,
        "weg_naam": road_name,
        "type": "weg",
    }

    # add to output
    output.add(items=items, geometry=road.geometry)


# show progress as well
for rail in Progress(rails, "Adding rails"):

    # extract data
    rail_name = rail["geocode_naam"]
    rail_id = rail["fuuid"]

    # transform to new names
    items = {"spoor_fuuid": rail_id, "spoor_geocode_naam": rail_name, "type": "spoor"}

    # add to output
    output.add(items=items, geometry=rail.geometry)

# write as a geopackeg to keep the names
output.write("levees/all_levees.gpkg")
