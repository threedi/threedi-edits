# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 10:39:47 2021

@author: chris.kerklaan

Steps for sampling/reading data from vector geometries and burning it back to 
the dem. 
"""
import os
import numpy as np
from threedi_edits import Vector, Raster, Progress, RasterGroup

# dem_path = r"data/09bn2.tif"
# buildings_path = r"data/Panden_nieuw.shp"
# tegel_path = r"data/AHN_HHNK_grenzen.shp"

# os.chdir(r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing")

# # load raster
# dem = Raster(dem_path)

# # get dem tile
# tegels = Vector(tegel_path)
# tegel = tegels.filter(name="09bn2")[0]

# # get buildings within tile.
# buildings = Vector(buildings_path, create_index=True)
# buildings_within_tile = buildings.spatial_filter(tegel.geometry)

# # make output
# output = buildings.copy(shell=True)
# output.add_field("level", float)

# # find percentile, add buffer
# for building in buildings_within_tile:

#     # do a 1 meter buffer
#     geometry_buffered = building.geometry.buffer(1)

#     # difference with the original to create a hull
#     difference = geometry_buffered.difference(building.geometry)

#     # read the data
#     data = dem.read(difference)

#     # find the percentile
#     percentile = np.percentile(data, 75)

#     # percentile + 5 cm
#     floorheight = percentile + 0.05

#     # write to output with new fields
#     items = building.items
#     items["level"] = floorheight
#     output.add(items=items, geometry=building.geometry)

# output.write("building_output.shp")

# # find out how to find functions
# dir(dem)

# # find arguments
# help(dem.push_vector)

# # push the data onto the raster
# dem_floorheight = dem.push_vector(output, field="level")
# dem_floorheight.write("dem_level.tif")


# # Use the above data as a function a loop over multiple tiles
# def burn_floorlevels_in_tile(dem_tile: Raster, buildings: Vector, output_raster_path):
#     # get buildings within tile. Slightly different then the above method
#     # we are now directory getting the extent from the raster with raster.extent_geometry
#     buildings_within_tile = buildings.spatial_filter(dem_tile.extent_geometry)

#     # make output
#     output = buildings.copy(shell=True)
#     output.add_field("level", float)

#     # find percentile, add buffer
#     for building in buildings_within_tile:

#         # do a 1 meter buffer
#         geometry_buffered = building.geometry.buffer(1)

#         # difference with the original to create a hull
#         difference = geometry_buffered.difference(building.geometry)

#         # read the data
#         data = dem_tile.read(difference)

#         # find the percentile
#         percentile = np.percentile(data, 75)

#         # percentile + 5 cm
#         floorheight = percentile + 0.05

#         # write to output with new fields
#         items = building.items
#         items["level"] = floorheight
#         output.add(items=items, geometry=building.geometry)

#     dem_floorheight = dem_tile.push_vector(output, field="level")
#     dem_floorheight.write(output_raster_path)


# # Now we are using this function to automate our processes
# dem_directory = r"\\utr-3fs-01.nens.local\WorkDir\J_vLieshout\HHNK_AHN_basiskaarten\1_AHN4_Interpolated\ahn4"

# # We can only show our progress using the Progress class as is show below
# for dem_file in Progress(os.listdir(dem_directory), "Generating floorlevels"):

#     # skip the file if it does not end with '.tif'
#     if not dem_file.endswith(".tif"):
#         continue

#     # load data
#     dem_tile = Raster(dem_directory + "/" + dem_file)
#     buildings = Vector(
#         r"\\utr-3fs-01.nens.local\WorkDir\J_vLieshout\HHNK_AHN_basiskaarten\Panden\panden_met_index.shp"
#     )

#     # make output path
#     output_raster_path = (
#         r"\\utr-3fs-01.nens.local\WorkDir\J_vLieshout\HHNK_AHN_basiskaarten\1_AHN4_Interpolated\ahn4_vloerpijl/"
#         + dem_file
#     )

#     # start processing
#     burn_floorlevels_in_tile(dem_tile, buildings, output_raster_path)


## bruggen

panden_path = r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\data\Panden_nieuw.shp"
buildings = Vector(panden_path)
dem_tiles_dir = r"K:\J_vLieshout\HHNK_AHN_basiskaarten\AHN4_DTM"

dem_tiles = [
    Raster(dem_tiles_dir + f"/{i}")
    for i in os.listdir(dem_tiles_dir)
    if i.endswith("tif")
]
dem = RasterGroup(dem_tiles)

# buffered = buildings.buffer(1, quiet=False)

# buffered = Vector( r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\data\panden_buffered.shp")
# clusters = buffered.dissolve(intersecting=True, quiet=False)
# clusters.write(r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\data\panden_clusters.shp")


buildings = Vector(
    r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\data\Panden_nieuw.shp"
)
clusters = Vector(
    r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\data\panden_clusters.shp"
)

output = buildings.copy(shell=True)
for cluster in Progress(clusters, "Getting the sample areas"):

    building_ids = cluster["diss_ids"].split("_")
    for building_id in building_ids:

        if building_id == "":
            continue
        building = buildings[int(building_id)]

        # do a 1 meter buffer
        geometry_buffered = building.geometry.buffer(1)

        # difference with the original to create a hull
        difference = geometry_buffered.difference(building.geometry)

        # # also it cannot intersect with other builindgs
        for intersecting_building_id in building_ids:
            if intersecting_building_id == building.id:
                continue
            if intersecting_building_id == "":
                continue

            intersecting_building = buildings[int(intersecting_building_id)]

            intersecting_buffer = intersecting_building.geometry.buffer(1)
            difference = difference.difference(intersecting_buffer)

        if difference.IsEmpty():
            continue
        if difference.GetGeometryType() == 3:

            output.add(geometry=difference, fid=int(building_id))

output.write(
    r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing/data/panden_sample_area2.shp"
)

buildings = Vector(
    r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\panden_sampled/Panden_nieuw.shp"
)
new = buildings.copy(shell=True)
new.add_field("vloer", float)
for building in Progress(buildings, "buildings"):
    new.add(
        items=building.items,
        vloer=building["level_75"] + 0.05,
        geometry=building.geometry,
    )


for i, dem_tile_name in enumerate(os.listdir(dem_tiles_dir)):
    if not dem_tile_name.endswith("tif"):
        continue

    output_path = rf"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\dtm_interpolated/{dem_tile_name}"
    if os.path.exists(output_path):
        print("Exists")
        continue

    print(dem_tile_name, f"{i}/{len(os.listdir(dem_tiles_dir))}")
    dem_tile = Raster(
        rf"G:\Projecten W (2021)\W0234 - Nieuwe basiskaarten 3Di voor HHNK\Gegevens\Resultaat\21-12-2021 - Oplevering\Geinterpoleerde DEM\{dem_tile_name}"
    )
    new_dem = dem_tile.push_vector(new, field="vloer")
    new_dem.write(output_path)


output = Vector.from_scratch("panden", 3, 28992)
output.add_field("fid", int)

output.add_field("gid", int)
output.add_field("bouwjaar", int)
output.add_field("identifica", str)
output.add_field("pandstatus", str)
output.add_field("vloerpeil_75", float)
output.add_field("vloerpeil_50", float)
output.add_field("vloerpeil_25", float)
output.add_field("vloerpeil_min", float)
output.add_field("vloerpeil_max", float)

for building in Progress(buildings, "Sampling panden"):

    # do a 1 meter buffer
    geometry_buffered = building.geometry.buffer(1)

    # difference with the original to create a hull
    difference = geometry_buffered.difference(building.geometry)

    # # also it cannot intersect with other builindgs
    for intersecting_building_id in building_ids:
        intersecting_building = buildings[int(intersecting_building_id)]
        if intersecting_building.id == building.id:
            continue
        intersecting_buffer = intersecting_building.geometry.buffer(1)
        difference = difference.difference(intersecting_buffer)

    # read the data
    data = dem.read(difference)

    if len(data) == 0:
        # do a 1 meter buffer
        geometry_buffered = building.geometry.buffer(3)

        # difference with the original to create a hull
        difference = geometry_buffered.difference(building.geometry)

        # read the data
        data = dem.read(difference)

    if len(data) == 0:
        data = [np.nan]

    # find the percentile
    percentile_75 = np.nanpercentile(data, 75) + 0.05
    percentile_50 = np.nanpercentile(data, 50) + 0.05
    percentile_25 = np.nanpercentile(data, 25) + 0.05
    percentile_min = np.nanmin(data) + 0.05
    percentile_max = np.nanmax(data) + 0.05

    # write to output with new fields
    items = building.items
    items["percentile_75"] = percentile_75
    items["percentile_50"] = percentile_50
    items["percentile_25"] = percentile_25
    items["percentile_min"] = percentile_min
    items["percentile_max"] = percentile_max

    output.add(items=items, geometry=building.geometry)
output.write(
    r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing\panden_sampled\panden_hhnk.gpkg"
)
