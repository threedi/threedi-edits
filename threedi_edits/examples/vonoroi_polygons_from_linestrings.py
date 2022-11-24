# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 14:08:57 2021

@author: chris.kerklaan
"""
import os
import sys

sys.path.append(r"C:\Users\chris.kerklaan\Documents\Github\threedi-raster-edits")
from threedi_edits import Vector

os.chdir(r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk")


def create(hydrolines, hydroparts, output_path):
    # first we write processing data
    # hydrolines.write("processing/tile/hydrolines.shp")
    # hydroparts.write("processing/tile/hydroparts.shp")

    # first buffer the lines
    hydrolines_buffered = hydrolines.buffer(1, quiet=False)
    # hydrolines_buffered.write("processing/tile/hydrolines_buffered.shp")

    # Then only select hydroparts which intersect the hydrolines
    hydroparts_filtered = hydroparts.spatial_filter(
        hydrolines, quiet=False, return_vector=True
    )
    # hydroparts_filtered.write("processing/tile/hydroparts_filtered.shp")

    # We merge to files
    merged = Vector.from_scratch("merged", 3, 28992)
    for hydroline in hydrolines_buffered:
        merged.add(hydroline)

    for hydropart in hydroparts_filtered:
        merged.add(hydropart)

    fixed = merged.fix(quiet=False)
    single = fixed.to_single(quiet=False)
    dissolved = single.dissolve(quiet=False)
    # dissolved.write("processing/tile/dissolved.shp")

    # We have to do vonoroi (influence area) of a linestring.
    # hence we make to points on a line, 10 cm from the start
    # and 10 cm from the end of the linestring to make sure
    # intersection work well.
    # Finally add some points halfway to ensure proper splits between adjacent
    # lines.
    voronoi_points = Vector.from_scratch("points", 1, 28992)
    voronoi_points.add_field("line_id", int)
    for hydroline in hydrolines:
        points_to_add = []
        for point in hydroline.geometry.points_on_line(
            custom_interval=[0.25, hydroline.geometry.length - 0.25],
            start=False,
            end=False,
        ):
            points_to_add.append(point)

        for point in hydroline.geometry.points_on_line(
            interval=5,
            start=False,
            end=False,
        ):
            points_to_add.append(point)

        for point in points_to_add:
            voronoi_points.add(geometry=point, items={"line_id": hydroline.fid})

    # voronoi_points.write("processing/tile/vonoroi_points.shp")

    # Create vonoroi polygons and dissolve the result
    # due to the fact that there are multiple points per line.
    vonoroi = voronoi_points.vonoroi()
    vonoroi = vonoroi.dissolve(field="line_id", quiet=False)
    # vonoroi.write("processing/tile/vonoroi4.shp")

    # Create the final vector by clipping it with the dissolved hydroparts.
    # We have to fix due to a 3D geometry, annoying, cannot be handled.
    singles = dissolved.to_single(quiet=False)
    singles_fixed = singles.fix(quiet=False)

    output = Vector.from_scratch("hydroparts", 3, 28992)
    output.add_field("line_id", int)

    for shape in vonoroi:
        for single in singles_fixed.spatial_filter(shape.geometry):
            clipped = single.geometry.clip(shape.geometry)
            if not clipped:
                print("Skipping shape", shape.id)
                continue
            output.add(geometry=clipped[0], items=shape.items)

    output.write(output_path, index=False)


hydrolines_path = "processing/data/hydroobjecten.shp"
hydroparts_path = "processing/data/bgt_waterdelen2.shp"
tiles_path = "processing/data/ahn_hhnk_grenzen_minor.shp"
tiles_major_path = "processing/data/AHN_HHNK_grenzen.shp"


hydrolines = Vector(hydrolines_path, create_index=True)
hydroparts = Vector(hydroparts_path, create_index=True)
tiles = Vector(tiles_path)
tiles_reversed = reversed(tiles)

for i, tile in enumerate(tiles_reversed):

    print(i, len(tiles), tile["name"])
    path = f"processing/tiles_minor/waterdeel_{tile['name']}.shp"
    linepath = f"processing/tiles_minor/hydrolijn_{tile['name']}.shp"
    if os.path.exists(path):
        print(path, "exists")
        continue

    if os.path.exists(linepath):
        tile_hydrolines = Vector(linepath)
    else:
        tile_hydrolines = hydrolines.spatial_filter(
            tile.geometry, return_vector=True, quiet=False, method="Extent"
        )

    if len(tile_hydrolines) == 0:
        if not os.path.exists(linepath):
            tile_hydrolines.write(linepath, overwrite=True, index=False)
        tile_hydrolines.write(path, overwrite=True, index=False)  # as a dummy
        continue

    tile_hydroparts = hydroparts.spatial_filter(
        tile.geometry, return_vector=True, quiet=False, method="Extent"
    )

    create(tile_hydrolines, tile_hydroparts, path)
    # tile_hydrolines.write(linepath, index=False)


# # # minor to major
# minor_tiles = Vector(tiles_path, create_table=True)
# major_tiles = Vector(tiles_major_path)

# for tile in major_tiles:
#     print(tile['name'])
#     waterdeel_path = f"processing/tiles/waterdeel_{tile['name']}.shp"
#     if os.path.exists(waterdeel_path):
#         print(f"{tile['name']} exists")
#         continue

#     skip = False
#     minors = []
#     for minor_tile in minor_tiles:
#         if tile['name'] in minor_tile['name']:
#             path = f"processing/tiles_minor/waterdeel_{minor_tile['name']}.shp"
#             if not os.path.exists(path):
#                 skip=True
#                 break
#             else:
#                 minors.append(minor_tile )

#     if skip == True:
#         print(f"Minor tile does not exist, skipping")
#         continue

#     output = Vector.from_scratch('merged', 3, 28992)
#     output.add_field("line_id", float)
#     for minor_tile in minors:
#         print("Adding", minor_tile['name'])
#         path = f"processing/tiles_minor/waterdeel_{minor_tile['name']}.shp"
#         minor = Vector(path)
#         minor = minor.fix(quiet=False)

#         minor_tile_vector = Vector.from_scratch("clip", 3, 28992)
#         minor_tile_vector.add(geometry=minor_tile.geometry.Buffer(0))
#         clipped = minor.clip(minor_tile_vector, quiet=False)
#         for feature in clipped:
#             output.add(feature)

#     print(len(output))
#     singles= output.to_single(quiet=False)
#     dissolved = singles.dissolve(field="line_id", quiet=False)

#     dissolved.write(waterdeel_path, index=False)


# # write hydrolines
# tiles_major_path = "processing/data/AHN_HHNK_grenzen.shp"
# hydrolines_path = "processing/data/hydroobjecten.shp"
# tiles = Vector(tiles_major_path)
# hydrolines = Vector(hydrolines_path)
# for i, tile in enumerate(tiles):
#     print(i, tile['name'])
#     linepath = f"processing/tiles/hydrolijn_{tile['name']}.shp"
#     tile_vector = Vector.from_scratch("clip", 3, 28992)
#     tile_vector.add(geometry=tile.geometry.Buffer(0))
#     hydrolines.set_spatial_filter(tile.geometry.Buffer(0))
#     tile_hydrolines = hydrolines.clip(tile_vector, quiet=False)
#     hydrolines.set_spatial_filter(None)

#     tile_hydrolines.write(linepath, index=False)

# merge all tiles into one

tiles = Vector(tiles_path)

output = Vector.from_scratch("all", 3, 28992)
output.add_field("line_id", int)

for tile in tiles:
    print(tile["name"])
    tile_vector = Vector.from_scratch("clip", 3, 28992)
    tile_vector.add(geometry=tile.geometry.Buffer(0))

    path = f"processing/tiles_minor/waterdeel_{tile['name']}.shp"
    waterdelen = Vector(path)
    waterdelen = waterdelen.fix(quiet=False)
    tile_hydroparts = waterdelen.clip(tile_vector, quiet=False, use_ogr=False)
    for waterdeel in tile_hydroparts:
        output.add(waterdeel)


singles = output.to_single(quiet=False)
fixed = singles.fix(quiet=False)
dissolved = fixed.dissolve(field="line_id", quiet=False)
dissolved.write("processing/processed_waterdeel4.shp")
