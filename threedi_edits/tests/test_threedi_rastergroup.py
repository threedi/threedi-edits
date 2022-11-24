# # -*- coding: utf-8 -*-
# """
# Created on Fri Feb 26 09:35:28 2021

# @author: chris.kerklaan

# #TODO functions:
#     1. Memory counter

# Currently supported functions
#     1. Check properties
#         1. nodata
#         2. projection
#         3. data_type
#         4. extreme_values
#         5. max_allowed_pixels --> takes a long time, so currently switched off
#         6. Not square pixels --> could not create this hence not present
#     2. Check null rasters
#     3. Check load csv
#     4. Checks generate functions
#         1. Friction
#         2. Interceptions
#         3. Infiltration (Permeability and max infiltration)
#         4. Crop type
#         5. Hydrauilic conductivty
#         6. Building interception

# """
# # First-party imports
# import pathlib

# # Third-party imports
# import numpy as np
# from osgeo import gdal
# import logging

# logger = logging.getLogger()


# # Local imports
# from threedi_edits.threedi.rastergroup import ThreediRasterGroup
# from threedi_edits.gis.raster import Raster
# from threedi_edits.gis.vector import Vector

# # Globals
# TEST_DIRECTORY = (
#     str(pathlib.Path(__file__).parent.absolute()) + "/data/threedi_raster_group/"
# )

# MANDATORY_TITLES_SOIL = [
#     "Title",
#     "Code",
#     "Soiltype",
#     "Max_infiltration_rate",
#     "Hydraulic_conductivity",
# ]
# MANDATORY_TITLES_LANDUSE = [
#     "Code",
#     "Landuse",
#     "Friction",
#     "Permeability",
#     "Interception",
#     "Crop_type",
# ]


# table_path = TEST_DIRECTORY + "conversietabellen/"
# csv_landuse_file = table_path + "Conversietabel_landgebruik_2020.csv"
# csv_soil_file = table_path + "Conversietabel_bodem.csv"
# csv_landuse_file_broken = table_path + "Conversietabel_landgebruik_2020_broken.csv"
# dem_file = TEST_DIRECTORY + "dem.tif"
# interception_file = TEST_DIRECTORY + "interception.tif"
# frict_coef_file = TEST_DIRECTORY + "friction.tif"
# infiltration_file = TEST_DIRECTORY + "infiltration.tif"
# soil_file = TEST_DIRECTORY + "soil.tif"
# landuse_file = TEST_DIRECTORY + "landuse.tif"
# crop_type_file = TEST_DIRECTORY + "crop_type.tif"
# hydraulic_conductivity_file = TEST_DIRECTORY + "hydraulic_conductivity.tif"
# buildings_path = TEST_DIRECTORY + "buildings.shp"
# building_interception_file = TEST_DIRECTORY + "building_interception.tif"


# def test_check_nodata_properties():
#     """tests if checks work for nodata value equals -9999"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )
#     group.dem.nodata_value = -8888
#     checks = group.check_properties()

#     assert len(checks["errors"]) > 0

#     if len(checks["errors"]) > 0:
#         assert checks["errors"][0][0] == "nodata_value"


# def test_check_projection_properties():
#     """tests if checks work for projection and units if it is 28992 and metere"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )

#     group.dem.spatial_reference = 4326
#     checks = group.check_properties()

#     errors_present = len(checks["errors"]) > 0
#     assert errors_present
#     if errors_present:
#         assert checks["errors"][0][0] == "unit"
#         assert checks["errors"][1][0] == "projection"


# def test_check_data_type_properties():
#     """tests if checks work for data_type"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )
#     group.insert(group.dem.copy(shell=True, data_type=gdal.GDT_Int16))

#     checks = group.check_properties()

#     errors_present = len(checks["errors"]) > 0
#     assert errors_present

#     if errors_present:
#         assert checks["errors"][-1][0] == "data_type"


# def test_check_extreme_values_properties():
#     """tests if checks work for extreme values"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )
#     group.dem.load_to_memory()
#     group.dem.array = group.dem.array + 1 * 10**8
#     checks = group.check_properties()

#     errors_present = len(checks["errors"]) > 0
#     assert errors_present

#     if errors_present:
#         assert checks["errors"][1][0] == "extreme_values"


# def test_check_max_allowed_pixels_properties():
#     """tests if checks work for maximum allowed pixels by making an empty copy"""
#     group = ThreediRasterGroup(dem_file=dem_file)
#     copy = Raster.from_scratch(
#         "new",
#         group.dem.geotransform,
#         100000,
#         10001,
#     )
#     group["dem"] = copy
#     checks = group.check_properties()
#     errors_present = len(checks["errors"]) > 0
#     assert errors_present
#     assert checks["errors"][-1][0] == "maximum_allowed_pixels"
#     group = None


# def test_null_rasters():
#     """tests null rasters"""
#     group = ThreediRasterGroup(dem_file=dem_file)
#     array = group.null_raster().array
#     assert np.nansum(array) == 0.0


# def test_load_csv_soil():
#     """tests load csv on keys, values and if exists in object"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         landuse_file=landuse_file,
#         soil_file=soil_file,
#     )

#     group.load_soil_conversion_table(csv_soil_file)
#     assert hasattr(group, "ct_soil")
#     if hasattr(group, "ct_soil"):
#         for title in MANDATORY_TITLES_SOIL:
#             assert title in group.ct_soil.keys()
#             values = group.ct_soil[title]
#             assert len(values) > 0


# def test_load_csv_lu():
#     """tests load csv on keys, values and if exists in object"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         landuse_file=landuse_file,
#         soil_file=soil_file,
#     )

#     group.load_landuse_conversion_table(csv_landuse_file)
#     assert hasattr(group, "ct_lu")

#     if hasattr(group, "ct_lu"):
#         for title in MANDATORY_TITLES_LANDUSE:
#             assert title in group.ct_lu.keys()

#             values = group.ct_lu[title]
#             assert len(values) > 0


# def test_generate_friction():
#     """tests if generated friction using total counts"""

#     group = ThreediRasterGroup(dem_file=dem_file, landuse_file=landuse_file)

#     group.load_landuse_conversion_table(csv_landuse_file)
#     group.generate_friction()
#     or_sum = np.nansum(Raster(frict_coef_file).array)
#     generated_sum = np.nansum(group.friction.array)
#     assert or_sum == generated_sum


# def test_generate_interception():
#     """tests if generated interceptions using total counts"""

#     group = ThreediRasterGroup(dem_file=dem_file, landuse_file=landuse_file)

#     group.load_landuse_conversion_table(csv_landuse_file)
#     group.generate_interception()

#     or_sum = np.nansum(Raster(interception_file).array)
#     generated_sum = np.nansum(group.interception.array)
#     assert or_sum == generated_sum


# def test_generate_infiltration():
#     """tests if generated infiltraion using total counts"""

#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         landuse_file=landuse_file,
#         soil_file=soil_file,
#     )

#     group.load_landuse_conversion_table(csv_landuse_file)
#     group.load_soil_conversion_table(csv_soil_file)
#     group.generate_infiltration()

#     or_sum = np.nansum(Raster(infiltration_file).array)
#     generated_sum = np.nansum(group.infiltration.array)
#     assert or_sum == generated_sum


# def test_generate_crop_type():
#     """tests if generated crop type using total counts"""

#     group = ThreediRasterGroup(dem_file=dem_file, landuse_file=landuse_file)

#     group.load_landuse_conversion_table(csv_landuse_file)
#     group.generate_crop_type()

#     or_sum = np.nansum(Raster(crop_type_file).array)
#     generated_sum = np.nansum(group.crop_type.array)
#     assert or_sum == generated_sum


# def test_generate_hydraulic_conductivity():
#     """tests if generated hydraulic conductivity using total counts"""

#     group = ThreediRasterGroup(dem_file=dem_file, soil_file=soil_file)

#     group.load_soil_conversion_table(csv_soil_file)
#     group.generate_hydraulic_conductivity()
#     or_sum = np.nansum(Raster(hydraulic_conductivity_file).array)
#     generated_sum = np.nansum(group.hydraulic_conductivity.array)
#     assert or_sum == generated_sum


# def test_generate_building_interception():
#     """tests if generated buildings interception using total counts"""

#     group = ThreediRasterGroup(dem_file, buildings_file=buildings_path)
#     group.generate_building_interception(10)
#     group.interception.write(TEST_DIRECTORY + "building_interception2.tif")
#     or_sum = np.nansum(Raster(building_interception_file).array)
#     generated_sum = np.nansum(group.interception.array)
#     assert or_sum == generated_sum


# def test_generate_building_interception_nodata():
#     """tests if generated buildings interception using 0"""

#     group = ThreediRasterGroup(dem_file, buildings_file=buildings_path)
#     group.generate_building_interception(-9999)

#     assert group.interception.read(Vector(buildings_path)[2].geometry)[0] == 0


# def test_correct():
#     """tests if property correcting is properly perfomed"""
#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )

#     group.interception.nodata_value = -8888
#     group.friction.spatial_reference = 28991
#     group.infiltration.change_data_type(int)

#     buildings = Vector(buildings_path)
#     group["infiltration"] = group.infiltration.push_vector(
#         buildings, set_to_nodata=True
#     )

#     checks = group.check_properties()
#     assert checks["errors"][1][0] == "nodata_value"
#     assert checks["errors"][2][0] == "projection"
#     assert checks["errors"][3][0] == "data_type"

#     group.correct()
#     checks = group.check_properties()

#     assert len(checks["errors"]) == 1 and checks["errors"][0][0] == "nodata_value"


# def test_write():
#     """tests write based on a dictionary"""
#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )
#     group.write(TEST_DIRECTORY)


# def test_write_w_dict():
#     """tests write based on a dictionary"""
#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )
#     group.write(TEST_DIRECTORY, output_names={"dem": "raster/dem.tif"})


# def test_clip():
#     """tests if the whole group is clipped"""
#     group = ThreediRasterGroup(
#         dem_file=dem_file,
#         interception_file=interception_file,
#         frict_coef_file=frict_coef_file,
#         infiltration_rate_file=infiltration_file,
#     )
#     group_sum = np.nansum(group.dem.array)
#     buildings = Vector(buildings_path)
#     group.clip(buildings)

#     assert group_sum < np.nansum(group.dem.array)


# def test_fix_table():
#     """tests if a broken csv is still loaded correctly
#     checks only if the lengths are still the same
#     """
#     group = ThreediRasterGroup()

#     group.load_landuse_conversion_table(csv_landuse_file_broken)
#     group.fix_table()

#     for c, (key, value) in enumerate(group.ct_lu.items()):

#         if c == 0:
#             base_len = len(value)
#         else:
#             if not key == "Landuse":
#                 assert len(value) == base_len


# # def test_memory_capacity():
# #     """ rastergroup should be able to handle three rasters of 500.000.000
# #         pixels each
# #     """
# #     group = ThreediRasterGroup(dem_file=dem_file,
# #                                 landuse_file=landuse_file,
# #                                 soil_file=soil_file)
# #     copy = group.dem.empty_copy(20000, 25000)
# #     group["dem"] = copy

# #     landuse = copy.empty_copy()
# #     for tile in copy:
# #         array = tile.array
# #         array[array==0]=1
# #         landuse.array = array, *tile.location
# #     landuse.name ="landuse"

# #     group["landuse"] = landuse

# #     soil = copy.empty_copy()
# #     for tile in copy:
# #         array = tile.array
# #         array[array==0]=1
# #         soil.array = array, *tile.location
# #     soil.name = "soil"
# #     group["soil"] = soil


# #     group.load_landuse_conversion_table(csv_landuse_file)
# #     group.load_soil_conversion_table(csv_soil_file)

# #     group.generate_friction()
# #     group.generate_permeability()
# #     group.generate_interception()
# #     group.generate_crop_type()
# #     group.generate_max_infiltration()
# #     group.generate_infiltration()
# #     group.generate_hydraulic_conductivity()

# #     assert True
# # if __name__ == "__main__":
# #     __file__ = "C:/Users/chris.kerklaan/Documents/Github/threedi-raster-edits/threedi_edits/tests/test_threedi_rastergroup.py"
# #     import sys

# #     sys.path.append("C:/Users/chris.kerklaan/Documents/Github/threedi-raster-edits")
