# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:40:16 2021

@author: chris.kerklaan


ThreediRasterGroup class checks for 3Di properties, alignment etc.

Input:
    1. Raster paths as string
    2. Optional panden

Checks:
    1. Check alignment --> ALready included in rastergroups
    2. Check 3Di properties  --> Done!

TODO:
    - None

Notes:
    1. Memory loading speeds up raster analysis, however is costly in memory
       Not all rasters are loaded into memory, e.g., it is not usefull to load
       a dem into memory since it is not used in conversion.
    2. Maximum memory is three rasters with 500.000.000 pixels each.

"""
# First-party imports
import csv
import shutil
import pathlib
import logging
from pathlib import Path

# Third-party imports
import numpy as np
from osgeo import gdal

# local imports
from threedi_edits.gis.rastergroup import RasterGroup
from threedi_edits.gis.raster import Raster
from threedi_edits.gis.vector import Vector
from threedi_edits.utils.project import Classes, Functions

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)

# GLOBALS
# Logger
logger = logging.getLogger(__name__)

# CSV path
FILE_PATH = str(pathlib.Path(__file__).parent.absolute()) + "/data/"
CSV_LANDUSE_PATH = FILE_PATH + "Conversietabel_landgebruik_2020.csv"
CSV_SOIL_PATH = FILE_PATH + "Conversietabel_bodem.csv"

# Output paths
OUTPUT_NAMES = {
    "friction": "rasters/friction.tif",
    "dem": "rasters/dem.tif",
    "infiltration": "rasters/infiltration.tif",
    "interception": "rasters/interception.tif",
    "intial_waterlevel": "rasters/initial_waterlevel.tif",
    "infiltration_capacity": "rasters/infiltration_capacity.tif",
    "permeability": "rasters/permeability.tif",
    "max_infiltration_rate": "rasters/max_infiltration_rate.tif",
    "crop_type": "rasters/crop_type.tif",
    "hydraulic_conductivity": "rasters/hydraulic_conductivity.tif",
}
THREEDI_RASTER_NAMES = {
    "dem_file": "dem",
    "interception_file": "interception",
    "frict_coef_file": "friction",
    "infiltration_rate_file": "infiltration",
    "initial_waterlevel_file": "initial_waterlevel",
    "max_infiltration_capacity_file": "infiltration_capacity",
}
THREEDI_RASTER_NAMES_REVERSED = {v: k for k, v in THREEDI_RASTER_NAMES.items()}


CT_LU_KEYS = [
    "Code",
    "Landuse",
    "Friction",
    "Permeability",
    "Interception",
    "Crop_type",
]
CT_LU_TYPES = [int, str, float, float, float, int]

CT_SOIL_KEYS = [
    "Title",
    "Code",
    "Soiltype",
    "Max_infiltration_rate",
    "Hydraulic_conductivity",
]
CT_SOIL_TYPES = [str, int, str, float, float]


class ThreediRasterGroup(RasterGroup):
    def __init__(
        self,
        dem_file: str = None,
        landuse_file: str = None,
        soil_file: str = None,
        interception_file: str = None,
        frict_coef_file: str = None,
        infiltration_rate_file: str = None,
        initial_waterlevel_file: str = None,
        wind_shielding_file: str = None,
        max_infiltration_capacity_file: str = None,
        buildings_file: str = None,
        nodata_value=-9999,
        data_type=gdal.GDT_Float32,
        np_data_type="f4",
        epsg=28992,
        memory=True,
    ):
        """
        ThreediRasterGroup -- info
        ThreediRasterGroup is an object to edit generate,
        check and correct threedi rasters
        General procedure is as follows:
            1. Load rasters by given up the path and create the class group = ThreediRasterGroup(landuse_file='path.tif')
            2. load a conversion csv group.load_landuse_conversion_table()
            3. generate a raster 'group.generate_friction()'
            4. check properties 'group.properties'
            5. check alignment (if you have multiple rasters) 'group.alignment'
            6. correct all mishaps 'group.correct()'
            7. write all rasters  'group.write(folder)'

        Please refer to 'help()' of the specific functions

        params:
            dem_file: str = None,
            landuse_file: str = None,
            soil_file: str = None,
            interception_file: str = None,
            frict_coef_file: str = None,
            infiltration_rate_file: str = None,
            initial_waterlevel_file: str = None,
            max_infiltration_capacity_file: str = None,
            buildings: Vector = None,
            nodata_value=-9999,
            data_type=gdal.GDT_Float32,
            np_data_type="f4",
            epsg=28992
            memory=True

        use self['name'] = Raster to add a raster if it is already an object
        use self.add(path) if you have only a path

        If memory if True, rasters will be loaded into memory.
        If false, they will not. This is slower but needed sometimes to
        preserve memory.
        """

        rasters = []
        self.original_names = {}

        if dem_file:
            logger.info("Loading dem")
            dem = Raster(dem_file, name="dem")
            if memory:
                dem.load_to_memory()
            rasters.append(dem)
            self.original_names["dem"] = Path(dem_file).stem

        if landuse_file:
            logger.info("Loading landuse")
            landuse = Raster(landuse_file, name="landuse")
            rasters.append(landuse)
            self.original_names["landuse"] = Path(landuse_file).stem

        if soil_file:
            logger.info("Loading soil")
            soil = Raster(soil_file, name="soil")
            rasters.append(soil)
            self.original_names["soil"] = Path(soil_file).stem

        if interception_file:
            logger.info("Loading interception")
            interception = Raster(interception_file, name="interception")
            if memory:
                interception.load_to_memory()
            rasters.append(interception)
            self.original_names["interception"] = Path(interception_file).stem

        if frict_coef_file:
            logger.info("Loading friction")
            friction = Raster(frict_coef_file, name="friction")
            if memory:
                friction.load_to_memory()
            rasters.append(friction)
            self.original_names["friction"] = Path(frict_coef_file).stem

        if infiltration_rate_file:
            logger.info("Loading infiltration")
            infiltration = Raster(infiltration_rate_file, name="infiltration")
            if memory:
                infiltration.load_to_memory()
            rasters.append(infiltration)
            self.original_names["infiltration"] = Path(infiltration_rate_file).stem

        if max_infiltration_capacity_file:
            logger.info("Loading infiltration capacity")
            infiltration_capacity = Raster(
                max_infiltration_capacity_file, name="infiltration_capacity"
            )
            if memory:
                infiltration_capacity.load_to_memory()
            rasters.append(infiltration_capacity)
            self.original_names["infiltration_capacity"] = Path(
                max_infiltration_capacity_file
            ).stem

        if initial_waterlevel_file:
            logger.info("Loading initial waterlevel")
            ini_wl = Raster(initial_waterlevel_file, name="initial_waterlevel")
            if memory:
                ini_wl.load_to_memory()
            rasters.append(ini_wl)
            self.original_names["initial_waterlevel"] = Path(
                initial_waterlevel_file
            ).stem

        if wind_shielding_file:
            logger.info("Loading wind shielding")
            wind_shield = Raster(wind_shielding_file, name="wind_shielding")
            if memory:
                wind_shield.load_to_memory()
            rasters.append(wind_shield)
            self.original_names["wind_shielding"] = Path(wind_shielding_file).stem

        RasterGroup.__init__(self, rasters)

        if buildings_file:
            logger.debug("Setting buildings")
            self.buildings = Vector(buildings_file)

        if hasattr(self, "dem"):
            self.epsg = self.dem.epsg
        else:
            self.epsg = epsg
        self.data_type = data_type
        self.no_data_type = np_data_type
        self.nodata_value = nodata_value

        self.retrieve_soil_conversion_table = retrieve_soil_conversion_table
        self.retrieve_landuse_conversion_table = retrieve_landuse_conversion_table

        self.memory = memory

    @property
    def properties(self):
        if not hasattr(self, "_properties"):
            self.check_properties()

        return self._properties

    @property
    def alignment(self):
        if not hasattr(self, "_alignment"):
            self.check_alignment()

        return self._alignment

    @property
    def base_raster(self):
        if hasattr(self, "dem"):
            logger.info("Using dem as base/alignment raster")
            return self.dem
        else:
            logger.info(f"Using the {self.rasters[0].name} as base/alignment raster")
            return self.rasters[0]

    def __str__(self):
        return f"3Di rasters {self.names}"

    def __repr__(self):
        return f"3Di rasters {self.names} \n{self.properties}"

    def add(self, name: str, path: str):
        """adds a raster to the object
        params:
            name: name of raster
            path: path of raster
        """
        if not (name in THREEDI_RASTER_NAMES.values() or name in THREEDI_RASTER_NAMES):
            raise ValueError("name should be object or threedi compatibale")

        if name in THREEDI_RASTER_NAMES:
            name = THREEDI_RASTER_NAMES[name]

        if type(path) == str:
            self[name] = Raster(path)
        else:
            self[name] = path

    def check_table(self, table="soil"):
        """Checks if attribute tables are present, raises an attribute error
        if not present.
        params:
            table: 'soil' or 'landuse'
        """
        logger.info("Checking tables")
        if table == "soil":
            if not hasattr(self, "ct_soil"):
                raise AttributeError(
                    """
                                     Please load soil csv using
                                     'load_soil_conversion_table'"""
                )

        elif table == "landuse":
            if not hasattr(self, "ct_lu"):
                raise AttributeError(
                    """
                                     Please load landuse csv using
                                     'load_landuse_conversion_table'"""
                )

    def fix_table(self):
        """
        Fixes table based on the presence of keys in the current dict
        and fixes the types of the values based on what should be present

        """
        logger.info("Fixing table")

        if hasattr(self, "ct_lu"):
            new_ct_lu = {}
            for key, value in self.ct_lu.items():

                # fix keys
                for fixed_key in CT_LU_KEYS:
                    if fixed_key in key:
                        key = fixed_key

                # fix value types
                key_type = CT_LU_TYPES[CT_LU_KEYS.index(key)]
                new_values = []
                for v in value:
                    if key_type in [float, int]:
                        if isnumber(v):
                            new_values.append(key_type(v))
                    else:
                        new_values.append(key_type(v))

                new_ct_lu[key] = new_values

            self.ct_lu = new_ct_lu

        if hasattr(self, "ct_soil"):
            new_soil_lu = {}
            for key, value in self.ct_soil.items():

                # fix keys
                for fixed_key in CT_LU_KEYS:
                    if fixed_key in key:
                        key = fixed_key

                # fix value types
                key_type = CT_SOIL_TYPES[CT_SOIL_KEYS.index(key)]
                new_values = []
                for v in value:
                    if key_type in [float, int]:
                        if isnumber(v):
                            new_values.append(key_type(v))
                    else:
                        new_values.append(key_type(v))

                new_soil_lu[key] = new_values

            self.ct_soil = new_soil_lu

    def check_properties(
        self, min_allow=-1000, max_allow=1000, max_pixel_allow=1000000000
    ):
        """
        Checks if properties of the rasters are adequate for 3Di
        E.g, epsg is in metre, nodata must be -9999

        params:
            min_allow: minimum allow value of the raster
            max_allow: maximum value of the raster
            max_pixel_allow: total allowed amount of pixels in all rasters
        """
        logger.info("Checking properties")
        self._properties = check_properties(
            self.rasters,
            nodata=self.nodata,
            projection=self.epsg,
            data_type=self.data_type,
            min_allow=min_allow,
            max_allow=max_allow,
            max_pixel_allow=max_pixel_allow,
        )
        return self._properties

    def null_raster(self):
        """Creates a rasters with zeros by using the dem"""
        logger.info("Creating null raster")
        copy = self.base_raster.copy(shell=True, memory=self.memory)
        null_array = np.zeros((int(copy.rows), int(copy.columns)))
        null_array[~self.dem.mask] = np.nan
        copy.array = null_array
        copy.name = "null"
        return copy

    def load_soil_conversion_table(self, csv_soil_path=CSV_SOIL_PATH):
        """loads the soil conversion table, if nothing is entered a standard
        soil table is loaded
        params:
            csv_soil_path
        """
        logger.info("Loading soil conversion table")
        self.ct_soil, self.ct_soil_info = load_csv_file(csv_soil_path, "soil")

    def load_landuse_conversion_table(self, csv_lu_path=CSV_LANDUSE_PATH):
        """loads the soil conversion table, if nothing is entered a standard
        soil table is loaded
        params:
            csv_soil_path
        """
        logger.info("Loading landuse conversion table")
        self.ct_lu, self.ct_lu_info = load_csv_file(csv_lu_path, "landuse")

    def generate_friction(self):
        """Generates a friction raster based on the landuse table"""
        logger.info("Generating friction")
        if hasattr(self, "friction"):
            return

        self.friction = classify(self.landuse, "Friction", self.ct_lu, self.memory)
        self["Friction"] = self.friction

    def generate_permeability(self):
        """Generates a permeability raster based on the landuse table"""
        logger.info("Generating permeability")
        if hasattr(self, "permeability"):
            return
        self.permeability = classify(
            self.landuse, "Permeability", self.ct_lu, self.memory
        )
        self["Permeability"] = self.permeability

    def generate_interception(self):
        """Generates an interception raster based on the landuse table"""
        logger.info("Generating interception")
        if hasattr(self, "interception"):
            return
        self.interception = classify(
            self.landuse, "Interception", self.ct_lu, self.memory
        )
        self["Interception"] = self.interception

    def generate_crop_type(self):
        """Generates a crop type raster based on the landuse table"""
        logger.info("Generating crop type")
        if hasattr(self, "crop_type"):
            return

        self.crop_type = classify(self.landuse, "Crop_type", self.ct_lu, self.memory)
        self["Crop_type"] = self.crop_type

    def generate_max_infiltration(self):
        """Generates a max infiltration raster based on the soil table"""
        logger.info("Generating max infiltration")
        if hasattr(self, "max_infiltration_rate"):
            return

        print("generating max infiltration")
        self.max_infiltration_rate = classify(
            self.soil, "Max_infiltration_rate", self.ct_soil, self.memory
        )
        self["Max_infiltration_rate"] = self.max_infiltration_rate

    def generate_infiltration(self):
        """Generates an infiltration raster based permeability and max infiltration"""
        logger.info("Generating infiltration")
        if hasattr(self, "infiltration"):
            return None

        self.generate_permeability()
        self.generate_max_infiltration()

        output = self.base_raster.copy(shell=True, memory=self.memory)
        for perm_tile, inf_tile in zip(self.permeability, self.max_infiltration_rate):
            perm_array = perm_tile.array
            inf_array = inf_tile.array

            infiltration_array = np.where(
                np.logical_and(
                    perm_array != self.nodata_value,
                    inf_array != self.nodata_value,
                ),
                perm_array * inf_array,
                self.nodata_value,
            ).astype(self.np_data_type)
            output.array = infiltration_array, *perm_tile.location

        self.infiltration = output
        self.infiltration.name = "Infiltration"
        self["Infiltration"] = self.infiltration

    def generate_hydraulic_conductivity(self):
        """Generates a hydraulic conductivity raster based on the soil table"""
        logger.info("Generating hydraulic conductivity")
        if hasattr(self, "hydraulic_conductivity"):
            return
        self.hydraulic_conductivity = classify(
            self.soil, "Hydraulic_conductivity", self.ct_soil, self.memory
        )
        self["Hydraulic_conductivity"] = self.hydraulic_conductivity

    def generate_building_interception(self, value):
        """Generates buildings interception
        params:
            value: value which is placed on top of the buildings
        """
        if hasattr(self, "interception"):
            return

        logger.info("Fixing building vector")
        buildings = self.buildings.fix()

        logger.info("Generating building interception")
        null = self.null_raster()
        value_raster = null.push_vector(buildings, value=value)
        self.interception = value_raster.align(
            self.dem, nodata_align=True, fill_value=0
        )
        self.interception.name = "Interception"
        self["Interception"] = self.interception

    def generate_all(self):
        """generates all rasters except building interception"""

        logger.info("Generating all possible rasters")
        self.generate_friction()
        self.generate_permeability()
        self.generate_interception()
        self.generate_crop_type()
        self.generate_max_infiltration()
        self.generate_infiltration()
        self.generate_hydraulic_conductivity()

    def correct(self, nodata=True, projection=True, data_type=True, alignment=True):
        """Corrects the rasters based on nodata value, epsg, data type and
        Data/Nodata pixels compared to the dem
        params:
            nodata: boolean, if false is skipped
            projection: boolean, if false is skipped
            data_type: boolean, if false is skipped
        """

        if self.count == 0:
            return

        logger.info("Correcting properties of ThreediGroup rasters")
        checks = self.check_properties()

        nodata_value = self.nodata_value
        spatial_reference = self.epsg
        data_type = self.data_type

        for error in checks["errors"]:
            error_type, raster_name, _ = error
            if error_type == "projection" and projection:
                logger.debug(
                    f"replacing projection of {raster_name} to {spatial_reference}"
                )
                raster = getattr(self, raster_name.lower())
                raster.spatial_reference = spatial_reference
                self[raster_name] = raster

            elif error_type == "data_type" and data_type:
                logger.debug(f"replacing data type of {raster_name} to {data_type}")
                raster = getattr(self, raster_name.lower())
                raster.change_data_type(data_type)
                self[raster_name] = raster

        if alignment:
            alignment_check = self.check_alignment()
            for raster in self:
                if raster.name in alignment_check["to_be_aligned"]:
                    logger.debug(f"Aligning {raster.name} to dem")
                    self[raster.name] = raster.align(
                        self.base_raster, idw=False, fill_value=True, nodata_align=True
                    )

        for error in checks["errors"]:
            error_type, raster_name, _ = error
            if error_type == "nodata_value" and nodata:
                logger.debug(f"replacing nodata of {raster_name} to {nodata_value}")
                raster = getattr(self, raster_name.lower())
                raster.replace_nodata(nodata_value, quiet=True)
                self[raster_name] = raster

    def write(self, folder_path, output_names: dict = OUTPUT_NAMES):
        """
           writes raster if present in output names and in object
           creates a folder if it has to...

        params:
           folder_path: path in which the rasters are written
           output_names: e.g. {
           "friction": "rasters/friction.tif",
           "dem":"rasters/dem.tif",
           "infiltration": "rasters/infiltation.tif",
           "interception": "rasters/interception.tif",
           "intial_waterlevel": "rasters/initial_waterlevel.tif"
           }
        """

        # make folder path
        p = pathlib.Path(folder_path)
        p.mkdir(parents=True, exist_ok=True)

        for k, v in output_names.items():
            if v is not None:
                (pathlib.Path(folder_path) / v.split("/")[0]).mkdir(
                    parents=True, exist_ok=True
                )

        for table_name, table_path in output_names.items():
            if table_name in THREEDI_RASTER_NAMES:
                table_name = THREEDI_RASTER_NAMES[table_name]

            if hasattr(self, table_name):
                logger.info(f"Writing {table_name}")
                table = getattr(self, table_name)
                table.write(str(p / table_path))


def isnumber(number):
    """It is quite hard to check if it is a number
    This is currently the best approach
    """
    try:
        float(number)
    except ValueError:
        return False
    else:
        return True


def retrieve_soil_conversion_table(output_path):
    logger.info("retrieving soil conversion table")
    shutil.copyfile(CSV_SOIL_PATH, output_path)


def retrieve_landuse_conversion_table(output_path):
    logger.info("retrieving landuse conversion table")
    shutil.copyfile(CSV_LANDUSE_PATH, output_path)


def check_properties(
    raster_list,
    nodata=-9999,
    projection=28992,
    max_pixel_allow=1000000000,
    data_type=gdal.GDT_Float32,
    unit="metre",
    min_allow=-1000,
    max_allow=1000,
):
    output = {
        "nodata": {},
        "unit": {},
        "projection": {},
        "data_type": {},
        "resolution": {},
        "square_pixels": {},
        "min_max": {},
        "total_pixels": 0,
        "errors": [],
    }

    total_pixels = 0
    for raster in raster_list:

        # nodata value check
        if raster.nodata_value != nodata:
            msg = (
                "nodata_value",
                raster.name,
                f"has a nodata value of {raster.nodata_value}",
            )
            logger.debug(msg)
            output["errors"].append(msg)
        output["nodata"][raster.name] = raster.nodata_value

        raw_nan_count = np.nansum(np.isnan(raster.band.ReadAsArray()))
        if raw_nan_count > 0:
            msg = (
                "nodata_value",
                raster.name,
                " has raw np.nan values",
            )
            logger.debug(msg)
            output["errors"].append(msg)
        output["nodata"][raster.name] = raster.nodata_value

        # unit check
        if raster.spatial_reference.unit != unit:
            msg = ("unit", raster.name, f"has not unit {unit}")
            logger.debug(msg)
            output["errors"].append(msg)
        output["unit"][raster.name] = raster.spatial_reference.unit

        # projection check
        if raster.spatial_reference.epsg != projection:
            msg = ("projection", raster.name, f"has not epsg {projection}")
            logger.debug(msg)
            output["errors"].append(msg)
        output["projection"][raster.name] = raster.spatial_reference.epsg

        # data type check
        if raster.data_type != data_type:
            msg = (
                "data_type",
                raster.name,
                f"is not a {gdal.GetDataTypeName(data_type)}",
            )
            logger.debug(msg)
            output["errors"].append(msg)
        output["data_type"][raster.name] = gdal.GetDataTypeName(data_type)

        # square pixel check
        if abs(raster.resolution["width"]) != abs(raster.resolution["height"]):
            msg = ("width/height", raster.name, "has not a square pixel")
            output["errors"].append(msg)
        output["square_pixels"][raster.name] = raster.resolution

        # extreme value check
        _max, _min = np.nanmax(raster.array), np.nanmin(raster.array)
        if not (min_allow < _max < max_allow and min_allow < _min < max_allow):
            msg = (
                "extreme_values",
                raster.name,
                f"has extreme values < {min_allow}, > {max_allow}",
            )
            logger.debug(msg)
            output["errors"].append(msg)
        output["min_max"][raster.name] = {"min": _min, "max": _max}

        total_pixels += raster.pixels

    # max pixel allowed check
    if total_pixels > max_pixel_allow:
        msg = f"Rasters combined pixels are larger than {max_pixel_allow}"
        logger.debug(msg)
        output["errors"].append(("maximum_allowed_pixels", msg, ""))

    output["total_pixels"] = total_pixels

    if len(output["errors"]) == 0:
        logger.debug("ThreediRasterGroup - Check properties found no problems")
    return output


def load_csv_file(csv_path, csv_type="landuse"):
    csv_data = {}
    csv_info = {}
    if csv_type == "landuse":
        csv_structure = {
            1: "description",
            2: "unit",
            3: "range",
            4: "type",
        }
        meta_list = [1, 2, 3, 4]
    elif csv_type == "soil":
        csv_structure = {
            1: "description",
            2: "source",
            3: "unit",
            4: "range",
            5: "type",
        }
        meta_list = [1, 2, 3, 4, 5]

    with open(csv_path) as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=";")
        for i, line in enumerate(csv_reader):

            # headers
            if i == 0:
                headers = line
                for column_value in line:
                    csv_data[column_value] = []
                    csv_info[column_value] = {}

            # units, descriptions, ranges etc.
            elif i in meta_list:

                for column_index, column_value in enumerate(line):
                    field = csv_structure[i]
                    csv_info[headers[column_index]][field] = column_value

            # csv data
            else:
                for column_index, column_value in enumerate(line):
                    column = headers[column_index]
                    column_type = csv_info[column]["type"]

                    if column_value == "":
                        csv_data[column].append(None)
                    else:
                        if column_type == "Integer":
                            column_value = int(column_value)
                        elif column_type == "Double" or column_type == "Real":
                            column_value = float(column_value)
                        elif column_type == "String":
                            column_value = str(column_value)
                        csv_data[column].append(column_value)

        return csv_data, csv_info


def classify(raster: Raster, table: str, ct: dict, memory=False):
    """input is a table to classify,
    raster is the template raster,
    ct is the conversion table
    returns a dictionary of classified rasters"""

    ct_codes = ct["Code"]
    ct_table = ct[table]
    output = raster.copy(shell=True, memory=memory)
    for tile in raster:
        array = tile.array

        output_array = classify_array(array, ct_table, ct_codes)

        if output_array is None:
            output.array = array, *tile.location
        else:
            output.array = output_array, *tile.location

        tile.close()

    output.name = table
    return output


def classify_array(array, ct_table, ct_codes):
    if np.isnan(array).all():
        return None

    codes = np.unique(array[~np.isnan(array)])
    output_array = np.copy(array)
    for code in codes:
        try:
            output_array[array == code] = ct_table[ct_codes.index(code)]
        except IndexError:
            raise IndexError(f"""{code} not found in table""")
        # array = None # dit levert problemen op

    return output_array
