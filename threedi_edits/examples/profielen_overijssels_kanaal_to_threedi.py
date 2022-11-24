# -*- coding: utf-8 -*-
"""
Created on Wed May 25 09:03:46 2022

@author: stijn.overmeen

Converts channel profiles to threedi cross section definitions and locations
Add these to the threedi model
Create vertex in channel

Bonus:
    Search for unused cross section defintions and delete these
"""

import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

plt.ioff()  # don't show plots in console
import numpy as np
import sys

sys.path.append(r"C:/Users/stijn.overmeen/Documents/GitHub/threedi-raster-edits")
import threedi_edits as tre
from osgeo import ogr, osr
from threedi_edits.threedi.tables.templates import Templates

template = Templates("integrated")

# packages for logging
import datetime, time, logging

# Define working dir
FOLDER = Path(
    "C:/Users/stijn.overmeen/Documents/Projecten_lokaal/Extern/overijssels_kanaal/scripts"
)

# Logging
LOG_FOLDER = os.path.join(FOLDER, "log")
if not os.path.exists(LOG_FOLDER):
    os.makedirs(LOG_FOLDER)
lines = ["Profielen Overijssels kanaal to threedi", ""]
log = (
    f"{datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d__%H-%M-%S')}"
    ".log"
)
with open(os.path.join(FOLDER, "log", log), "w") as f:
    for line in lines:
        f.write(line)
        f.write("\n")

logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.FileHandler(os.path.join(FOLDER, "log", log)),
        logging.StreamHandler(),
    ],
    format="%(message)s",
    force=True,
)


def strings_to_floats(data):
    """
    converts strings to floats in the dataframe (specific actions per column name)
    """
    for col in columns_to_convert:
        if col in ["Afstand", "z"]:
            data[col] = pd.to_numeric(data[col].str.replace(",", "."))
        elif col in ["X", "Y"]:
            temp = data[col].str.replace(".", "")
            for i, v in temp.items():
                temp[i] = float(v[:6] + "." + v[6:])
            data[col] = pd.to_numeric(temp)
    return data


def yz_array_to_threedi_tabulated_definitions(yz_array, step_size):
    """
    convert yz array to threedi tabulated trapezium definitions
    threedi requires symmetrical profiles and continuously vertically increasing volume
    calculates flow profile for every vertical step and adapts definition accordingly
    """
    bottom_level = min(yz_array[:, 1])
    ind_middle = np.argmin(yz_array[:, 1])
    left_bank_level = max(yz_array[: ind_middle + 1, 1])
    right_bank_level = max(yz_array[ind_middle:, 1])
    left_bank_above_bottom = left_bank_level - bottom_level
    right_bank_above_bottom = right_bank_level - bottom_level
    left_bank = min(yz_array[:, 0])
    yz_array[:, 1] = yz_array[:, 1] - bottom_level
    yz_array[:, 0] = yz_array[:, 0] - left_bank
    min_z = min(yz_array[:, 1])
    max_z = min(left_bank_above_bottom, right_bank_above_bottom)  # max(yz_array[:,1])
    steps = np.arange(min_z, max_z, step_size) + step_size
    channel_area_under_step = np.empty(steps.shape)
    width = np.empty(steps.shape)
    for i, step in enumerate(steps):
        y = yz_array[:, 0]
        z = yz_array[:, 1]
        z_clip = np.clip(z, 0, step)
        area_under_clip = np.trapz(z_clip, y, dx=0.01)
        area_under_step = step * (max(y) - min(y))
        channel_area_under_step[i] = area_under_step - area_under_clip
    area_increase = channel_area_under_step[1:] - channel_area_under_step[:-1]
    width = np.insert(area_increase / step_size, 0, 0, axis=0)
    height = np.linspace(0, max_z, width.shape[0])
    middle = max(width) / 2
    yz_tab_list = []
    for h, w in zip(height, width):
        if h == 0:
            yz_tab_list.append([middle, h])
        else:
            z = h
            yl = middle - w / 2
            yr = middle + w / 2
            yz_tab_list.append([yl, z])
            yz_tab_list.append([yr, z])
    yz_tab = np.array(yz_tab_list)
    yz_tab_left = yz_tab[yz_tab[:, 0] < middle]
    yz_tab_right = yz_tab[yz_tab[:, 0] >= middle]
    yz_tab_left_sorted = yz_tab_left[
        np.lexsort((yz_tab_left[:, 0], -yz_tab_left[:, 1]))
    ]
    yz_tab_right_sorted = yz_tab_right[
        np.lexsort((yz_tab_right[:, 0], yz_tab_right[:, 1]))
    ]
    yz_tab_sorted = np.concatenate((yz_tab_left_sorted, yz_tab_right_sorted), axis=0)
    np.set_printoptions(precision=2)
    widthstring = np.array2string(
        width, formatter={"float_kind": lambda width: "%.2f" % width}
    )
    heightstring = np.array2string(
        height, formatter={"float_kind": lambda height: "%.2f" % height}
    )

    return (
        left_bank_level,
        right_bank_level,
        widthstring,
        heightstring,
        bottom_level,
        yz_tab_sorted,
    )


def plot_profile_comparisons(profile_name, yz_tab_sorted, yz_array):
    """
    plots and save figures of profile comparisons, original yz against threedi tabulated definition
    """
    plt.figure(figsize=(10, 5))
    plt.title(
        "Comparison between YZ and tabulated trapezium definitions for "
        f"{profile_name}"
    )
    plt.plot(yz_tab_sorted[:, 0], yz_tab_sorted[:, 1], label="Tabulated trapezium")
    plt.plot(yz_array[:, 0], yz_array[:, 1], label="YZ")
    plt.xlabel("Distance from bank [m]", fontdict=None, labelpad=None)
    plt.ylabel("Depth [m]", fontdict=None, labelpad=None)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.savefig(
        os.path.join(FOLDER, "profiel_figuren", "slib", f"{profile_name}"),
        bbox_inches="tight",
    )


def store_information_cross_sections_in_dict(
    profile_name,
    profile_df,
    left_bank_level,
    right_bank_level,
    widthstring,
    heightstring,
    bottom_level,
):
    """
    returns dict with relevant information for further processing of information to 3Di model
    """
    cross_section = {}
    cross_section["code"] = profile_name
    cross_section["height"] = heightstring
    cross_section["width"] = widthstring
    cross_section["bank_level"] = min(left_bank_level, right_bank_level)
    cross_section["code"] = profile_name
    cross_section["reference_level"] = bottom_level
    cross_section["X"] = profile_df["X"].iloc[0]
    cross_section["Y"] = profile_df["Y"].iloc[0]

    return cross_section


def add_cross_section_definition_to_model(model, cs_def_profile, shape):
    """
    add cross section definition

    needed:
        code, width string, height string from profile defintion dict
        shape (integer)
    """
    cs_def_threedi = template.cross_section
    cs_def_threedi["code"] = cs_def_profile["code"]
    cs_def_threedi["shape"] = shape  # tabulated trapezium
    cs_def_threedi["width"] = cs_def_profile["width"].replace("[", "").replace("]", "")
    cs_def_threedi["height"] = (
        cs_def_profile["height"].replace("[", "").replace("]", "")
    )
    model.cross_section_definitions.add(cs_def_threedi)

    return model


def add_cross_section_location_to_model(
    model, cs_loc_profile, friction_type, friction_value, EPSG_number
):
    """
    add cross section location, geometry based on X, Y coordinates
    finds closest channel based on spatial search, and adds this channel id to the table

    needed:
        code, reference level, bank level, X and Y coordinates from profile defintion dict
        EPSG code beloging to X and Y coordinates
        friction type and friction value
    """
    cs_loc_threedi = template.cross_section_location
    cs_loc_threedi["code"] = cs_loc_profile["code"]
    # definition id based on profile code
    for cs_def_threedi in subselection_of_cs_definitions:
        if cs_def_threedi["code"] == cs_loc_profile["code"]:
            cs_loc_threedi["definition_id"] = cs_def_threedi.fid
            break
    cs_loc_threedi["reference_level"] = cs_loc_profile["reference_level"]
    cs_loc_threedi["bank_level"] = cs_loc_profile["bank_level"]
    cs_loc_threedi["friction_type"] = friction_type
    cs_loc_threedi["friction_value"] = friction_value

    # TODO
    # Tip!
    # zo zou ik het doen:
    point = tre.Point.from_point(
        (float(cs_loc_profile["X"]), float(cs_loc_profile["Y"])), epsg=28992
    )
    point = point.reproject(4326)

    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(float(cs_loc_profile["X"]), float(cs_loc_profile["Y"]))

    # transform point
    source = osr.SpatialReference()
    source.ImportFromEPSG(EPSG_number)
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)

    # The line below sets the TRADITIONAL_GIS_ORDER.
    target.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    transform = osr.CoordinateTransformation(source, target)
    point.Transform(transform)

    cs_loc_threedi["the_geom"] = tre.gis.point.Point(point)

    # TODO
    # Tip!, je kunt spatial filteren op channels. Ik kan me voorstellen dat de onderste code wat traag is.
    # Zo zou ik het doen met een prefilter:
    convert = tre.threedi.utils.conversions.METER_TO_DEGREE
    buffer = convert * 100  # meter
    output = []
    while len(output) == 0:
        output = model.channels.spatial_filter(point.Buffer(buffer), return_vector=True)
        buffer = buffer + convert * 100  # extra 100 meter zoekafstand

    # find closest channel
    smallest_distance = 99999999  # initialise big
    for channel in output:
        distance = point.Distance(channel.geometry)
        if distance < smallest_distance:
            smallest_distance = distance
            closest_channel = channel

    # find closest channel
    smallest_distance = 99999999  # initialise big
    for channel in model.channels:
        distance = point.Distance(channel.geometry)
        if distance < smallest_distance:
            smallest_distance = distance
            closest_channel = channel

    # add channel_id to cross section location
    cs_loc_threedi["channel_id"] = closest_channel.fid
    model.cross_section_locations.add(cs_loc_threedi)

    return model, closest_channel, point


if __name__ == "__main__":

    file_hoofdtak = Path(
        os.path.join(FOLDER, "data", "profielen_ovkanaal (bewerking2)_afvoer_slib.csv")
    )  # hoofdkanaal profielen

    file_zijtak_profiles = Path(
        os.path.join(FOLDER, "data", "inmeting 2020 - lettele.xlsx")
    )  # zijtak profielen
    file_zijtak_coord = Path(
        os.path.join(FOLDER, "data", "PROFILES_Locations.csv")
    )  # bijbehorende coordinaten

    hoofdtak_df = pd.read_csv(file_hoofdtak, sep=";", decimal=",")  # read file

    columns_to_convert = ["Afstand", "z", "X", "Y"]
    hoofdtak_df_improved = strings_to_floats(hoofdtak_df)

    zijtak_profiles_df = pd.read_excel(
        file_zijtak_profiles, sheet_name="bagger", decimal=","
    )  # read file
    zijtak_coord_df = pd.read_csv(file_zijtak_coord, sep=";", decimal=".")

    rename_dict_coord = {"ID": "profielnaam", "x": "X", "y": "Y"}

    rename_dict_profiles = {
        "Profiel": "profielnaam",
        "x ": "Afstand",
        "Type": "type",
        "Bron": "bron",
    }

    zijtak_profiles_df.rename(columns=rename_dict_profiles, inplace=True)
    zijtak_coord_df.rename(columns=rename_dict_coord, inplace=True)

    zijtak_df = zijtak_profiles_df.merge(zijtak_coord_df, on="profielnaam")

    data = pd.concat([hoofdtak_df_improved, zijtak_df], ignore_index=True)

    yz_profiles = []  # create empty list for YZ-profiles
    cross_sections = []  # create empty list for storing cross section information
    profiles = data.groupby(data["profielnaam"])  # group data by profile name

    logging.info("converting profiles to threedi tabulated definitions...")

    # TODO
    # Tip!
    # ik zou loopen met tre.Progress(profiles, "converting profiles to threedi tabulated defintions")
    for profile in profiles:

        profile_name = profile[0]
        logging.info(profile_name)
        profile_df = profile[1]
        yz_array = profile_df[["Afstand", "z"]].to_numpy()

        (
            left_bank_level,
            right_bank_level,
            widthstring,
            heightstring,
            bottom_level,
            yz_tab_sorted,
        ) = yz_array_to_threedi_tabulated_definitions(yz_array=yz_array, step_size=0.1)
        # plot_profile_comparisons(profile_name,yz_tab_sorted,yz_array)

        cross_section = store_information_cross_sections_in_dict(
            profile_name=profile_name,
            profile_df=profile_df,
            left_bank_level=left_bank_level,
            right_bank_level=right_bank_level,
            widthstring=widthstring,
            heightstring=heightstring,
            bottom_level=bottom_level,
        )
        cross_sections.append(cross_section)

    sqlite_path = r"C:/Users/stijn.overmeen/Documents/Projecten_lokaal/Extern/overijssels_kanaal/speeltuin/zandkasteel.sqlite"
    model = tre.ThreediEdits(sqlite_path, mode="memory")

    logging.info("\nadding cross section definitions...")
    codes = []  # appending relevant codes for later use in for-loop
    for cs_def_profile in cross_sections:
        codes.append(cs_def_profile["code"])
        model = add_cross_section_definition_to_model(
            model=model, cs_def_profile=cs_def_profile, shape=6
        )
        logging.info(f"added cross section definition for {cs_def_profile['code']}")

    # get subselection of cross section definitions you just added for looping over while adding cross section locations
    subselection_of_cs_definitions = []

    for cs_def_threedi in model.cross_section_definitions:
        if cs_def_threedi["code"] in codes:
            subselection_of_cs_definitions.append(cs_def_threedi)

    logging.info("\nadding cross section locations...")
    for cs_loc_profile in cross_sections:
        logging.info("\n" + cs_loc_profile["code"])
        model, closest_channel, point = add_cross_section_location_to_model(
            model=model,
            cs_loc_profile=cs_loc_profile,
            friction_type=2,
            friction_value=0.033,
            EPSG_number=28992,
        )
        logging.info("added cross section location")

        # delete existing cross section locations of this channel, but not the ones we might have added previously
        # so search for code name starting with any of the following prefixes, keep these
        prefixes = ["prof_D2018", "AS1", "AS2"]
        for cs_loc_threedi_for_deleting_old_ones in model.cross_section_locations:
            if (
                cs_loc_threedi_for_deleting_old_ones["channel_id"]
                == closest_channel.fid
            ):
                if not cs_loc_threedi_for_deleting_old_ones["code"].startswith(
                    tuple(prefixes)
                ):
                    model.cross_section_locations.delete(
                        cs_loc_threedi_for_deleting_old_ones
                    )
                    logging.info(
                        "deleted old cross section location "
                        + cs_loc_threedi_for_deleting_old_ones["code"]
                    )

        logging.info("done, no old cross section locations left to delete")

        # add vertex to channel
        closest_channel.geometry = closest_channel.geometry.add_vertice(
            tre.gis.point.Point(point)
        )
        logging.info(
            f"added vertex to channel {closest_channel.fid} of cross section location"
        )

    # delete unused cross section definitions
    # sample used
    used_cs_def_ids = []
    for cs_loc in model.cross_section_locations:
        used_cs_def_ids.append(cs_loc["definition_id"])
    for orifice in model.orifices:
        used_cs_def_ids.append(orifice["cross_section_definition_id"])
    for culvert in model.culverts:
        used_cs_def_ids.append(culvert["cross_section_definition_id"])
    for weir in model.weirs:
        used_cs_def_ids.append(weir["cross_section_definition_id"])

    for cs_df in model.cross_section_definitions:
        if not cs_df.fid in used_cs_def_ids:  # delete if not used
            model.cross_section_definitions.delete(cs_df)

    # write features
    # model.cross_section_locations.write("cs_loc_test.gpkg")
    # model.channels.write("channels_test.gpkg")

    # write model
    model.write(
        sqlite_path.replace("zandkasteel", "model_met_ingemeten_profielen_slib"),
        check=True,
        rasters=True,
    )
