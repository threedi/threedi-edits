# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 15:34:04 2022

@author: chris.kerklaan

@docs: Command line tool used for the downloads and conversion of 3Di rasters on a linux server.

"""
# find the package
import sys

sys.path.append(r"C:\Users\chris.kerklaan\Documents\Github\threedi-raster-edits")

# First imports
import pathlib
import threedi_edits as tre
import argparse

# Globals
CELLSIZE = 0.5
THREADS = 1
FIELD_NAME = None
AHN_VERSION = "ahn4"


def rextract_convert(
    vector_path: str,
    lizard_api_key: str,
    output_directory,
    field_name: str,
    cellsize,
    threads: int,
    friction: bool,
    infiltration: bool,
    interception: bool,
    landuse: bool,
    soil: bool,
    ahn: str,
):

    # create a dict
    threedi_uuids = {}
    if landuse:
        threedi_uuids["landuse"] = tre.UUID.THREEDI_LANDUSE

    if soil:
        threedi_uuids["soil"] = tre.UUID.THREEDI_SOIL

    get_ahn = {
        "ahn4": tre.UUID.THREEDI_AHN4,
        "ahn3": tre.UUID.THREEDI_AHN3,
        "ahn4_ahn3_combined": "6a75a859-e644-4498-a4c6-dbf367464ffd",
    }

    # select correct version
    threedi_uuids["ahn"] = get_ahn[ahn.lower()]

    # load you vector
    vector = tre.Vector(vector_path)

    # Create a rextract object
    rextract = tre.RasterExtraction(lizard_api_key)

    if field_name:
        if not field_name in vector.fields:
            raise ValueError(f"Field name {field_name} not found in {vector.fields}")

    for feature in vector:

        # retrieve the correct field name
        if field_name:
            output_folder_name = feature[field_name]
        else:
            output_folder_name = f"feature_{feature.fid}"

        # check if the folder exists, create otherwise
        folder = pathlib.Path(output_directory + "/" + output_folder_name)
        folder.mkdir(parents=True, exist_ok=True)

        # do download, store locations
        threedi_downloads = {}
        for name, threedi_uuid in threedi_uuids.items():
            path = str(folder / f"{name}.tif")
            rextract.run(
                path, threedi_uuid, feature.geometry, cellsize=cellsize, threads=threads
            )
            threedi_downloads[name] = path

        # do conversion
        if landuse and soil:
            group = tre.ThreediRasterGroup(
                # dem_file=threedi_downloads["ahn"],
                landuse_file=threedi_downloads["landuse"],
                soil_file=threedi_downloads["soil"],
                memory=False,
            )

            # load table
            group.load_landuse_conversion_table()
            group.load_soil_conversion_table()

            if infiltration:
                group.generate_infiltration()  # soil, landuse
                group.infiltration.spatial_reference = 28992  # TODO
                group.infiltration.write(str(folder / "infiltration.tif"))

            if friction:
                group.generate_friction()  # landuse
                group.friction.spatial_reference = 28992  # TODO
                group.friction.write(str(folder / "friction.tif"))

            if interception:
                group.generate_interception()  # landuse
                group.interception.spatial_reference = 28992  # TODO
                group.interception.write(str(folder / "interception.tif"))


def get_parser():
    """Return argument parser."""

    parser = argparse.ArgumentParser(description=__doc__)

    # main
    parser.add_argument("vector_path", help="Features of this vector are downloaded.")
    parser.add_argument("lizard_api_key", help="Lizard api key.")
    parser.add_argument("output_directory", help="Output rasters are placed here.")

    # options
    parser.add_argument(
        "-f", "--friction", action="store_true", help="Generates friction."
    )
    parser.add_argument(
        "-ic", "--interception", action="store_true", help="Generates interception"
    )
    parser.add_argument(
        "-if", "--infiltration", action="store_true", help="Generates infiltration"
    )
    parser.add_argument("-soil", "--soil", action="store_true", help="Downloads soil")
    parser.add_argument(
        "-landuse", "--landuse", action="store_true", help="Downloads landuse"
    )
    parser.add_argument(
        "-fn",
        "--field_name",
        default=FIELD_NAME,
        help="Field name used to place the rasters in a directory.",
    )
    parser.add_argument(
        "-cs", "--cellsize", default=CELLSIZE, type=float, help="Cellsize."
    )
    parser.add_argument(
        "-tr",
        "--threads",
        default=THREADS,
        type=int,
        help="Amount of threads used. More threads increases speed. Cannot be above 4.",
    )

    parser.add_argument(
        "-ahn",
        "--ahn",
        default=AHN_VERSION,
        type=str,
        help="AHN version: ahn4, ahn3 or ahn4_ahn3_combined",
    )

    parser.set_defaults(friction=False)
    parser.set_defaults(interception=False)
    parser.set_defaults(infiltration=False)
    parser.set_defaults(soil=False)
    parser.set_defaults(landuse=False)

    return parser


def main():
    """Call extract_all with args from parser."""
    return rextract_convert(**vars(get_parser().parse_args()))


if __name__ == "__main__":
    exit(main())
