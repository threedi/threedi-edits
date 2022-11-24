# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:52:16 2021

@author: chris.kerklaan

#TODO:
    
# The current threedi_model_checker has needs an 'migrated sqlite'. However,
when you update the migrated sqlite, you increase
    
"""

# First-party imports
import shutil
import pathlib
import logging
import types

# Third-party imports
from osgeo import ogr
from osgeo import osr

# Local imports
import threedi_edits as te
from threedi_edits.gis.vectorgroup import VectorGroup
from threedi_edits.gis.vector import Vector
from threedi_edits.threedi.vector import ThreediVector
from threedi_edits.globals import SUPPORTED_THREEDI_VERSIONS
from threedi_edits.threedi.constants import get_version

# Globals
FILE_PATH = pathlib.Path(__file__)
EMPTY_SQLITE_PATH = str(
    FILE_PATH.parent / "data" / SUPPORTED_THREEDI_VERSIONS[0] / "empty_klondike.sqlite"
)
V = get_version(SUPPORTED_THREEDI_VERSIONS[0])

# Drivers
OGR_SQLITE_DRIVER = ogr.GetDriverByName("SQLite")

# Logger
logger = logging.getLogger(__name__)


class ThreediVectorGroup(VectorGroup):
    """
    ThreediVectorGroup's makes it easier to read, write and copy threedi-objectes.
    params:
        ogr_ds: ogr.DataSource
        mode:
            'write': Used to write a file
            'empty': Returns an empty in-mem ogr threedi model
            'memory': The threedimodel is copied into memory, makes speedy edits
            'read': Used for reading a model only


    """

    instance = "threedi.vectorgroup.VectorGroup"

    def __init__(self, ogr_ds: ogr.DataSource, mode="memory"):

        self.ds = ogr_ds  # must be here to ensure correct copying.
        if mode == "write":
            super().__init__(ogr_ds=ogr_ds, memory=False, write=True)
        elif mode == "read":
            super().__init__(ogr_ds=ogr_ds, memory=False, write=False)
        elif mode == "empty":
            super().__init__(ogr_ds=create_empty_model().ds, memory=False)
        elif mode == "memory":
            super().__init__(ogr_ds=self._copy(ogr_ds), memory=False)

        logger.debug(f"Opened threedi model in {mode} mode")

        self.srs = osr.SpatialReference()
        self.srs.ImportFromEPSG(4326)
        self.mode = mode

    def __str__(self):
        return "3Di model"

    def __repr__(self):
        return f"({self.instance}) 3Di model"

    def __iter__(self):
        for table_name in V.TABLES["all"]:
            if table_name in self.layers:
                yield self._table(table_name)

    def _tables(self):
        tables = types.SimpleNamespace()
        for table in V.TABLES["all"]:
            if table in V.TRANSLATE:
                setattr(tables, V.TRANSLATE[table], self._table(table))

        return tables

    def _table(self, table_name):
        return ThreediVector(self.ds, table_name)

    def __getitem__(self, table):
        return self._table(table)

    def __setitem__(self, table_name, ogr_layer):
        """replaces an entire table"""
        self.ds.DeleteLayer(table_name)
        self.ds.CreateLayer(ogr_layer)

    def _copy(
        self, memory: bool = True, path: str = None, quiet=True
    ) -> ogr.DataSource:
        """Copies the current either to memory or to a path"""
        if memory:
            output = create_empty_model()
        else:
            output = create_output_model(path)

        current_layers = output.layers
        for table in self:
            # Some layers are remove, so we are skipping some things.
            if table.name not in current_layers:
                continue
            if table.name in V.TABLES["skip"]:
                continue

            table_output = output[table.name]
            subject = f"Copying {table.name} ({table.count} rows)"
            for feature in te.Progress(table, subject, quiet):
                table_output.add(
                    items=feature.evaluate_types(),
                    geometry=feature.geometry,
                    fid=feature.fid,
                )
            table_output = None

        return output.ds

    def _write(self, path, quiet=False):
        """
        Whilst copy layer is faster than iterating over every feature,
        we want check the feature on field types and use the sqlite constraints
        therefore iteration is chosen here
        """
        self._copy(False, path, quiet)


def create_output_model(path):
    """Creates an empty sqlite"""
    shutil.copyfile(EMPTY_SQLITE_PATH, path)
    output = ThreediVectorGroup(ogr.Open(path, 1), mode="write")
    empty_tables = [
        "v2_global_settings",
        "v2_aggregation_settings",
        "v2_numerical_settings",
        "schema_version",
    ]

    for empty_table in empty_tables:
        table = output[empty_table]
        table.table
        for fid in table.fids:
            table.delete(fid)
        table = None

    return output


def create_empty_model():
    """creates a ThreediModelBase from scratch"""
    group = VectorGroup.from_scratch(name="Threedi Model")
    for name, threedi_property in V.Properties():

        # create vector, with the right geometry and fields
        if "the_geom" in threedi_property:
            geometry_type = threedi_property["the_geom"]["type"]
        else:
            geometry_type = ogr.wkbUnknown

        vector = Vector.from_scratch(name, geometry_type, 4326)

        for i, (field_name, field_values) in enumerate(threedi_property.items()):
            if field_name in ["connected_tables", "the_geom", "id"]:
                continue
            vector.add_field(field_name, field_values["type"])

        # add to the group
        group.add(vector, name)

    return group
