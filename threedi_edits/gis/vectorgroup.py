# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:45:07 2021

@author: chris.kerklaan
"""

# First-party imports
import os
import pathlib

# Third-party imports
from osgeo import ogr

# Local imports
from .vector import Vector
from .utils import DRIVERS
from threedi_edits.utils.project import Classes, Functions

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)


# globals
OGR_MEM_DRIVER = ogr.GetDriverByName("Memory")
OPTIONS = ["OVERWRITE=YES", "FID=fid"]


class VectorGroup:
    instance = "gis.vector.VectorGroup"

    def __init__(
        self,
        path=None,
        ogr_ds=None,
        options={},
        memory=True,
        write=False,
    ):
        """vector groups"""
        if path is not None:
            self.ext = pathlib.Path(str(path)).suffix
            self.ds = ogr.Open(str(path), write, **options)

        elif ogr_ds is not None:
            self.ds = ogr_ds

        if memory:
            self.ds = self._copy()

    @classmethod
    def from_pg(cls, host, port, user, password, dbname, write=False, memory=False):
        pg_string = ("PG:host={} port={} user='{}'" "password='{}' dbname='{}'").format(
            host, port, user, password, dbname
        )

        return cls(ogr_ds=ogr.Open(pg_string, write), memory=memory)

    @classmethod
    def from_scratch(cls, name="scratch"):
        return cls(ogr_ds=OGR_MEM_DRIVER.CreateDataSource(name), memory=False)

    def __getitem__(self, i):
        return Vector.from_ds(self.ds, layer_name=i)

    def __iter__(self):
        for layer in self.layers:
            yield self[layer]

    @property
    def layers(self):
        return [layer.GetName() for layer in self.ds]

    @property
    def drivers(self):
        return DRIVERS

    def add(self, vector, name):
        vector.layer_name = name
        self.ds.CopyLayer(vector.layer, name, options=["OVERWRITE=YES"])

    def add_styling(self, vector, name="layer_styles"):
        """add styling to the gpkg by adding the layer_styles column"""
        self.add(vector, name)

    def clear(self):
        for table in self:
            table.delete_all()

    def _copy(self):
        return OGR_MEM_DRIVER.CopyDataSource(self.ds, "")

    def write(self, path, overwrite=False):
        path = str(path)  # pathlib

        driver_name = self.drivers[os.path.splitext(path)[-1]]
        driver = ogr.GetDriverByName(driver_name)

        if overwrite:
            if os.path.exists(path):
                driver.DeleteDataSource(path)

        # create output file
        out_ds = driver.CreateDataSource(path)

        for layer in self:
            if "fid" in layer.fields:
                layer.delete_field("fid")
            out_ds.CopyLayer(layer.layer, layer.name, options=OPTIONS)

        out_ds = None

    def execute_sql(self, sql, dialect=None):
        self.ds.ExecuteSQL(sql, dialect)

    def clip(self, vector):
        output = VectorGroup.from_scratch()
        for layer in self:
            output.add(layer.clip(vector), layer.name)
        return output

    def __repr__(self):
        return str(self.layers)
