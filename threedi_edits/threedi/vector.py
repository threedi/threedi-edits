# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:57:13 2021

@author: chris.kerklaan
    
"""
# Third-party imports
import logging
from osgeo import ogr

# Local imports
from threedi_edits.gis.vector import Vector
from threedi_edits.threedi.constants.constants import Properties
from threedi_edits.threedi.constants.templates import ThreediTemplateTable
from threedi_edits.gis.feature import Feature
from threedi_edits.threedi.feature import ThreediFeature

# Globals
PROPERTIES = Properties()
logger = logging.getLogger(__name__)


class ThreediVector(Vector):
    """
    The threedivector class retrieves special threedi features
    and is able to add data via a template.
    """

    instance = "threedi.vector.ThreediVector"

    def __init__(self, ogr_ds, table_name):
        Vector.__init__(self, path=None, datasource=ogr_ds, layer_name=table_name)
        self.name = table_name
        self.properties = PROPERTIES[table_name]

        if not hasattr(self, "layer"):
            logger.info(f"Table {table_name} not present.")
            return

    def __iter__(self):
        fids = [i.GetFID() for i in self.layer]
        for fid in fids:
            yield self.get_feature(fid)

    def __getitem__(self, i):
        return self.get_feature(i)

    def __repr__(self):
        fids = self.fids
        if len(fids) > 5:
            fid_print = f"{fids[:5]} ... {fids[-5:]}"
        else:
            fid_print = f"{fids}"

        return f"({self.instance}) 3Di table: {self.name}, row ids: {fid_print}"

    def get_feature(self, i):
        """helper for retrieving features"""
        return ThreediFeature(
            self.layer.GetFeature(i),
            self.layer,
            self.fields,
            self.layer_name,
        )

    def template_data(self, template: dict):
        """input is a ThreediTemplate, output is geometry, fid, attributes
        None if not present
        """
        # dervive id from template
        fid = template.get("id")
        if fid == "None":
            fid = None
        else:
            fid = int(fid)

        # derive geometry from template
        geometry = template.get("the_geom")

        # derive content from template
        attributes = {}
        for key in template.keys():
            if key in ["id", "the_geom", "table"]:
                continue
            value = template[key]
            try:
                content = value.content
            except AttributeError:
                content = value
            attributes[key] = content

        return fid, attributes, geometry

    def add(self, feature=None, geometry=None, fid=None, **items):

        # Get the last one or set to the first.

        if isinstance(feature, ThreediTemplateTable):
            fid, items, geometry = self.template_data(feature)

        elif isinstance(feature, Feature) or isinstance(feature, ogr.Feature):
            items = feature.items
            geometry = feature.geometry

        else:
            assert len(items) != 0

        if fid is None:
            if self.fids[-1] > 0:
                fid = self.fids[-1] + 1
            else:
                fid = 1

        return self._add(geometry=geometry, fid=fid, **items)
