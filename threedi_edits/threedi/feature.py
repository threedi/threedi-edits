# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:58:57 2021

@author: chris.kerklaan

We wrap a feature for threedimodels here.

"""
# system imports
import logging

# Third-party import
import warnings
from osgeo import ogr

# Local imports
from threedi_edits.gis.feature import Feature
from threedi_edits.gis.fields import OGR_TYPES_INVERT, Fields
from threedi_edits.threedi.utils.exceptions import ThreediValueTypeError
from threedi_edits.threedi.utils.exceptions import ThreediValueTypeWarning
from threedi_edits.threedi.utils.exceptions import _HAPPY
from threedi_edits.threedi.utils.exceptions import custom_formatwarning

from threedi_edits.globals import SUPPORTED_THREEDI_VERSIONS
from threedi_edits.threedi.constants import get_version

V = get_version(SUPPORTED_THREEDI_VERSIONS[0])

# Globals
PROPERTIES = V.Properties()
warnings.formatwarning = custom_formatwarning


logger = logging.getLogger(__name__)


class ThreediFeature(Feature):

    instance = "threedi.feature.ThreediFeature"

    def __init__(
        self,
        feature: ogr.Feature,
        layer: ogr.Layer,
        fields: Fields,
        table_name: str,
        custom_set_function=None,
    ):
        Feature.__init__(self, feature, layer, fields)
        self.table_name = table_name
        self.custom_set_function = custom_set_function
        self.properties = PROPERTIES[self.table_name]

    def __str__(self):
        return f"3Di row, items: {self.all_items}"

    def __repr__(self):
        return f"({self.instance}) 3Di row, items: {self.all_items}"

    def __getitem__(self, field_name):
        return self.correct_type(field_name, self.get_item(field_name))

    def __setitem__(self, field_name, value):
        """setitem uses a function which is define by the model below"""
        if self.custom_set_function:
            key, value = self.custom_set_function(field_name, value)
        self.set_item(field_name, value)

    def set_geometry(self, geometry):
        """sets and converts the geometry"""
        if "the_geom" in self.geometry_fields.items:
            field_id = self.geometry_fields.translate["the_geom"]
            self.feature.SetGeomField(field_id, geometry)
        else:
            self.feature.SetGeometry(geometry)

        self.layer.SetFeature(self.feature)

    def set_item(self, field_name, value):
        """
        General routine for __setitem__,
        specific routines are added in the models themselves
        """

        # set the geometry
        if field_name in ["geometry", "the_geom"]:
            self.set_geometry(value)
            return

        self.feature.SetField(field_name, self.correct_type(field_name, value))
        self.layer.SetFeature(self.feature)

    def evaluate_types(self) -> dict:
        """
        Checks and changes field types are corresponding to 3Di types.
        In memory types can be slightley adjusted.
        returns a new dictionary.
        """

        return {k: self.correct_type(k, v) for k, v in self}

    def correct_type(self, field_name, value):
        """
        Loading in memory causes some type changes. Therefore we need to
        correct the types to get a good 3Di model.

        If the field is optional and is an 'id' field and is 0
        it should be None.

        Other types are compared with a standard types in constants.py

        Returns a correct value.
        """

        # This is what we should have.
        output_type = OGR_TYPES_INVERT[self.properties[field_name]["type"]]

        # First we check if the value we have should be None.
        # Integers
        is_optional = self.properties[field_name]["optional"]
        if is_optional and type(value) == int and value == 0:
            return None

        # floats # TODO this is pretty tricky since sometimes we want 0.0?
        if is_optional and type(value) == float and value == 0.0:
            return None

        # strings
        if is_optional and type(value) == str and value == "None":
            return None

        # Second we check if None is allowed.
        if value == None:
            if is_optional == True:  # none allowed
                return value
            else:
                raise ThreediValueTypeError(
                    "Value None is not allowed here",
                    self.table_name,
                    field_name,
                    value,
                )

        # Third we correct to a string, integer, float or bool.
        type_value = type(value)
        correct_type = output_type == type(value)

        if correct_type:
            return value

        if not correct_type:
            logger.debug(
                ThreediValueTypeWarning(
                    f"incorrect type, trying to correct {type_value} to {output_type}",
                    self.table_name,
                    field_name,
                    value,
                )
            )
            try:
                if output_type == bool:
                    value = int(value)
                value = output_type(value)
            except Exception:
                raise ThreediValueTypeError(
                    f"""Found incorrect type of {type_value} where  {output_type} is expected""",
                    self.table_name,
                    field_name,
                    value,
                )
            else:
                logger.debug(
                    ThreediValueTypeWarning(
                        f"{_HAPPY}, corrected input to {type(value)}: {value}",
                        self.table_name,
                        field_name,
                        value,
                    )
                )
                return value

    # lines = ['{']
    # for key, value in hoorn.global_setting.items.items():
    #     lines.append('{}:{}'.format(key, value))
    # lines.append(['}'])
