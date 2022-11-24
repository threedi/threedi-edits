# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:45:53 2021

@author: chris.kerklaan

"""
# Third-party imports
from osgeo import ogr, _ogr

# ogr type mapping
OGR_TYPES = {
    bool: ogr.OFTBinary,
    int: ogr.OFTInteger,
    float: ogr.OFTReal,
    str: ogr.OFTString,
}

OGR_TYPES_INVERT = {
    ogr.OFTInteger: int,  # 0
    ogr.OFSTBoolean: bool,  # 1
    ogr.OFTReal: float,  # 2
    ogr.OFTRealList: float,  # 3
    ogr.OFTString: str,  # 4
    ogr.OFTStringList: str,
    ogr.OFTWideString: str,  # 6
    ogr.OFTBinary: bool,  # 8
    ogr.OFTDate: str,  # 9
    ogr.OFTTime: str,  # 10
    ogr.OFTDateTime: str,  # 11
    ogr.OFTInteger64: int,  # 12
    ogr.OFTInteger64List: int,  # 13
}
GET_FIELD = {
    ogr.OFTInteger: _ogr.Feature_GetFieldAsInteger,  # 0
    ogr.OFSTBoolean: _ogr.Feature_GetFieldAsInteger,  # 1
    ogr.OFTReal: _ogr.Feature_GetFieldAsDouble,  # 2
    ogr.OFTRealList: _ogr.Feature_GetFieldAsDouble,  # 3
    ogr.OFTString: _ogr.Feature_GetFieldAsString,  # 4
    ogr.OFTStringList: _ogr.Feature_GetFieldAsString,
    ogr.OFTWideString: _ogr.Feature_GetFieldAsString,  # 6
    ogr.OFTBinary: _ogr.Feature_GetFieldAsInteger,  # 8
    ogr.OFTDate: _ogr.Feature_GetFieldAsString,  # 9
    ogr.OFTTime: _ogr.Feature_GetFieldAsString,  # 10
    ogr.OFTDateTime: _ogr.Feature_GetFieldAsString,  # 11
    ogr.OFTInteger64: _ogr.Feature_GetFieldAsInteger,  # 12
    ogr.OFTInteger64List: _ogr.Feature_GetFieldAsInteger,  # 13
}


class Fields:
    """ogr fields has to act like a dicionary"""

    def __init__(self, layer):
        self.layer = layer

    def __getitem__(self, field_name):
        return self.type(field_name)

    def __setitem__(self, field_name, field_type):
        self.set_type(field_name, field_type)

    def __repr__(self):
        return str(self.items)

    def __contains__(self, field_name):
        return field_name in self.items

    def __dict__(self):
        return self.items

    @property
    def defn(self):
        if not hasattr(self, "_defn"):
            self._defn = _ogr.Layer_GetLayerDefn(self.layer)
        return self._defn

    @property
    def names(self):
        if not hasattr(self, "_names"):
            self._names = get_field_names(self.defn)
        return self._names

    @property
    def types(self):
        if not hasattr(self, "_types"):
            self._types = get_field_types(self.defn)

        return self._types

    @property
    def nullables(self):
        if not hasattr(self, "_nullables"):
            self._nullables = get_is_nullables(self.defn)
        return self._nullables

    @property
    def type_names(self):
        if not hasattr(self, "_type_names"):
            self._type_names = [ogr.GetFieldTypeName(n) for n in self.types]
        return self._type_names

    @property
    def items(self):
        if not hasattr(self, "_items"):
            self._items = {
                n: [t, tn, ll]
                for n, t, tn, ll in zip(
                    self.names, self.types, self.type_names, self.nullables
                )
            }
        return self._items

    @property
    def keys(self):
        if not hasattr(self, "_keys"):
            self._keys = list(self.items.keys())
        return self._keys

    @property
    def indices(self):
        if not hasattr(self, "_indices"):
            self._indices = {
                name: self.layer.FindFieldIndex(name, 0) for name in self.names
            }
        return self._indices

    def get_field_func(self, field_name):
        """returns a function and an index"""
        return GET_FIELD[self.items[field_name][0]]

    def add(self, field_name, ftype, has_none=True):
        if type(ftype) != int:
            ftype = OGR_TYPES[ftype]

        if field_name not in self.names:
            defn = ogr.FieldDefn(field_name, ftype)
            defn.SetNullable(has_none)
            self.layer.CreateField(defn)

    def delete(self, field_name):
        self.layer.DeleteField(self.layer.FindFieldIndex(field_name, 0))

    def set_name(self, old_field_name, new_field_name):
        i = self.defn.GetFieldIndex(old_field_name)
        field_defn = self.layer.GetFieldDefn(i)
        field_defn.SetName(new_field_name)
        self.layer.AlterFieldDefn(i, field_defn, ogr.ALTER_NAME_FLAG)

    def set_type(self, field_name: str, field_type: type):
        """field type is change by changing the data"""
        ogr_type = OGR_TYPES[field_type]
        data = []

        for feature in self.layer:
            data.append(feature[field_name])

        self.layer.DeleteField(self.layer.FindFieldIndex(field_name, 0))

        defn = ogr.FieldDefn(field_name, ogr_type)
        self.layer.CreateField(defn)

        self.layer.ResetReading()
        for i, feature in enumerate(self.layer):
            value = data[i]
            if value:
                value = field_type(value)

            feature[field_name] = value
            self.layer.SetFeature(feature)


class GeometryFields:
    def __init__(self, feature):
        self.feature = feature

    def __getitem__(self, key):
        return self.items[key]

    @property
    def count(self):
        return _ogr.Feature_GetGeomFieldCount(self.feature)

    @property
    def translate(self):
        return {
            self.ogr_feature.GetGeomFieldDefnRef(i).name: i for i in range(self.count)
        }

    @property
    def names(self):
        return [self.ogr_feature.GetGeomFieldDefnRef(i).name for i in range(self.count)]

    @property
    def items(self):
        geom_items = {}
        for i in range(self.count):
            geometry = self.ogr_feature.GetGeomFieldRef(i)
            geom_items[self.ogr_feature.GetGeomFieldDefnRef(i).name] = geometry
        return geom_items


def get_field_defn(defn):
    return [_ogr.FeatureDefn_GetFieldDefn(defn, i) for i in range(defn.GetFieldCount())]


def get_field_names(defn):
    """returns of list of column names for this field"""
    return [_ogr.FieldDefn_GetName(i) for i in get_field_defn(defn)]


def get_field_types(defn):
    """returns of list of column names for this field"""
    return [_ogr.FieldDefn_GetType(i) for i in get_field_defn(defn)]


def get_is_nullables(defn):
    return [_ogr.FieldDefn_IsNullable(i) for i in get_field_defn(defn)]
