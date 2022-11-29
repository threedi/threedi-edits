# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:01:54 2021

@author: chris.kerklaan

# settings: https://github.com/nens/threedi-turtle-scripts/blob/master/sql-scripts/scripts-2B/process-3di/settings_for_2d_model.sql
"""
# system imports
import copy

# Local imports
from .constants import Properties

PROPERTIES = Properties()


template_integrated = {
    "global_settings": {
        # General
        "id": None,
        "name": "test",
        "use_2d_flow": 1,
        "use_1d_flow": 1,
        "use_0d_inflow": 1,
        "use_2d_rain": 1,
        # Grid
        "grid_space": 10,
        "kmax": 1,
        "table_step_size": 0.01,
        # Terrain info
        # Dem
        "dem_file": None,
        "epsg_code": 28992,
        # Friction
        "frict_coef_file": None,
        "frict_coef": 0.026,
        "frict_type": 2,
        "frict_avg": 0,
        # Groundwater
        "initial_groundwater_level_file": None,
        "initial_groundwater_level": None,
        "initial_groundwater_level_type": None,
        # Initial waterlevel
        "initial_waterlevel_file": None,
        "initial_waterlevel": -10,
        "water_level_ini_type": None,
        # Interception
        "interception_file": None,
        "interception_global": None,
        # Wind
        "wind_shielding_file": None,
        # Time
        "start_date": "2014-05-14",
        "start_time": "2014-05-14 17:42:18",
        "sim_time_step": 5,
        "timestep_plus": True,
        "minimum_sim_time_step": 1,
        "maximum_sim_time_step": 60,
        "nr_timesteps": 2880,
        "output_time_step": 1,
        # Settings id's
        "interflow_settings_id": None,
        "groundwater_settings_id": None,
        "numerical_settings_id": 1,
        "simple_infiltration_settings_id": 1,
        "control_group_id": None,
        # Extra options 1D
        "advection_1d": 0,
        "dist_calc_points": 1000,
        "manhole_storage_area": 100,
        "max_angle_1d_advection": 90,
        "table_step_size_1d": None,
        "advection_2d": 0,
        "dem_obstacle_detection": False,
        "guess_dams": False,
        "dem_obstacle_height": None,
        "embedded_cutoff_threshold": None,
        "flooding_threshold": 0.001,
        "table_step_size_volume_2d": None,
    },
    "numerical_settings": {
        # General
        "id": None,
        "limiter_grad_1d": None,
        "limiter_grad_2d": None,
        "limiter_slope_crossectional_area_2d": None,
        "limiter_slope_friction_2d": None,
        # Matrix
        "convergence_cg": None,
        "convergence_eps": 0.00005,
        "use_of_cg": 20,
        "use_of_nested_newton": 1,
        "max_degree": 7,
        "max_nonlin_iterations": 20,
        "precon_cg": 1,
        "integration_method": 0,
        # Thresholds
        "flow_direction_threshold": None,
        "general_numerical_threshold": None,
        "thin_water_layer_definition": None,
        "minimum_friction_velocity": None,
        "minimum_surface_area": None,
        # Miscellaneous
        "cfl_strictness_factor_1d": 1,
        "cfl_strictness_factor_2d": 1,
        "frict_shallow_water_correction": 0,
        "pump_implicit_ratio": 1,
        "preissmann_slot": None,
    },
    "simple_infiltration": {
        "id": None,
        "display_name": "Infiltration",
        "infiltration_rate": 0.01,
        "infiltration_rate_file": None,
        "max_infiltration_capacity_file": None,
        "infiltration_surface_option": 0,
    },
    "connection_nodes": {
        "id": None,
        "code": None,
        "initial_waterlevel": None,
        "storage_area": 0.64,
        "the_geom": None,
    },
    "manhole": {
        # general
        "id": None,
        "display_name": None,
        "code": None,
        "connection_node_id": None,
        "calculation_type": None,
        "bottom_level": None,
        "calculation_type": None,
        "drain_level": None,
        "length": None,
        "width": None,
        "zoom_category": 1,
        "manhole_indicator": None,
        "shape": None,
        "surface_level": None,
        "sediment_level": None,
    },
    "channel": {
        "id": None,
        "calculation_type": 102,
        "code": None,
        "connection_node_end_id": None,
        "connection_node_start_id": None,
        "display_name": None,
        "dist_calc_points": 10,
        "zoom_category": 5,
        "the_geom": None,
    },
    "cross_section_definition": {
        "id": None,
        "code": None,
        "height": None,
        "shape": None,
        "width": None,
    },
    "cross_section_location": {
        "id": None,
        "bank_level": None,
        "channel_id": None,
        "code": None,
        "definition_id": None,
        "friction_type": 2,
        "friction_value": 0.026,
        "reference_level": None,
        "the_geom": None,
    },
    "grid_refinement_area": {
        "id": None,
        "display_name": None,
        "refinement_level": None,
        "the_geom": None,
        "code": None,
    },
    "grid_refinement": {
        "id": None,
        "display_name": None,
        "refinement_level": None,
        "the_geom": None,
        "code": None,
    },
    "1d_boundary_conditions": {
        "id": None,
        "boundary_type": None,
        "connection_node_id": None,
        "timeseries": None,
    },
    "pumpstation": {
        "id": None,
        "capacity": None,
        "zoom_category": 5,
        "code": None,
        "connection_node_end_id": None,
        "connection_node_start_id": None,
        "classification": None,
        "display_name": None,
        "lower_stop_level": None,
        "sewerage": False,
        "start_level": None,
        "type": 1,
        "upper_stop_level": None,
    },
    "levee": {
        "id": None,
        "crest_level": None,
        "material": 2,
        "max_breach_depth": 10,
        "code": None,
        "the_geom": None,
    },
    "obstacle": {
        "id": None,
        "crest_level": None,
        "code": None,
        "the_geom": None,
    },
    "interflow": {
        "id": None,
        "display_name": None,
        "interflow_type": None,
        "porosity": None,
        "porosity_file": 10,
        "porosity_layer_thickness": None,
        "hydraulic_conductivity_file": None,
        "hydraulic_conductivity": None,
        "impervious_layer_elevation": None,
    },
    "orifice": {
        "display_name": None,
        "code": None,
        "zoom_category": None,
        "connection_node_start_id": None,
        "connection_node_end_id": None,
        "crest_level": None,
        "crest_type": 4,
        "cross_section_definition_id": None,
        "discharge_coefficient_negative": 0.8,
        "discharge_coefficient_positive": 0.8,
        "friction_type": 2,
        "friction_value": 0.011,
        "sewerage": 1,
    },
}

template_2d_model = {
    "global_settings": {
        # General
        "id": None,
        "name": "test",
        "use_2d_flow": 1,
        "use_1d_flow": 0,
        "use_0d_inflow": 0,
        "use_2d_rain": 1,
        # Grid
        "grid_space": 10,
        "kmax": 1,
        "table_step_size": 0.01,
        # Terrain info
        # Dem
        "dem_file": "rasters/dem.tif",
        "epsg_code": 28992,
        # Friction
        "frict_coef_file": None,
        "frict_coef": 0.026,
        "frict_type": 2,
        "frict_avg": 0,
        # Groundwater
        "initial_groundwater_level_file": None,
        "initial_groundwater_level": None,
        "initial_groundwater_level_type": None,
        # Initial waterlevel
        "initial_waterlevel_file": None,
        "initial_waterlevel": -10,
        "water_level_ini_type": None,
        # Interception
        "interception_file": None,
        "interception_global": None,
        # Wind
        "wind_shielding_file": None,
        # Time
        "start_date": "2014-05-14",
        "start_time": "2014-05-14 17:42:18",
        "sim_time_step": 30,
        "timestep_plus": True,
        "minimum_sim_time_step": 30,
        "maximum_sim_time_step": 300,
        "nr_timesteps": 240,
        "output_time_step": 300,
        # Settings id's
        "interflow_settings_id": None,
        "groundwater_settings_id": None,
        "numerical_settings_id": 1,
        "simple_infiltration_settings_id": 1,
        "control_group_id": None,
        # Extra options 1D
        "advection_1d": 0,
        "dist_calc_points": 1000,
        "manhole_storage_area": 0,
        "max_angle_1d_advection": 90,
        "table_step_size_1d": None,
        "advection_2d": 0,
        "dem_obstacle_detection": False,
        "guess_dams": False,
        "dem_obstacle_height": None,
        "embedded_cutoff_threshold": None,
        "flooding_threshold": 0.01,
        "table_step_size_volume_2d": None,
    },
    "numerical_settings": {
        "id": None,
        # General
        "limiter_grad_1d": None,
        "limiter_grad_2d": None,
        "limiter_slope_crossectional_area_2d": None,
        "limiter_slope_friction_2d": None,
        # Matrix
        "convergence_cg": None,
        "convergence_eps": 0.00001,
        "use_of_cg": 20,
        "use_of_nested_newton": 0,
        "max_degree": 5,
        "max_nonlin_iterations": 20,
        "precon_cg": 1,
        "integration_method": 0,
        # Thresholds
        "flow_direction_threshold": None,
        "general_numerical_threshold": None,
        "thin_water_layer_definition": None,
        "minimum_friction_velocity": None,
        "minimum_surface_area": None,
        # Miscellaneous
        "cfl_strictness_factor_1d": 1,
        "cfl_strictness_factor_2d": 1,
        "frict_shallow_water_correction": 0,
        "pump_implicit_ratio": 1,
        "preissmann_slot": None,
    },
    "simple_infiltration": {
        "id": None,
        "display_name": "Infiltration",
        "infiltration_rate": 0,
        "infiltration_rate_file": None,
        "max_infiltration_capacity_file": None,
        "infiltration_surface_option": 0,
    },
    "connection_nodes": {
        "id": None,
        "code": None,
        "initial_waterlevel": None,
        "storage_area": 0.64,
        "the_geom": None,
    },
    "manhole": {
        "id": None,
        # general
        "display_name": None,
        "code": None,
        "connection_node_id": None,
        "calculation_type": None,
    },
    "grid_refinement": {
        "id": None,
        "display_name": None,
        "refinement_level": None,
        "the_geom": None,
        "code": None,
    },
    "grid_refinement_area": {
        "id": None,
        "display_name": None,
        "refinement_level": None,
        "the_geom": None,
    },
    "1d_boundary_conditions": {
        "id": None,
        "boundary_type": None,
        "connection_node_id": None,
        "timeseries": None,
    },
    "levee": {
        "id": None,
        "crest_level": None,
        "material": None,
        "code": None,
        "max_breach_depth": 10,
        "the_geom": None,
    },
    "obstacle": {
        "id": None,
        "crest_level": None,
        "code": None,
        "the_geom": None,
    },
    "interflow": {
        "id": None,
        "display_name": None,
        "interflow_type": None,
        "porosity": None,
        "porosity_file": 10,
        "porosity_layer_thickness": None,
        "hydraulic_conductivity_file": None,
        "hydraulic_conductivity": None,
        "impervious_layer_elevation": None,
    },
    "orifice": {
        "display_name": None,
        "code": None,
        "zoom_category": None,
        "connection_node_start_id": None,
        "connection_node_end_id": None,
        "crest_level": None,
        "crest_type": 4,
        "cross_section_definition_id": None,
        "discharge_coefficient_negative": 0.8,
        "discharge_coefficient_positive": 0.8,
        "friction_type": 2,
        "friction_value": 0.011,
        "sewerage": 1,
    },
}


template_mapping = {
    "integrated": copy.deepcopy(template_integrated),
    "2d_model": copy.deepcopy(template_2d_model),
}


class ThreediString(str):
    """
    Class acts (therefore we use __new__) as a string, with multiple inputs
    Class has info (__init__)
    Class prints info

    """

    name = "ThreediString"

    def __new__(cls, content, info, table):
        return super().__new__(cls, content)

    def __init__(self, content, info, table):
        self.table = table
        self.finfo = info
        self.content = content

    @property
    def info(self):
        print(self.finfo)

    def __str__(self):
        if self.info is not None:
            return f"Value: {self.content}\nDescription: {self.info}"
        elif self.content is not None:
            return str(self.content)
        else:
            return "(ThreediString) Template not filled"

    def __repr__(self):
        if self.info is not None:
            return f"Value: {self.content}\nDescription: {self.info}"
        elif self.content is not None:
            return str(self.content)
        else:
            return "(ThreediString) Template not filled"


class ThreediInteger(int):
    """integer or bool"""

    name = "ThreediInteger"

    def __new__(cls, content, info, table):
        return super().__new__(cls, content)

    def __init__(self, content, info, table):
        self.table = table
        self.finfo = info
        self.content = content

    @property
    def info(self):
        print(self.finfo)

    def __str__(self):
        if self.info is not None:
            return f"Value: {self.content}\nDescription: {self.info}"
        elif self.content is not None:
            return str(self.content)
        else:
            return "(ThreediInteger) Template not filled"

    def __repr__(self):
        if self.info is not None:
            return f"Value: {self.content}\nDescription: {self.info}"
        elif self.content is not None:
            return str(self.content)
        else:
            return "(ThreediInteger) Template not filled"


class ThreediFloat(float):
    """Float"""

    name = "ThreediFloat"

    def __new__(cls, content, info, table):
        return super().__new__(cls, content)

    def __init__(self, content, info, table):
        self.table = table
        self.finfo = info
        self.content = content

    @property
    def info(self):
        print(self.finfo)

    def __str__(self):
        if self.info is not None:
            return f"Value: {self.content}\nDescription: {self.info}"
        elif self.content is not None:
            return str(self.content)
        else:
            return "(ThreediFloat) Template not filled"

    def __repr__(self):
        if self.info is not None:
            return f"Value: {self.content}\nDescription: {self.info}"
        elif self.content is not None:
            return str(self.content)
        else:
            return "(ThreediFloat) Template not filled"


def threedi_value(value, info, table):
    value_type = type(value)
    if value_type == str or value_type == type(None):
        return ThreediString(value, info, table)
    elif value_type == int or value_type == bool:
        return ThreediInteger(value, info, table)
    elif value_type == float:
        return ThreediFloat(value, info, table)
    else:
        return value


class ThreediTemplateTable(dict):
    """table template"""

    def __init__(self, in_dict):
        self.__dict__ = in_dict

    def __setitem__(self, key, value):
        d = dict(self.__dict__)
        table = self.table
        # if type(d[key]) != type(value):
        #     raise TypeError(f"Type mismatch {type(d[key]) } {type(value)}")

        d[key] = threedi_value(
            value, PROPERTIES["v2_" + table][key]["description"], table
        )
        self.__dict__ = d

    def __getitem__(self, key):
        return self.__dict__[key]

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def items(self):
        for key, value in self.__dict__.items():
            if not key == "table":
                yield key, value

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def get(self, key):
        try:
            value = self[key]
        except KeyError:
            value = None
        return value


class Templates(object):
    def __init__(self, template_type="integrated"):
        self.type = template_type
        self.template = template_mapping[self.type]

    def get_template(self, table):
        values = self.template[table]
        new_dict = {}
        for key, value in values.items():
            v = threedi_value(
                value, PROPERTIES["v2_" + table][key]["description"], table
            )
            new_dict[key] = v

        table_template = ThreediTemplateTable(new_dict)
        table_template.table = table
        return table_template

    @property
    def global_setting(self):
        return self.get_template("global_settings")

    @property
    def numerical_setting(self):
        return self.get_template("numerical_settings")

    @property
    def simple_infiltration(self):
        return self.get_template("simple_infiltration")

    @property
    def node(self):
        return self.get_template("connection_nodes")

    @property
    def manhole(self):
        return self.get_template("manhole")

    @property
    def channel(self):
        return self.get_template("channel")

    @property
    def cross_section(self):
        return self.get_template("cross_section_definition")

    @property
    def cross_section_location(self):
        return self.get_template("cross_section_location")

    @property
    def grid_refinement_area(self):
        return self.get_template("grid_refinement_area")

    @property
    def grid_refinement(self):
        return self.get_template("grid_refinement")

    @property
    def boundary_condition_1d(self):
        return self.get_template("1d_boundary_conditions")

    @property
    def pumpstation(self):
        return self.get_template("pumpstation")

    @property
    def levee(self):
        return self.get_template("levee")

    @property
    def interflow(self):
        return self.get_template("interflow")

    @property
    def obstacle(self):
        return self.get_template("obstacle")

    @property
    def orifice(self):
        return self.get_template("orifice")
