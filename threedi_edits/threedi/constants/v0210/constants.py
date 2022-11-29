# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:32:51 2020

@author: chris.kerklaan

Version 0209.

"""
from osgeo import ogr

TABLES = {
    "all": [
        "v2_simple_infiltration",
        "v2_global_settings",
        "v2_aggregation_settings",
        "v2_interflow",
        "v2_groundwater",
        "v2_1d_boundary_conditions",
        "v2_2d_boundary_conditions",
        "v2_1d_lateral",
        "v2_2d_lateral",
        "v2_calculation_point",
        "v2_connection_nodes",
        "v2_cross_section_location",
        "v2_culvert",
        "v2_channel",
        "v2_connected_pnt",
        "v2_pipe",
        "v2_weir",
        "v2_orifice",
        "v2_manhole",
        "v2_pumpstation",
        "v2_windshielding",
        "v2_cross_section_definition",
        "v2_obstacle",
        "v2_levee",
        "v2_grid_refinement",
        "v2_grid_refinement_area",
        "v2_numerical_settings",
        "v2_dem_average_area",
        "v2_impervious_surface_map",
        "v2_impervious_surface",
        "v2_surface_parameters",
        "v2_surface_map",
        "v2_surface",
        "v2_control",
        "v2_control_delta",
        "v2_control_group",
        "v2_control_measure_group",
        "v2_control_measure_map",
        "v2_control_memory",
        "v2_control_pid",
        "v2_control_table",
        "v2_control_timed",
        "v2_floodfill",
        "v2_pumped_drainage_area",
        "v2_initial_waterlevel",
        "schema_version",
    ],
    "csd": ["v2_pipe", "v2_weir", "v2_culvert"],
    "geom": [
        "v2_connection_nodes",
        "v2_obstacle",
        "v2_culvert",
        "v2_channel",
        "v2_cross_section_location",
    ],
    "order": [
        "v2_pipe",
        "v2_weir",
        "v2_culvert",
        "v2_pumpstation",
        "v2_orifice",
        "v2_cross_section_location",
        "v2_channel",
        "v2_impervious_surface_map",
        "v2_impervious_surface",
        "v2_manhole",
        "v2_1d_boundary_conditions",
        "v2_connection_nodes",
        "v2_cross_section_definition",
        "v2_grid_refinement_area",
        "v2_obstacle",
        "v2_aggregation_settings",
        "v2_global_settings",
        "v2_numerical_settings",
        "v2_simple_infiltration",
    ],
    "single": [
        "v2_manhole",
        "v2_impervious_surface_map",
        "v2_1d_boundary_conditions",
    ],
    "start_end": [
        "v2_weir",
        "v2_pipe",
        "v2_culvert",
        "v2_pumpstation",
        "v2_orifice",
        "v2_channel",
    ],
    "skip": ["v2_aggregation_settings"],
}

raster_fields = [
    "dem_file",
    "frict_coef_file",
    "max_interception_file",
    "interception_file",
    "initial_waterlevel_file",
    "wind_shielding_file",
]

extra_sqlite_tables = [
    "layer_statistics",
    "views_layer_statistics",
    "spatialite_history",
    "views_geometry_columns",
]


connected_tables = {
    "global_settings": ["numerical_settings"],
    "boundary_conditions": [],
    "1d_lateral": ["connection_nodes"],
    "2d_boundary_conditions": [],
    "2d_lateral": [],
    "aggregation_settings": ["global_settings"],
    "channel": ["connection_nodes"],
    "connection_nodes": [],
    "cross_section_definition": [],
    "cross_section_location": ["channel", "cross_section_definition"],
    "culvert": ["cross_section_definition", "connection_nodes"],
    "numercial_settings": [],
    "grid_refinement": [],
    "grid_refinement_area": [],
    "impervious_surface": [],
    "impervious_surface_map": ["connection_nodes", "impervious_surface"],
    "levee": [],
    "manhole": ["connection_nodes"],
    "obstacle": [],
    "orifice": ["cross_section_definition", "connection_nodes"],
    "pipe": ["cross_section_definition", "connection_nodes"],
    "pumpstation": ["connection_nodes"],
    "surface": ["surface_parameters"],
    "surface_map": ["surface", "connection_nodes"],
    "surface_parameters": [],
    "weir": ["cross_section_definition", "connection_nodes"],
    "wind": [],
    "windshielding": ["channel"],
}


connected_tables_reversed = {}
for key1, values in connected_tables.items():
    connected_tables_reversed[key1] = []

    for key2, values2 in connected_tables.items():
        if key2 == key1:
            continue
        else:
            if key1 in values2:
                connected_tables_reversed[key1].append(key2)


class Properties(object):
    def __getitem__(self, key):
        return getattr(self, key)

    def __iter__(self):
        for i in dir(self):
            if "__" not in i:
                yield i, getattr(self, i)

    @property
    def schema_version(self):
        return {
            "version_num": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Version of your schematisation",
                "unit": "id",
            },
        }

    @property
    def v2_1d_boundary_conditions(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
            },
            "boundary_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2, 3, 5],
                "description": """
                                Boundary type.
                                1: waterlevel boundary
                                2: velocity boundary
                                3: discharge boundary 
                                5: Sommerfeld boundary (waterlevel slope). 

                                For types 2, 3 and 5 the channel direction determines sign of the input value. 
                                If the boundary is placed on the channel endpoint, 
                                positive values mean for example water is being extracted from the model.""",
                "unit": "1 M above datum, 2: m/s, 3: m3/s",
            },
            "connection_node_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique connection node id",
                "unit": None,
            },
            "timeseries": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Format has to be 'min,value min,value'",
                "unit": "min,value min,value",
            },
        }

    @property
    def v2_1d_lateral(self):
        return {
            "connected_tables": {"v2_connection_nodes": ["connection_node_id"]},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
            },
            "connection_node_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique connection node id",
                "unit": None,
            },
            "timeseries": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Lateral for 1D connection nodes",
                "unit": "min, m3/s",
            },
        }

    @property
    def v2_2d_boundary_conditions(self):
        return {
            "connected_tables": {"v2_connection_nodes": ["connection_node_id"]},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "boundary_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1],
                "description": "Boundary condition for 2d model edge (must be on edge of DEM file)",
                "unit": "1:m",
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": "1:m",
                "default": None,
            },
            "timeseries": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Lateral for 1D connection nodes",
                "unit": "min,value min,value",
                "default": None,
            },
        }

    @property
    def v2_2d_lateral(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "timeseries": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Discharge for lateral",
                "unit": "min, m3/s",
                "default": None,
            },
            "type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1],
                "description": "lateral on 2d",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_aggregation_settings(self):
        return {
            "connected_tables": {"v2_global_settings": ["global_settings_id"]},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "aggregation_method": {
                "type": ogr.OFTString,
                "optional": False,
                "options": [
                    "avg",
                    "min",
                    "max",
                    "cum",
                    "med",
                    "cum_negative",
                    "cum_positive",
                    "current",
                ],
                "description": "Method of aggregation, only use 'current' for volume and interception",
                "unit": None,
                "default": None,
            },
            "flow_variable": {
                "type": ogr.OFTString,
                "optional": True,
                "options": [
                    "discharge",
                    "flow_velocity",
                    "pump_discharge",
                    "rain",
                    "waterlevel",
                    "wet_cross-section",
                    "wet_surface",
                    "lateral_discharge",
                    "volume",
                    "simple_infiltration",
                    "leakage",
                    "interception",
                ],
                "description": "The name of output variable that is aggregated.",
                "unit": None,
                "default": None,
            },
            "global_settings_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "v2_global_settings scenario id",
                "unit": None,
                "default": None,
            },
            "timestep": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Timestep size for aggregation",
                "unit": "seconds",
                "default": None,
            },
            "var_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Method of aggregation, only use 'current' for volume and interception",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_channel(self):
        return {
            "connected_tables": {
                "v2_connection_nodes": [
                    "connection_node_end_id",
                    "connection_node_start_id",
                ]
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "calculation_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [100, 101, 102, 105],
                "description": """100 = embedded channel
                101 = stand-alone channel
                102 = connected channel
                105 = double connected channel
                Embedded or connected can only be used where a DEM is present. Any start-, end- or calculation node along a channel with these types may not lay outside the DEM.""",
                "unit": None,
                "default": False,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": None,
            },
            "connection_node_end_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "End node of channel line, must be present in v2_channel and the channel geometry endpoint must be snappen on the given connection node",
                "unit": None,
                "default": None,
            },
            "connection_node_start_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "End node of channel line, must be present in v2_channel and the channel geometry endpoint must be snappen on the given connection node",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name fields, no constraints",
                "unit": None,
                "default": None,
            },
            "dist_calc_points": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Distance between calculation points on line segments",
                "unit": "Meter",
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Visibility in livesite",
                "unit": None,
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbLineString,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_connection_nodes(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": False,
            },
            "initial_waterlevel": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Initial water level at connection node",
                "unit": "Meter above datum",
                "default": False,
            },
            "storage_area": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Storage area, e.g., for manholes in sewerage calculations",
                "unit": "Meter squared",
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbPoint,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_cross_section_definition(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": False,
            },
            "height": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "For tabulated fill in space-separated heights of profile",
                "unit": "Meter",
                "default": None,
            },
            "shape": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2, 3, 5, 6],
                "description": """
                1 = rectangle; specify width and height (profile/upper side is not automatically closed)
                2 = circle; specify width (profile/upper side is automatically closed)
                3 = egg; specify only 1 width. From this 3Di creates an egg-shaped profile with heigth = 1.5*width
                5 = tabulated rectangle; specify space-seperated width and height intervals. Between intervals the profile is defined straight. Can be closed by stating width 0 at heighest heigth
                6 = tabulated trapezium; specify space-sperated width and height intervals. Between intervals the profile is interpolated. Can be closed by stating width 0 at heighest heigth""",
                "unit": None,
                "default": None,
            },
            "width": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "For tabulated fill in space-separated widths of profile. Fill in diameter for circle.",
                "unit": "Meter",
                "default": None,
            },
        }

    @property
    def v2_cross_section_location(self):
        return {
            "connected_tables": {
                "v2_channel": ["channel_id"],
                "v2_cross_section_definition": ["definition_id"],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "bank_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "For connected channels only. Reference lavel for exchange between 1D and 2D.",
                "unit": "Meter above datum",
                "default": False,
            },
            "channel_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Reference to channel id.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": None,
            },
            "definition_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Reference to v2_cross_section_definition_id",
                "unit": None,
                "default": None,
            },
            "friction_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2],
                "description": """1 = Chezy !not yet implemented
                                  2 = Manning""",
                "unit": "1: m1/2/s 2: s/m1/3",
                "default": None,
            },
            "friction_value": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Friction or rougness value of profile.",
                "unit": "1: m1/2/s 2: s/m1/3",
                "default": None,
            },
            "reference_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Reference level or bottom level of profile",
                "unit": "Meter above datum",
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbPoint,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_culvert(self):
        return {
            "connected_tables": {
                "v2_cross_section_definition": ["cross_section_definition_id"],
                "v2_connection_nodes": [
                    "connection_node_start_id",
                    "connection_node_end_id",
                ],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "calculation_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [100, 101, 102, 105],
                "description": """
                100 = embedded channel
                101 = stand-alone channel
                102 = connected channel
                105 = double connected channel
                Embedded or connected can only be used where a DEM is present. Any start-, end- or calculation node along a channel with these types may not lay outside the DEM.""",
                "unit": "Meter above datum",
                "default": False,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": None,
            },
            "connection_node_end_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "End node of culvert, must be present in v2_connection_nodes.",
                "unit": None,
                "default": None,
            },
            "connection_node_start_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "End node of culvert, must be present in v2_connection_nodes.",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name fields, no constraints",
                "unit": None,
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Visibility in livesite",
                "unit": None,
                "default": None,
            },
            "cross_section_definition_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Reference to v2_cross_section_definition_id",
                "unit": None,
                "default": None,
            },
            "discharge_coefficient_negative": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Discharge coefficient for negative flow (from end node to start node). Can be set to 0 when closed",
                "unit": None,
                "default": None,
            },
            "discharge_coefficient_positive": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Discharge coefficient for positive flow (from end node to start node). Can be set to 0 when closed",
                "unit": None,
                "default": None,
            },
            "dist_calc_points": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Distance between calcultion points of linesegments",
                "unit": "Meter",
                "default": None,
            },
            "friction_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2],
                "description": "Friction type 1 = Chezy (not yet implemented) 2 = Manning",
                "unit": None,
                "default": None,
            },
            "friction_value": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Friction or roughness value for profile",
                "unit": "1: m1/2/s 2: s/m1/3",
                "default": None,
            },
            "invert_level_end_point": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Invert level at culvert endpoint",
                "unit": "Meter above datum",
                "default": None,
            },
            "invert_level_start_point": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Invert level at culvert endpoint",
                "unit": "Meter above datum",
                "default": None,
            },
        }

    @property
    def v2_global_settings(self):
        return {
            "connected_tables": {"v2_numerical_settings": ["numerical_settings_id"]},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "advection_1d": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1],
                "description": "0: off, 1: on, Use advection in 1D, other schemes 2-7 are in experimental fase.",
                "unit": None,
                "default": 0,
            },
            "advection_2d": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1],
                "description": "0: off, 1: on, Use advection in 2D.",
                "unit": None,
                "default": 0,
            },
            "control_group_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Reference to id that contains control settings for this scenario.",
                "unit": None,
                "default": None,
            },
            "dem_file": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Relative to dem file (.tif)",
                "unit": "Meter above datum",
                "default": None,
            },
            "dem_obstacle_detection": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Automatically detect obstacles based on DEM-file.Works only in combination with dem_obstacle_height (has no relation with v2_obstacle).",
                "unit": None,
                "default": None,
            },
            "dem_obstacle_height": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Relative height (above lowest pixel of calc cell) for obstacle detection.",
                "unit": "Above lowest pixel of calc cell",
                "default": None,
            },
            "dist_calc_points": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Global distance between calculation points for line elements",
                "unit": "Meter",
                "default": None,
            },
            "embedded_cutoff_threshold": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Relative length of cell size. When embedded channel intersects 2D grid size over length shorter than the cellsize * cutoff threshold, the embedded channel skips this 2D cell. Usefull for preventing very short embedded channel segments (which slow down your model).",
                "unit": "Factor between 0 and 1",
                "default": 0.05,
            },
            "epsg_code": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Define map projection for study area. Must match raster projection.",
                "unit": "EPSG projection",
                "default": 28992,
            },
            "flooding_threshold": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Water depth threshold for flow between 2D cells. Depth relative to lowest bathymetry pixel at the edge between two 2D cell.",
                "unit": "Meter",
                "default": 0.0,
            },
            "frict_avg": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1],
                "description": "The roughness coefficient will be average within one cell.",
                "unit": None,
                "default": 0,
            },
            "frict_coef": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Constant friction coefficient for 2D.",
                "unit": "1: m1/2/s 2: s/m1/3",
                "default": None,
            },
            "frict_coef_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Relative path to friction file",
                "unit": "1: m1/2/s 2: s/m1/3",
                "default": None,
            },
            "frict_type": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [1, 2],
                "description": "1: Chezy for 2D, 2: Manning for 2D.",
                "unit": "1: m1/2/s 2: s/m1/3",
                "default": None,
            },
            "grid_space": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Smallest grid size in quadtree, k=1, must be at least twice the cellsize of the raster",
                "unit": "Meters",
                "default": None,
            },
            "groundwater_settings_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Reference to id that contains groundwater settings for this scenario.",
                "unit": None,
                "default": None,
            },
            "initial_groundwater_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Initial groundwater level",
                "unit": "Meter above datum",
                "default": None,
            },
            "initial_groundwater_level_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Relative path to groundwater file",
                "unit": "Meter above datum",
                "default": None,
            },
            "initial_groundwater_level_type": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2],
                "description": """0=max, 1=min, 2=avg
                                  **Mandatory when using initial water level file.
                                  """,
                "unit": None,
                "default": None,
            },
            "initial_waterlevel": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Global initial waterlevel",
                "unit": "Meter above datum",
                "default": None,
            },
            "initial_waterlevel_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Relative path to initial waterlevel file file",
                "unit": "Meter above datum",
                "default": None,
            },
            "interflow_settings_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Reference to interflow settings id",
                "unit": "Meter above datum",
                "default": None,
            },
            "interception_global": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Global value for interception.",
                "unit": "Meter",
                "default": None,
            },
            "interception_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Relative path to interception file.",
                "unit": "Meter",
                "default": None,
            },
            "kmax": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Maximum multitude of smallest grid size in quadtree starting from grid_space at k=1. Grid size increases according to 2^(k-1) * grid_space.",
                "unit": "Meter",
                "default": None,
            },
            "manhole_storage_area": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """ Default manhole storage area. This is the surface area that each manhole is given when water reaches above the surface level.
                    **Mandatory when using only 1d flow (no dem) manhole area must be larger than 0 (and an INTEGER)
                    Must be NULL when using only 2d.""",
                "unit": "Meter squared",
                "default": None,
            },
            "max_angle_1d_advection": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Maximum angle at which advection is considered.",
                "unit": "degrees [0-90]",
                "default": None,
            },
            "maximum_sim_time_step": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """Maximum timestep during simulation.
                                 ** Mandatory when using timestep plus.""",
                "unit": "Seconds",
                "default": None,
            },
            "minimum_sim_time_step": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Minimum timestep during simulation, keep it small for stability",
                "unit": "Seconds",
                "default": 0.01,
            },
            "maximum_table_step_size": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Maximum step size for tables.",
                "unit": "meters",
                "default": None,
            },
            "name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name",
                "unit": None,
                "default": None,
            },
            "nr_timesteps": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Maximum nr of timesteps",
                "unit": None,
                "default": None,
            },
            "numerical_settings_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Reference to id that contains numerical settings for this scenario.",
                "unit": None,
                "default": None,
            },
            "output_time_step": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Timestep written in the output file.",
                "unit": "Seconds",
                "default": None,
            },
            "sim_time_step": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Simulation time step.",
                "unit": "Seconds",
                "default": None,
            },
            "simple_infiltration_settings_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Reference to id that contains settings for simple infiltration for this scenario.",
                "unit": "Seconds",
                "default": None,
            },
            "start_date": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Start datum",
                "unit": "Format: 2017-01-01",
                "default": None,
            },
            "start_time": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Starttime of simulation.",
                "unit": "Format: 00:00:00",
                "default": None,
            },
            "table_step_size": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "User-defined table stepsize/increment (m).",
                "unit": "Meters",
                "default": None,
            },
            "table_step_size_1d": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "User-defined table stepsize/increment (m) for 1d cross-sections and volumes.  default value = table_step_size.",
                "unit": "Meters",
                "default": None,
            },
            "table_step_size_volume_2d": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "User-defined table stepsize/increment (m) for defining 2D volumes. Can increase speed when this is set larger than  table_step_size. default value = table_step_size.",
                "unit": "Meters",
                "default": None,
            },
            "timestep_plus": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1],
                "description": "Allow switching to larger timestep when simulation is steady.",
                "unit": None,
                "default": 0,
            },
            "use_0d_inflow": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Include 0D inflow (NWRW/impervious surfaces) in simulation.",
                "unit": None,
                "default": None,
            },
            "use_1d_flow": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Include 1D flow in simulation.",
                "unit": None,
                "default": None,
            },
            "use_2d_flow": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Include 2D flow in simulation.",
                "unit": None,
                "default": None,
            },
            "use_2d_rain": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Use rainfall via 2D surface for this scenario.",
                "unit": None,
                "default": None,
            },
            "water_level_ini_type": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Use rainfall via 2D surface for this scenario.",
                "unit": None,
                "default": None,
            },
            "wind_shielding_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "NOT YET IMPLEMENTED",
                "unit": None,
                "default": None,
            },
            "guess_dams": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "NOT YET IMPLEMENTED",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_groundwater(self):
        return {}  # TODO

    @property
    def v2_grid_refinement(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": False,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "refinement_level": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Local refinement level. Starting from 1. Values above kmax (v2_global_settings) are ignored.",
                "unit": None,
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbLineString,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_grid_refinement_area(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": False,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "refinement_level": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Local refinement level. Starting from 1. Values above kmax (v2_global_settings) are ignored.",
                "unit": None,
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbPolygon,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_impervious_surface(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "area": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "area of surface",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "dry_weather_flow": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Dry weather flow per inhabitant",
                "unit": "L/day per inhabitant",
                "default": None,
            },
            "nr_of_inhabitants": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Number of inhabitant used for dry wheather flow.",
                "unit": "Nr of inhabitants",
                "default": None,
            },
            "surface_class": {
                "type": ogr.OFTString,
                "optional": False,
                "options": [
                    "impervious paving",
                    "pervious paving",
                    "semi-pervious paving",
                    "unpaved",
                    "building",
                ],
                "description": """gesloten verharding -- impervious paving
                                  open verharding       -- pervious paving
                                  half verhard               -- semi-pervious paving
                                  onverhard                  -- unpaved
                                  pand                         -- building""",
                "unit": None,
                "default": None,
            },
            "surface_sub_class": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": """unkown""",
                "unit": None,
                "default": None,
            },
            "surface_inclination": {
                "type": ogr.OFTString,
                "optional": False,
                "options": ["level", "inclined", "elongated"],
                "description": """vlak              -- level
                                hellend        -- inclined
                                uitgestrekt  -- elongated """,
                "unit": None,
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_impervious_surface_map(self):
        return {
            "connected_tables": {
                "v2_connection_nodes": ["connection_node_id"],
                "v2_impervious_surface": ["impervious_surface_id"],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "connection_node_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of connection node",
                "unit": None,
                "default": None,
            },
            "impervious_surface_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of impervious surface feature",
                "unit": None,
                "default": None,
            },
            "percentage": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Percentage of imperivous surface are places on connection nodes/",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_interflow(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "hydraulic_conductivity": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """
                            Global hydraulic conductivity (Darcy)
                            **When interflow_type > 0 then hydraulic_conductivity OR hydraulic_conductivity_file is mandatory""",
                "unit": "m/day",
                "default": None,
            },
            "hydraulic_conductivity_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": """
                            Global hydraulic conductivity (Darcy)
                            **When interflow_type > 0 then hydraulic_conductivity OR hydraulic_conductivity_file is mandatory""",
                "unit": "m/day",
                "default": None,
            },
            "impervious_layer_elevation": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """
                                When using interflow: Depth of interflow layer defined below lowest pixel (so always positive). Imaginary bottom of interflow layer. 
                                For interflow types 1 and 2 it is ignored for the volume in the interflow layer (but still it must be filled in when using interflow). The volume in these types is determined by the porosity and the porosity layer thickness. 
                                For interflow types 3 and 4 it is used to determine the volume in the interflow layer. 
                                In all types the waterlevel in the interflow layer starts at this level. It does not influence flow.
                """,
                "unit": "Meter",
                "default": None,
            },
            "interflow_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": """
                                Include interflow in simulation.
                                0: No Interflow
                                1: define 1 porosity value for model. This porsity will be rescaled per pixel (to lowest pixel per cell), so (interflow) volume is the same for each pixel within 1 cell --> define porosity, hydraulic_conductivity, porosity_layer_thickness and impervious_layer_elevation (can be used in combination with groundwater)
                                2: define 1 porosity value for model. This porsity will be rescaled per pixel (to lowest pixel whole model), so (interflow) volume is the same for each pixel in whole model --> define porosity, hydraulic_conductivity, porosity_layer_thickness and impervious_layer_elevation (cannot be used in combination with groundwater)
                                3: define 1 porosity value for model. This porsity will not be rescaled, but each pixel in the model has the same porosity. The (interflow) volume for each pixel depends on the impervious_layer_elevation, which is below lowest pixel of cell --> define porosity, hydraulic_conductivity and impervious_layer_elevation (can be used in combination with groundwater)
                                4: define 1 porosity value for model. This porsity will not be rescaled, but each pixel in the model has the same porosity. The (interflow) volume for each pixel depends on the impervious_layer_elevation, which is below lowest pixel of whole model --> define porosity, hydraulic_conductivity and impervious_layer_elevation (cannot be used in combination with groundwater)
                                
                                * NOT NULL
                                ** do not use interflow in combination with limiter_slope_crossectional_area_2d >0 AND/OR limiter_slope_friction_2d>0
                """,
                "unit": None,
                "default": 0,
            },
            "porosity": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """
                                Porosity (between 0 and 1) of interflow layer.
                                ** Mandatory when using interflow
                """,
                "unit": None,
                "default": None,
            },
            "porosity_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": " Relative path to porosity file",
                "unit": None,
                "default": None,
            },
            "porosity_layer_thickness": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """Thickness of porosity layer relative to DEM.
                                 **Mandatory for interflow_type 1 and 2.""",
                "unit": "Meter",
                "default": None,
            },
        }

    @property
    def v2_levee(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "crest_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Crest level of levee segment",
                "unit": "Meter above datum",
                "default": None,
            },
            "material": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [1, 2],
                "description": """** Mandatory when you want to use a levee breach during your calculation 
                                     Material used for breach growth.
                                     1: sand
                                     2: clay""",
                "unit": None,
                "default": None,
            },
            "max_breach_depth": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """** Mandatory when you want to use a levee breach during your calculation 
                                    Maximum breach depth relative to crest level (thus a positive value must be filled in).""",
                "unit": "Meter below crest level",
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbLineString,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_manhole(self):
        return {
            "connected_tables": {"v2_connection_nodes": ["connection_node_id"]},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "bottom_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Manhole bottom level.",
                "unit": "Meter above datum",
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "calculation_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1, 2],
                "description": """ 
                                Manhole calculation type for 1D-2D connection.
                                0: embedded
                                1: isolated
                                2: connected""",
                "unit": "L/day per inhabitant",
                "default": None,
            },
            "connection_node_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of connection ndoe on which manhole is placed",
                "unit": None,
                "default": None,
            },
            "drain_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Manhole drain level (**for connected manholes). If there is a connected manhole without drain level, 3Di will take the top of the pipe from the connection pipes as drain level.",
                "unit": "Meter above datum",
                "default": None,
            },
            "length": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Manhole length. This value is for administrative purposes only and has no effect on the storage area of the connection node. **Mandatory when shape = 02",
                "unit": "Meter",
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "manhole_indicator": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1, 2],
                "description": """  0: inspection (inspectieput) 1: outlet 2: pump""",
                "unit": None,
                "default": None,
            },
            "shape": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": """Manhole shape. This value is for administrative purposes only and has no effect on the storage area of the connection node. To add storage to a connection node, adjust the 'storage_area' in the v2_connection_nodes table.
                                00: square
                                01: round
                                02: rectangle """,
                "unit": None,
                "default": None,
            },
            "surface_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Manhole surface level.",
                "unit": "Meter above datum",
                "default": None,
            },
            "sediment_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Manhole surface level.",
                "unit": "Meter above datum",
                "default": None,
            },
            "width": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Manhole width or diameter. This value is for administrative purposes only and has no effect on the storage area of the connection node.",
                "unit": "Meter",
                "default": None,
            },
        }

    @property
    def v2_numerical_settings(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "cfl_strictness_factor_1d": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Strictness of CFL for 1D",
                "unit": None,
                "default": 1,
            },
            "cfl_strictness_factor_2d": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Strictness of CFL for 2D",
                "unit": None,
                "default": 1,
            },
            "convergence_cg": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """For numerical computation several tresholds are needed in the code, 
                    to avoid deficiencies due to a limited numerical accuracy. Generally this is to keep the behaviour consistent.""",
                "unit": None,
                "default": 0.000000001,
            },
            "convergence_eps": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """Minimal residual for convergence of newton iteration.""",
                "unit": None,
                "default": 0.0001,
            },
            "flow_direction_threshold": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """For numerical computation several tresholds are needed in the code, 
                    to avoid deficiencies due to a limited numerical accuracy. Generally this is to keep the behaviour consistent.""",
                "unit": "Meter per second",
                "default": 0.000001,
            },
            "frict_shallow_water_correction": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": """In case the friction assumptions based on the dominant friction balance gives a structurally underestimation of the friction, one can switch this setting on. 
                                0 is off, 1 is maximum between averaged friction and divided channel based friction, 2 is always linearized, 3 linearizes the depth based on a weighed averaged.  In this case the maximum depth of a thin layer needs to be defined. Do not use in combination with interflow.""",
                "unit": None,
                "default": 0,
            },
            "general_numerical_threshold": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """For numerical computation several tresholds are needed in the code, 
                    to avoid deficiencies due to a limited numerical accuracy. Generally this is to keep the behaviour consistent.""",
                "unit": None,
                "default": 0.00000001,
            },
            "integration_method": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": """Time integration method: 0=Euler implicit.""",
                "unit": 0,
                "default": 0.000000001,
            },
            "limiter_grad_1d": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": """The limiter on the water level gradient allows the model to deal with unrealistically steep gradients.""",
                "unit": None,
                "default": 1,
            },
            "limiter_grad_2d": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": """The limiter on the water level gradient allows the model to deal with unrealistically steep gradients. When field is left empty, it is switched on!!!""",
                "unit": None,
                "default": 0,
            },
            "limiter_slope_crossectional_area_2d": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": """This limiter starts working in case the depth based on the downstream water level is zero and may be useful in sloping areas.
                                  0 is off, and 1 is a limiter which ends in a a higher order scheme, but is sensitive too instabilities, 2, treats the cross-sections as an upwind method volume/surface area under the assumption that the flow acts like a thin layer,  3 makes a combination of the traditional method in combination with the  thin layer apporach. In this case the maximum depth of a thin layer needs to be defined.
                                  Do not use in combination with interflow""",
                "unit": None,
                "default": 0,
            },
            "limiter_slope_friction_2d": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": """This limiter starts working in case the depth based on the downstream water level is zero and may be useful in sloping areas.  0 is off, and 1 is on. 
                            This limiter is obliged in combination with limiter_slope_crossectional_area_2d>0.. l Do not use in combination with interflow""",
                "unit": None,
                "default": 0,
            },
            "max_degree": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [700, 7, 5, 7, 70],
                "description": """Setting for matrix solver. Values below are adviced for different model types
                        700 for 1D flow
                        7 for 1D and 2D flow
                        5 for surface 2D flow only
                        7 for surface and groundwater flow
                        70 for 1D, 2D surface and groundwater flow or higher. Play around with this value in case of groundwater, can speed up your model significantly""",
                "unit": None,
                "default": 0,
            },
            "max_nonlin_iterations": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": """Maximum number of nonlinear iterations in single time step.""",
                "unit": None,
                "default": 20,
            },
            "minimum_friction_velocity": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """For numerical computation several tresholds are needed in the code, to avoid deficiencies due to a limited numerical accuracy. Generally this is to keep the behaviour consistent..""",
                "unit": "m/s",
                "default": 0.05,
            },
            "minimum_surface_area": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """Maximum number of nonlinear iterations in single time step.""",
                "unit": "m2",
                "default": 20,
            },
            "precon_cg": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1],
                "description": """Use preconditioner for matrix solver. Increases simulation speed in most cases, Set to 0 or 1 (default).""",
                "unit": None,
                "default": 1,
            },
            "preissmann_slot": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """A conceptual vertical and narrow slot providing a conceptual free surface condition for the flow when the water level is above the top of a closed conduit. Often used to guarentee stability, in 3Di unnessary unless used for presurized pipe flow. Works only for circular profiles.""",
                "unit": "m2",
                "default": 0,
            },
            "pump_implicit_ratio": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": [0, 1],
                "description": """Determines whether pump discharge is always maximum capacity (0) or discharge is limited to available inflow (1). The latter ensures a smooth discharge. Value between 0 and 1.""",
                "unit": None,
                "default": 1,
            },
            "thin_water_layer_definition": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": """** mandatory when using friction shallow water correction option 3 or limiter_slope_crossectional_area_2d on option 3""",
                "unit": "m",
                "default": 0.05,
            },
            "use_of_cg": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": """Number of iteration of conjugate gradient method, before swtiching to another method.""",
                "unit": "m",
                "default": 20,
            },
            "use_of_nested_newton": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": """1 for 1D calculation with closed profiles to handle non-linearity in volume-waterlevel relation. When using 0 nested newton is switch off by default but will be used when calculations become non-linear. For sewerage systems 1 is advised.""",
                "unit": "m",
                "default": 0,
            },
        }

    @property
    def v2_obstacle(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints",
                "unit": None,
                "default": False,
            },
            "crest_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Crest level of obstacle segment.",
                "unit": "Meter above datum",
                "default": None,
            },
        }

    @property
    def v2_orifice(self):
        return {
            "connected_tables": {
                "v2_cross_section_definition": ["cross_section_definition_id"],
                "v2_connection_nodes": [
                    "connection_node_start_id",
                    "connection_node_end_id",
                ],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "connection_node_start_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Start node for orifice",
                "unit": None,
                "default": None,
            },
            "connection_node_end_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "End node for orifice",
                "unit": None,
                "default": None,
            },
            "crest_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Crest or bottom level",
                "unit": "Meter above datum",
                "default": None,
            },
            "crest_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [3, 4],
                "description": """Type of weir formulation.
                                  3: broad crested
                                  4: short crested""",
                "unit": None,
                "default": 4,
            },
            "cross_section_definition_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of cross section definiton in v2_cross_section_definition",
                "unit": None,
                "default": None,
            },
            "discharge_coefficient_negative": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Discharge coefficient for negative flow (from end to start node). Can be set to 0 when closed.",
                "unit": None,
                "default": 1,
            },
            "discharge_coefficient_positive": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Discharge coefficient for postive flow (from start to end node). Can be set to 0 when closed.",
                "unit": None,
                "default": 1,
            },
            "friction_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2],
                "description": """Friction Type.
                                1: Chezy (not yet implemented)
                                 2: Manning""",
                "unit": None,
                "default": 2,
            },
            "friction_value": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Friction or roughness value for profile",
                "unit": "1:[m1/2/s], 2:[s/m1/3]",
                "default": None,
            },
            "sewerage": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "For internal book keeping. Can be used for statistics in QGIS plugin.",
                "unit": None,
                "default": None,
            },
            # TODO
            # Has been removed after klondike
            # "max_capacity": {
            #     "type": ogr.OFTReal,
            #     "optional": True,
            #     "options": None,
            #     "description": "Not Needed",
            #     "unit": None,
            #     "default": None,
            # },
        }

    @property
    def v2_pipe(self):
        return {
            "connected_tables": {
                "v2_cross_section_definition": ["cross_section_definition_id"],
                "v2_connection_nodes": [
                    "connection_node_end_id",
                    "connection_node_start_id",
                ],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "connection_node_start_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Start node for pipe",
                "unit": None,
                "default": None,
            },
            "connection_node_end_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "End node for pipe",
                "unit": None,
                "default": None,
            },
            "calculation_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1, 2, 3, 4],
                "description": """Calculation type for pipe. When start en end connection nodes are manholes only used for calculation points half-way pipe.
                                    0 = embedded
                                    1 = isolated
                                    2 = connected
                                    3 = broad crest
                                    4 = short crest""",
                "unit": None,
                "default": 1,
            },
            "dist_calc_points": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": [3, 4],
                "description": """Type of weir formulation.
                                  3: broad crested
                                  4: short crested""",
                "unit": None,
                "default": 4,
            },
            "cross_section_definition_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of cross section definiton in v2_cross_section_definition",
                "unit": None,
                "default": None,
            },
            "invert_level_start_point": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Invert level at pipe startpoint. ",
                "unit": None,
                "default": None,
            },
            "invert_level_end_point": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Invert level at pipe endpoint. ",
                "unit": None,
                "default": None,
            },
            "friction_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2],
                "description": """Friction Type.
                                1: Chezy (not yet implemented)
                                 2: Manning""",
                "unit": None,
                "default": 2,
            },
            "friction_value": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Friction or roughness value for profile",
                "unit": "1:[m1/2/s], 2:[s/m1/3]",
                "default": None,
            },
            "material": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5, 6, 7, 8],
                "description": """Material of pipe, used for internal bookkeeping only.
                                  0: concrete
                                 1: pvc
                                 2: gres
                                 3: cast iron
                                 4: brickwork
                                 5: HPE
                                 6: HDPE
                                 7: plate iron
                                 8: steel""",
                "unit": None,
                "default": None,
            },
            "original_length": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "For internal book keeping.",
                "unit": "Meter",
                "default": None,
            },
            "profile_num": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "For internal book keeping.",
                "unit": None,
                "default": None,
            },
            "sewerage_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1, 2, 3, 4, 5, 6, 7],
                "description": """ 
                            Pipe type. 3Di requires the sewerage_type to be one of the following:
                            0: gemengd - mixed
                             1: rwa - rain water
                             2: dwa - dry wheather flow
                             3: transport
                             4: overstort - spillway
                             5: zinker
                             6: berging - storage
                             7: bergbezinkbak - storage tank
                             """,
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_pumpstation(self):
        return {
            "connected_tables": {
                "v2_connection_nodes": [
                    "connection_node_start_id",
                    # "connection_node_end_id",
                ]
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "capacity": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Pump capacity.",
                "unit": "Liters per second",
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "connection_node_start_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Start node for pumpstation",
                "unit": None,
                "default": None,
            },
            "connection_node_end_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "End node for pumpstation",
                "unit": None,
                "default": None,
            },
            "classification": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "TODO",
                "unit": None,
                "default": None,
            },
            "lower_stop_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Level at pump start or end node at which pump stops pumping.",
                "unit": "Meter above datum",
                "default": None,
            },
            "sewerage": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "For internal book keeping. Specifies if the pumpstation is in sewerage",
                "unit": None,
                "default": None,
            },
            "start_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Level at pump start or end node from from which it starts pumping.",
                "unit": "Meter above datum",
                "default": None,
            },
            "type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2],
                "description": """
                            Type that determines pump function.
                            1 Pump behaviour is based on water levels on the suction-side of the pump
                            2 Pump behaviour is based on water levels on the delivery-side of the pump  """,
                "unit": None,
                "default": None,
            },
            "upper_stop_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Level at pump start or end node at which pump stops pumping.Must be above start level.",
                "unit": "Meter above datum",
                "default": None,
            },
        }

    @property
    def v2_simple_infiltration(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "infiltration_rate": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "global infiltration raste",
                "unit": "mm/day",
                "default": 0,
            },
            "infiltration_rate_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Relative path to infiltration file (.tif).",
                "unit": "mm/day",
                "default": None,
            },
            "infiltration_surface_option": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2],
                "description": """Option that sets how the infiltration works in calculation cells.
                                0: rain (whole surface when raining, only wet pixels when dry)
                                1: whole surface (always whole surface)
                                2: only wet surface (always only wet pixels)
                                in case not defined then option 0 is used
                                """,
                "unit": None,
                "default": 0,
            },
            "max_infiltration_capacity_file": {
                "type": ogr.OFTString,
                "optional": True,
                "options": [1, 2],
                "description": """Relative path to max infiltration file (.tif). Maximum infiltration uses the sum of pixel values per calculation cell.""",
                "unit": "Meter",
                "default": None,
            },
        }

    @property
    def v2_surface(self):
        return {
            "connected_tables": {"v2_surface_parameters": ["surface_parameters_id"]},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "area": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "area of surface",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "dry_weather_flow": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Dry weather flow per inhabitant",
                "unit": "L/day per inhabitant",
                "default": None,
            },
            "nr_of_inhabitants": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Number of inhabitant used for dry wheather flow.",
                "unit": "Nr of inhabitants",
                "default": None,
            },
            "function": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": """For your own administration.""",
                "unit": None,
                "default": None,
            },
            "surface_parameters_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": """Reference to v2_surface_parameters. The id filled in here must be present in this (v2_surface_parameters) table """,
                "unit": None,
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_surface_map(self):
        return {
            "connected_tables": {
                "v2_surface": ["surface_id"],
                "v2_connection_nodes": ["connection_node_id"],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "connection_node_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of connection_node. Connecting the surface area from v2_surface to an connection_node",
                "unit": None,
                "default": None,
            },
            "percentage": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "percentage of area to the connection_node",
                "unit": "%",
                "default": None,
            },
            "surface_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of surface feature",
                "unit": None,
                "default": None,
            },
            "surface_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": ["v2_surface", "v2_impervious_surface"],
                "description": "choice to use the 'v2_surface' or 'v2_impervious_surface'",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_surface_parameters(self):
        return {
            "connected_tables": {},
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "infiltration": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1],
                "description": "0 or 1",
                "unit": None,
                "default": None,
            },
            "infiltration_decay_constant": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "time factor decay infiltration capacity of the surface",
                "unit": "/hour",
                "default": None,
            },
            "infiltration_recovery_constant": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "time factor recovery infiltration capacity of the surface",
                "unit": "/hour",
                "default": None,
            },
            "max_infiltration_capacity": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Maximum infiltration capacity of the surface",
                "unit": "mm/h",
                "default": None,
            },
            "min_infiltration_capacity": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Minimum infiltration capacity of the surface",
                "unit": "mm/h",
                "default": None,
            },
            "outflow_delay": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "delay of outflow",
                "unit": "/min",
                "default": None,
            },
            "surface_layer_thickness": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "mm storage on the surface",
                "unit": "mm",
                "default": None,
            },
        }

    @property
    def v2_weir(self):
        return {
            "connected_tables": {
                "v2_cross_section_definition": ["cross_section_definition_id"],
                "v2_connection_nodes": [
                    "connection_node_start_id",
                    "connection_node_end_id",
                ],
            },
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "display_name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "code": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "zoom_category": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": [0, 1, 2, 3, 4, 5],
                "description": "Name field, no constraints.",
                "unit": None,
                "default": None,
            },
            "connection_node_start_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Start node for weir",
                "unit": None,
                "default": None,
            },
            "connection_node_end_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "End node for weir",
                "unit": None,
                "default": None,
            },
            "crest_level": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Crest or bottom level",
                "unit": "Meter above datum",
                "default": None,
            },
            "crest_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [3, 4],
                "description": """Type of weir formulation.
                                  3: broad crested
                                  4: short crested""",
                "unit": None,
                "default": 4,
            },
            "cross_section_definition_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of cross section definiton in v2_cross_section_definition",
                "unit": None,
                "default": None,
            },
            "discharge_coefficient_negative": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Discharge coefficient for negative flow (from end to start node). Can be set to 0 when closed.",
                "unit": None,
                "default": 1,
            },
            "discharge_coefficient_positive": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Discharge coefficient for postive flow (from start to end node). Can be set to 0 when closed.",
                "unit": None,
                "default": 1,
            },
            "friction_type": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [1, 2],
                "description": """Friction Type.
                                1: Chezy (not yet implemented)
                                 2: Manning""",
                "unit": None,
                "default": 2,
            },
            "friction_value": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Friction or roughness value for profile",
                "unit": "1:[m1/2/s], 2:[s/m1/3]",
                "default": None,
            },
            "sewerage": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "For internal book keeping. Can be used for statistics in QGIS plugin.",
                "unit": None,
                "default": None,
            },
            "external": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "For internal book keeping. Can be used for statistics in QGIS plugin.",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_calculation_point(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "content_type_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Content type",
                "unit": None,
                "default": None,
            },
            "user_ref": {
                "type": ogr.OFTString,
                "optional": True,
                "options": None,
                "description": "Reference",
                "unit": None,
                "default": None,
            },
            "calc_type": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Calculation type",
                "unit": None,
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbPoint,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_connected_pnt(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "exchange_level": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "exchange level for leveees",
                "unit": None,
                "default": None,
            },
            "calculation_pnt_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "Caculation point id",
                "unit": None,
                "default": None,
            },
            "levee_id": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "id of levee",
                "unit": None,
                "default": None,
            },
            "the_geom": {
                "type": ogr.wkbPoint,
                "optional": False,
                "options": None,
                "description": "Geometry",
                "unit": "Depend on projection",
                "default": None,
            },
        }

    @property
    def v2_control(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "control_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": ["table", "memory", "timed"],
                "description": "Type of control, options are amongst others: table, memory",
                "unit": None,
                "default": None,
            },
            "control_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "id in the v2_control_table (in case of table control) or v2_control_memory (in case of memory control)",
                "unit": None,
                "default": None,
            },
            "control_group_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "id of the v2_control_group this control is part of",
                "unit": None,
                "default": None,
            },
            "measure_group_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "id of the v2_measure_group",
                "unit": None,
                "default": None,
            },
            "start": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "Start time of the control in seconds since beginning of the simulation. Can be used to link multiple control tables to one control structure.",
                "unit": "seconds",
                "default": None,
            },
            "end": {
                "type": ogr.OFTReal,
                "optional": True,
                "options": None,
                "description": "End time of the control in seconds since beginning of the simulation",
                "unit": None,
                "default": None,
            },
            "measure_frequency": {
                "type": ogr.OFTInteger,
                "optional": True,
                "options": None,
                "description": "NOT YET IMPLEMENTED",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_control_group(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "name": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Name",
                "unit": None,
                "default": None,
            },
            "description": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Description",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_control_measure_group(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            }
        }

    @property
    def v2_control_measure_map(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "measure_group_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "ID of the v2_measure_group this measure station is part of.",
                "unit": None,
                "default": None,
            },
            "object_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Type of object to measure at, for example:'v2_connection_nodes'",
                "unit": None,
                "default": None,
            },
            "object_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "id of the object (of type defined in object_type)",
                "unit": None,
                "default": None,
            },
            "weight": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "weight of measuring station in group, use 1.0 for groups with single measuring station. Combined weight should be 1.0.",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_control_table(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "action_table": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": """Semicolumn seperated table with action values, use # for newline. For Example: -1.7;-1.4#-1.6;-1.3#-1.5;-1.2
                                  When controlling set_discharge_coefficients you need to supply 2 values. One for the positive discharge coefficient and one for the negative discharge coefficient. Example: -1.7;0 0#-1.6;0.5 0.8#-1.5;1 1
                               """,
                "unit": None,
                "default": None,
            },
            "action_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Type of action; For instance: 'set_crest_level' or 'set_discharge_coefficients'",
                "unit": None,
                "default": None,
            },
            "measure_variable": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Measure variable in action table. For instance: 'waterlevel'",
                "unit": None,
                "default": None,
            },
            "measure_operator": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Operator for direction the action table is read. '<' or '>'",
                "unit": None,
                "default": None,
            },
            "target_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Structure type the control is applied to. For instance: 'v2_weir' or 'v2_culvert' or 'v2_orifice'",
                "unit": None,
                "default": None,
            },
            "target_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Id of structure the control is applied to.",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_control_memory(self):
        return {
            "id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Unique identifier",
                "unit": None,
                "default": None,
            },
            "action_value": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": """Value that the measure_variable is set to when memory control becomes active""",
                "unit": None,
                "default": None,
            },
            "action_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Type of action; For instance: set_pump_capacity",
                "unit": None,
                "default": None,
            },
            "is_active": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1],
                "description": """0: control is inactive when initializing the model
                                  1: control is active when initializing the model""",
                "unit": None,
                "default": None,
            },
            "is_inverse": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": [0, 1],
                "description": """
                                0: normal functioning of the control
                                1: inverting the lower and upper threshold
                            """,
                "unit": None,
                "default": None,
            },
            "lower_threshold": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Lower threshold of measure_variable. Control becomes inactive when value drops below this value (unless is_inverse = 1)",
                "unit": None,
                "default": None,
            },
            "measure_variable": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Measure variable in action table. For instance: waterlevel",
                "unit": None,
                "default": None,
            },
            "target_type": {
                "type": ogr.OFTString,
                "optional": False,
                "options": None,
                "description": "Structure type the control is applied to. For instance: 'v2_weir' or 'v2_culvert' or 'v2_orifice'",
                "unit": None,
                "default": None,
            },
            "target_id": {
                "type": ogr.OFTInteger,
                "optional": False,
                "options": None,
                "description": "Id of structure the control is applied to.",
                "unit": None,
                "default": None,
            },
            "upper_threshold": {
                "type": ogr.OFTReal,
                "optional": False,
                "options": None,
                "description": "Upper threshold of measure_variable. Control becomes active when value rises above this value (unless is_inverse = 1)",
                "unit": None,
                "default": None,
            },
        }

    @property
    def v2_control_delta(self):
        return {}  # TODO

    @property
    def v2_control_pid(self):
        return {}  # TODO

    @property
    def v2_control_timed(self):
        return {}  # TODO

    @property
    def v2_windshielding(self):
        return {}  # TODO

    @property
    def v2_dem_average_area(self):
        return {}  # TODO

    @property
    def v2_floodfill(self):
        return {}  # TODO

    @property
    def v2_pumped_drainage_area(self):
        return {}  # TODO

    @property
    def v2_initial_waterlevel(self):
        return {}  # TODO

    @property
    def v2_outlet(self):
        return {}  # TODO


TRANSLATE = {
    "v2_global_settings": "global_settings",
    "v2_connection_nodes": "nodes",
    "v2_manhole": "manholes",
    "v2_pipe": "pipes",
    "v2_simple_infiltration": "simple_infiltration",
    "v2_numerical_settings": "numerical_settings",
    "v2_aggregation_settings": "aggregation_settings",
    "v2_obstacle": "obstacles",
    "v2_channel": "channels",
    "v2_cross_section_definition": "cross_section_definitions",
    "v2_cross_section_location": "cross_section_locations",
    "v2_grid_refinement_area": "grid_refinement_areas",
    "v2_grid_refinement": "grid_refinements",
    "v2_1d_boundary_conditions": "boundary_conditions_1d",
    "v2_1d_lateral": "laterals_1d",
    "v2_pumpstation": "pumpstations",
    "v2_levee": "levees",
    "v2_weir": "weirs",
    "v2_orifice": "orifices",
    "v2_culvert": "culverts",
    "v2_calculation_point": "calculation_points",
    "v2_connected_pnt": "connection_points",
    "v2_control": "control",
    "v2_control_group": "control_group",
    "v2_control_measure_group": "control_measure_group",
    "v2_control_measure_map": "control_measure_map",
    "v2_control_table": "control_table",
    "v2_control_memory": "control_memory",
    "v2_impervious_surface": "impervious_surface",
    "v2_interflow": "interflow",
}
