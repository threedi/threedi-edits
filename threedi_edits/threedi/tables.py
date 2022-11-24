# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 08:14:27 2020

@author: chris.kerklaan

#TODO:
    1. Set data constraints to look like the empty template
    1. Checks
        1. Data-checks is something allow or not-allowed in the field,
            These checks must be done to run smoothly through inpy            
        2. 
    General:
        Add functions
    3. Create extension of ogr feature here, to acces set_items and perform checks there
    4. Specify input and output of different variables.
    5. Connect get items field to fields of constants
    
Implemented checks:
    1. value type check in input feature at __setitem__ of feature, and custom in model
    
"""

# system imports
import operator
import logging

# Third-party imports
import json
import datetime as dt

# local imports
from threedi_edits import LineString
from threedi_edits.threedi.vector import ThreediVector
from threedi_edits.threedi.constants.templates import Templates
from threedi_edits.utils.project import Classes, Functions
from threedi_edits.threedi.utils.conversions import METER_TO_DEGREE

# globals
TEMPLATE = Templates()

# logger
logger = logging.getLogger(__name__)

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)


class SimpleInfiltration(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        ThreediVector.__init__(self, ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class GlobalSettings(ThreediVector):
    def __init__(self, ds, layer_name):
        ThreediVector.__init__(self, ds, layer_name)

    def check_scenario_exists(self):
        try:
            self[1]
        except IndexError:
            raise IndexError(
                f"""Scenario {self.scenario} does not exist, use self.set_scenario if a different one is needed"""
            )

    def custom_set_function(self, key, value):
        if key == "start_time":
            """start time must be in iso format"""
            if value is not None:
                try:
                    value = value.replace("/", "-")
                    dt_object = dt.datetime.fromisoformat(value)
                    value = dt_object.isoformat()
                except Exception:
                    raise ValueError(
                        "Object did not successfully convert to isoformat for start_time",
                        "is the input something like: 2014-05-14 17:42:18?",
                    )
        if key == "grid_space":
            """gridspace must be at least twice the cellsize of a raster"""
            if self.dem is not None:
                try:
                    resolution = self.dem.resolution["width"]
                    grid_space = float(self[self.scenario]["grid_space"])
                    if not ((grid_space / resolution) >= 2):
                        raise ValueError(
                            f"""
                                         Grid space must be at least twice 
                                         the resolution of the dem
                                         current grid_space {grid_space},
                                         current resolution of dem {resolution} 
                                         """
                        )
                except Exception:
                    pass
        return key, value


class AggregationSettings(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Interflow(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Groundwater(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class BoundaryConditions1d(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        """also accepts a list of tuplees"""

        if key == "timeseries":

            if type(value) == list:
                value = "\n".join([f"{_min},{ts}" for _min, ts in value])

        return key, value

    def to_dict(self, interpolate=False, path=None):
        json_bcs = []
        for bc in self:
            json_bc = {"interpolate": interpolate, "values": []}
            for values in bc["timeseries"].split("\n"):
                json_bc["values"] = [float(i) for i in values.split(",")]
            json_bcs.append((bc["connection_node_id"], json_bc))

        json_bcs_sorted = sorted(json_bcs, key=lambda x: x[0])
        json_bcs = [json_bc[1] for json_bc in json_bcs_sorted]

        if path:
            json_string = json.dumps(json_bcs)
            json_file = open(path, "w")
            json_file.write(json_string)
            json_file.close()

        return json_bcs


class BoundaryConditions2d(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Lateral1d(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Lateral2d(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class CalculationPoint(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ConnectionNodes(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value

    def custom_add_function(self, fid, items, geometry):
        if geometry:
            if len(self.spatial_filter(geometry)) > 0:
                raise ValueError(
                    f""" 
                      Trying to add a double geometry
                      feature 1: {fid}
                      feature 2: {self.spatial_filter(geometry)[0].fid}
                                 """
                )

        return fid, items, geometry

    def add_linestring(self, branches, template, NODES_STORAGE_AREA):
        """add connection nodes based on the start and end vertices of the channels"""

        for branch in branches:
            geometry = branch.geometry
            for geom, name in [
                (geometry.start_point, "start"),
                (geometry.end_point, "end"),
            ]:

                node = TEMPLATE.node
                node["the_geom"] = geom
                node["code"] = None
                node["storage_area"] = NODES_STORAGE_AREA
                geom.spatial_reference = branches.epsg

                try:
                    self.add(node)
                except ValueError:  # double connection
                    raise ValueError("found double connection")

    def create_line(self, node_id1: int, node_id2: int):
        """returns a geometry between two nodes"""

        return LineString.from_points(
            [self[node_id1].geometry.point, self[node_id2].geometry.point], epsg=4326
        )

    def reproject(self, epsg):
        return ConnectionNodes(self.reproject(epsg).ds, "v2_connection_nodes")


class ConnectedPnt(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class CrossSectionLocation(ThreediVector):
    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        # if self.channels != None:
        #     if key == "the_geom":
        #         filtered = self.channels.spatial_filter(value)
        #         if len(filtered) == 0:
        #             raise ValueError(
        #                 "Cross section location does not intersect a channel"
        #             )

        #         on_vertice = False
        #         for point in filtered[0].geometry.points:
        #             if point == value.point:
        #                 on_vertice = True
        #                 break
        #         if not on_vertice:
        #             raise ValueError(
        #                 "Cross section location geometry not on vertice"
        #             )

        return key, value

    def check_levels(self):
        """returns differences between the bank level and the reference level"""
        differences = []
        for location in self:
            differences.append(
                float(location["bank_level"]) - float(location["reference_level"])
            )
        return differences


class Culvert(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Channel(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)
        # self.csl = cross_section_locations

    def custom_set_function(self, key, value):
        return key, value

    def check_csl_presence(self):
        """checks if all channels have a cross section locations,
        returns the ids without one
        """
        channels_wo_locations = []
        for channel in self:
            if len(self.csl.filter(channel_id=channel.id)) == 0:
                channels_wo_locations.append(channel.id)
        return channels_wo_locations

    def check_start_end_point(self, minimum_distance=0.26):
        """checks if location is near the start or end point of channels"""
        # cannot be on start or end point of channel
        near_points = []
        distance = minimum_distance * METER_TO_DEGREE
        for location in self.csl:
            channel_id = int(location["channel_id"])
            channel = self[channel_id]
            near_start = (
                location.geometry.Distance(channel.geometry.start_point) < distance
            )
            near_end = location.geometry.Distance(channel.geometry.end_point) < distance

            if near_start or near_end:
                near_points.append(location.id)

        return near_points

    def check_vertices(self):
        """checks if all cross locations are on the channel vertices,
        returns the ids of locations which aer not
        """
        faulty_locations = []
        for location in self.csl:
            filtered = self.spatial_filter(location.geometry)
            if len(filtered) == 0:
                faulty_locations.append(location.id)
                logger.debug("Cross section location does not intersect a channel")
                continue

            on_vertice = False
            for point in filtered[0].geometry.points:
                if point == location.geometry.point:
                    on_vertice = True
                    break

            if not on_vertice:
                faulty_locations.append(location.id)
                logger.debug("Cross section location geometry not on vertice")

        return faulty_locations

    def closest_point_on_channel(
        self, channel_geometry, cross_geometry, minimum_distance
    ):
        """
        Checks if the geometry is located at the very start or end point of a channel,
        then it will be moved to a minimum distance
        """
        distance = minimum_distance * METER_TO_DEGREE
        near_start = cross_geometry.Distance(channel_geometry.start_point) < distance
        near_end = cross_geometry.Distance(channel_geometry.end_point) < distance
        on_vertice = cross_geometry.point in channel_geometry.points

        if near_start or near_end or not on_vertice:
            points = channel_geometry.points_on_line(
                interval=distance, start=False, end=False
            )

            closest = []
            for idx, point in enumerate(points):
                closest.append((idx, point.Distance(cross_geometry)))

            return points[sorted(closest, key=lambda x: x[1])[0][0]]
        else:
            return cross_geometry

    def non_existing_cross_section_locations(self):
        """returns a list of not existing locations"""
        channel_fids = self.fids
        return [l.id for l in self.csl if l["channel_id"] not in channel_fids]

    def predict_calculation_points(self, feature_id=None, projection=28992):
        """returns calculation points per feature"""
        if feature_id:
            feature = self[feature_id]
            geometry = feature.geometry
            interval = int(feature["dist_calc_points"])

            projected = LineString(geometry.reproject(projection))

            line_points = projected.points_on_line(interval, start=False, end=False)
            if len(line_points) == 0:
                line_points = [projected.start_point]

            calc_points = []
            for i, point in enumerate(line_points):
                fid = i + 1
                code = f"{feature_id}_{fid}_{feature['calculation_type']}"
                calc_points.append(
                    {
                        "content_type_id": feature_id,
                        "user_ref": f"{code}#{feature_id}#v2_channel#{fid}",
                        "calc_type": feature["calculation_type"] - 100,
                        "the_geom": point.reproject(4326),
                    }
                )
            return calc_points


class Pipe(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Weir(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Orifice(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Manhole(ThreediVector):
    constraints = {"drain_level": {"comparison": (operator.ge, "bottom_level")}}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value

    def set_drain_level(self, below_surface_level=0.15, none_only=True):
        """sets all drain levels with None to 0.15 m below surface level"""

        for feature in self:
            drain_level = feature["drain_level"]
            if none_only and feature["drain_level"] != None:
                continue

            feature["drain_level"] = feature["surface_level"] - below_surface_level

            self.edits.append(
                (
                    feature.fid,
                    "drain_level",
                    drain_level,
                    feature["drain_level"],
                )
            )


class Pumpstation(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Windshielding(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class CrossSectionDefinition(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    # def custom_plot(self, feature):
    #     """this is the only way to alter this function"""

    #     def plot(filename=None):
    #         width = [float(i) for i in feature.items["width"].split(" ")]
    #         height = [float(i) for i in feature.items["height"].split(" ")]
    #         middle = max(width) / 2

    #         # step 6: make the profile ready to be plotted
    #         yz_tab_list = []
    #         for h, w in zip(height, width):
    #             if h == 0:
    #                 yz_tab_list.append([middle, h])
    #             else:
    #                 yz_tab_list.append([middle - w / 2, h])
    #                 yz_tab_list.append([middle + w / 2, h])

    #         yz_tab = np.array(yz_tab_list)
    #         yz_tab_left = yz_tab[yz_tab[:, 0] < middle]
    #         yz_tab_right = yz_tab[yz_tab[:, 0] >= middle]
    #         yz_tab_left_sorted = yz_tab_left[
    #             np.lexsort((yz_tab_left[:, 0], -yz_tab_left[:, 1]))
    #         ]
    #         yz_tab_right_sorted = yz_tab_right[
    #             np.lexsort((yz_tab_right[:, 0], yz_tab_right[:, 1]))
    #         ]
    #         yz_tab_sorted = np.concatenate(
    #             (yz_tab_left_sorted, yz_tab_right_sorted), axis=0
    #         )

    #         plt.figure(figsize=(10, 5))
    #         plt.plot(
    #             yz_tab_sorted[:, 0],
    #             yz_tab_sorted[:, 1],
    #             label="Tabulated trapezium",
    #         )
    #         plt.xlabel("Distance from bank [m]", fontdict=None, labelpad=None)
    #         plt.ylabel("Depth [m]", fontdict=None, labelpad=None)
    #         plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    #         if filename:
    #             plt.savefig(filename)
    #         plt.show()

    #     return plot

    def custom_set_function(self, key, value):
        if (key == "width" or key == "height") and type(value) == list:
            value_str = ""
            for i, v in enumerate(value):
                if i == len(value) - 1:
                    value_str = value_str + f"{v}"
                else:
                    value_str = value_str + f"{v} "
            value = value_str

        return key, value


class Obstacle(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Levee(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class GridRefinement(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class GridRefinementArea(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        # if key == "refinement_level":
        #     if self.global_setting is not None:
        #         if value > self.global_setting['kmax']:
        #             raise ValueError("Refinment level above kmax")

        return key, value


class NumericalSettings(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class DemAverageArea(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ImperviousSurfaceMap(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ImperviousSurface(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class SurfaceParameters(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class SurfaceMap(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Surface(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(key, value):
        return key, value


class Control(ThreediVector):

    constraints = {}
    required = {
        "control_group": None,
        "control_table": None,
        "control_measure_map": None,
        "control_measure_group": None,
    }

    def __init__(
        self,
        ds,
        layer_name,
        control_group=None,
        control_table=None,
        control_measure_group=None,
        control_measure_map=None,
    ):
        super().__init__(ds, layer_name)
        self.control_group = control_group
        self.control_table = control_table
        self.control_measure_group = control_measure_group
        self.control_measure_map = control_measure_map

    def custom_set_function(self, key, value):
        return key, value

    @property
    def cview(self):
        if not hasattr(self, "control_view"):
            control_view = {}
            for control in self:
                control_view[control.id] = {}

                control_view[control.id]["control"] = control.items
                control_id = control["control_id"]

                # measure groups
                mappings = self.control_measure_map.filter(
                    measure_group_id=control["measure_group_id"]
                )
                control_view[control.id]["measure_group"] = []
                for mapping in mappings:
                    control_view[control.id]["measure_group"].append(
                        {mapping.id: mapping.items}
                    )
                control_view[control.id]["control_table"] = {
                    control_id: self.control_table[control_id].items
                }
                self.control_view = control_view
        return self.control_view


class ControlDelta(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlGroup(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlMeasureGroup(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlMeasureMap(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlMemory(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlPid(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlTable(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class ControlTimed(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class FloodFill(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class Outlet(ThreediVector):

    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class PumpedDrainageArea(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


class InitialWaterlevel(ThreediVector):
    constraints = {}

    def __init__(self, ds, layer_name):
        super().__init__(ds, layer_name)

    def custom_set_function(self, key, value):
        return key, value


MODEL_MAPPING = {
    "v2_simple_infiltration": SimpleInfiltration,
    "v2_global_settings": GlobalSettings,
    "v2_aggregation_settings": AggregationSettings,
    "v2_interflow": Interflow,
    "v2_initial_waterlevel": InitialWaterlevel,
    "v2_groundwater": Groundwater,
    "v2_1d_boundary_conditions": BoundaryConditions1d,
    "v2_2d_boundary_conditions": BoundaryConditions2d,
    "v2_1d_lateral": Lateral1d,
    "v2_2d_lateral": Lateral2d,
    "v2_connection_nodes": ConnectionNodes,
    "v2_cross_section_location": CrossSectionLocation,
    "v2_connected_pnt": ConnectedPnt,
    "v2_calculation_point": CalculationPoint,
    "v2_culvert": Culvert,
    "v2_channel": Channel,
    "v2_pipe": Pipe,
    "v2_weir": Weir,
    "v2_orifice": Orifice,
    "v2_manhole": Manhole,
    "v2_pumpstation": Pumpstation,
    "v2_windshielding": Windshielding,
    "v2_cross_section_definition": CrossSectionDefinition,
    "v2_obstacle": Obstacle,
    "v2_levee": Levee,
    "v2_grid_refinement": GridRefinement,
    "v2_grid_refinement_area": GridRefinementArea,
    "v2_numerical_settings": NumericalSettings,
    "v2_dem_average_area": DemAverageArea,
    "v2_impervious_surface_map": ImperviousSurfaceMap,
    "v2_impervious_surface": ImperviousSurface,
    "v2_surface_parameters": SurfaceParameters,
    "v2_surface_map": SurfaceMap,
    "v2_surface": Surface,
    "v2_control": Control,
    "v2_control_delta": ControlDelta,
    "v2_control_group": ControlGroup,
    "v2_control_measure_group": ControlMeasureGroup,
    "v2_control_measure_map": ControlMeasureMap,
    "v2_control_memory": ControlMemory,
    "v2_control_pid": ControlPid,
    "v2_control_table": ControlTable,
    "v2_control_timed": ControlTimed,
    "v2_floodfill": FloodFill,
    "v2_outlet": Outlet,
    "v2_pumped_drainage_area": PumpedDrainageArea,
    "v2_calculation_point": CalculationPoint,
}
