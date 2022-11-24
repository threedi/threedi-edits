# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:13:41 2022

@author: chris.kerklaan

This scripts creates a plot along a model of a channel.

For connection node the results ore sampled on a waterlevel, then the watersampled are put on  aline.



1. Create a view 


"""

from threedigrid.admin.gridresultadmin import GridH5ResultAdmin
from threedigrid.admin.gridadmin import GridH5Admin
import threedi_edits as tre
import numpy as np
import matplotlib.pyplot as plt

# globals
MIN_LEVEE_HOOGTE = 6.31


class Grid:
    def __init__(self, results_folder, start_time=0, end_time=None):
        self.h5_file_path = results_folder + "/gridadmin.h5"
        self.netcdf_file_path = results_folder + "/results_3di.nc"
        self.gr = GridH5ResultAdmin(self.h5_file_path, self.netcdf_file_path)
        self.ga = GridH5Admin(self.h5_file_path)

        self.start_time = start_time
        self.end_time = end_time
        if not end_time:
            self.end_time = self.gr.nodes.timestamps[-1]

    @property
    def nodes(self):
        if not hasattr(self, "_nodes"):
            self.set_nodes()

        return self._nodes

    def set_nodes(self):
        self._nodes = tre.Vector.from_scratch("nodes", 1, 28992)
        coordinates = np.vstack(self.gr.nodes.coordinates.T)
        ids = self.gr.nodes.id
        for id, coordinate in zip(ids, coordinates):
            self._nodes.add(geometry=tre.Point.from_point(coordinate), fid=id)

    def level(self, grid_id=[], model_id=[], start_time=0, end_time=None):
        """retrieves the waterlevel for different either model ids or grid ids
        params:
            grid_id: grid id, list or integer
            model_id: model id, list or integer
            start_time: start time, as a datetime object
            end_time: end time, as a datetime object
        """

        if end_time:
            self.end_time = end_time

        if type(grid_id) == int:
            return list(
                self.gr.nodes.filter(id=grid_id)
                .timeseries(self.start_time, self.end_time)
                .s1.flatten()
            )

        if type(model_id) == int:
            return list(
                self.gr.nodes.filter(content_pk=model_id)
                .timeseries(self.start_time, self.end_time)
                .s1.flatten()
            )

        output_dict = {}
        if len(grid_id) > 0:
            for gid in grid_id:
                output_dict[gid] = list(
                    self.gr.nodes.filter(id=gid)
                    .timeseries(self.start_time, self.end_time)
                    .s1.flatten()
                )

        if len(model_id) > 0:
            for mid in model_id:
                output_dict[mid] = list(
                    self.gr.nodes.filter(content_pk=mid)
                    .timeseries(self.start_time, self.end_time)
                    .s1.flatten()
                )

        return output_dict


def create(
    subject,
    result_folder,
    schematisation_path,
    selection_shape,
    output_name,
    start_time=43200,
    start_node_id=74507,
):

    results = Grid(result_folder, start_time=start_time, end_time=start_time + 350)

    # First select connection nodes with shape
    shape = tre.Vector(selection_shape).reproject(4326)[0].geometry

    # Open schematisation
    schema = tre.ThreediEdits(schematisation_path)

    node_fids = schema.nodes.spatial_filter(shape, "Within", return_fid=True)

    output = {"nodes": node_fids, "waterstand": [], "oeverhoogte": []}
    levels = results.level(model_id=node_fids)

    for node_id, level in levels.items():
        output["waterstand"].append(level[0])

    # now we find the bank levels
    channel = schema.channels.spatial_filter(shape, "Within", return_vector=True)

    # reproject
    channel = channel.copy().reproject(28992)
    levees = schema.levees.copy().reproject(28992)

    for node_id in node_fids:
        buffered_node = schema.nodes[node_id].geometry.buffer(0.00001)
        buffered_node = buffered_node.reproject(28992)
        clipped_channel = channel.clip(geometry=buffered_node).dissolve(quiet=True)
        line = clipped_channel[0]

        if line.geometry.IsEmpty():
            output["oeverhoogte"].append([None, None])
            continue

        perpendicular_line = line.geometry.perpendicular_lines(0.5, 20)[0]
        levee = levees.spatial_filter(perpendicular_line, return_vector=True)
        if levee.count > 0:
            if levee.count == 2:
                output["oeverhoogte"].append(
                    [
                        levee[levee.fids[0]]["crest_level"],
                        levee[levee.fids[1]]["crest_level"],
                    ]
                )
                print(levee.fids, node_id)
            else:
                output["oeverhoogte"].append(
                    [levee[levee.fids[0]]["crest_level"], None]
                )

        else:
            output["oeverhoogte"].append([None, None])

    # reorder from distance to gemaal
    start_node = schema.nodes[start_node_id].geometry

    distances = []
    for node_id in output["nodes"]:
        distances.append((schema.nodes[node_id].geometry.Distance(start_node), node_id))

    ordered_output = {"nodes": [], "waterstand": [], "oeverhoogte": [], "afstand": []}
    for (d, node_id) in sorted(distances):
        index = output["nodes"].index(node_id)
        ordered_output["nodes"].append(str(node_id))
        ordered_output["oeverhoogte"].append(output["oeverhoogte"][index])
        if output["waterstand"][index] is not None:
            ordered_output["waterstand"].append(output["waterstand"][index])
        else:
            ordered_output["waterstand"].append(output["waterstand"][index])

        ordered_output["afstand"].append(d / 0.00001 / 1000)

    # create data
    # delete nodata
    deletes = []
    for i, ws in enumerate(ordered_output["waterstand"]):
        if ws == -9999:
            deletes.append(i)

    for did in deletes:
        del ordered_output["nodes"][did]
        del ordered_output["waterstand"][did]
        del ordered_output["oeverhoogte"][did]
        del ordered_output["afstand"][did]

    oeverhoogtes = []
    for i in ordered_output["oeverhoogte"]:
        for hoogte in i:
            if hoogte is not None:
                oeverhoogtes.append(hoogte)

    minimale_oeverhoogte = min(oeverhoogtes)

    # plot lines
    fig, axs = plt.subplots(2, figsize=(10, 10))

    fig.suptitle(f"Verhanglijn en oeverhoogtes van het Overijssels Kanaal {subject}")

    axs[0].plot(
        ordered_output["afstand"],
        ordered_output["waterstand"],
        label=f"Waterstand na {int(start_time/3600)} uur, ~{int(start_time/3600/24)} dagen",
    )
    axs[0].plot(
        ordered_output["afstand"],
        [minimale_oeverhoogte for i in range(len(ordered_output["afstand"]))],
        label="Laagste oeverhoogte",
    )

    axs[0].set_xlabel("Afstand vanaf gemaal Ankersmit [kilometer] ")
    axs[0].set_ylabel("Waterstand [mNAP]")

    axs[0].scatter(
        ordered_output["afstand"],
        [i[0] for i in ordered_output["oeverhoogte"]],
        color="red",
        label="Oeverhoogte Zuid Oost",
        s=5,
    )

    axs[0].scatter(
        ordered_output["afstand"],
        [i[1] for i in ordered_output["oeverhoogte"]],
        color="blue",
        label="Oeverhoogte Noord West",
        s=5,
    )
    axs[0].legend()
    axs[0].grid()

    axs[1].plot(
        ordered_output["afstand"],
        ordered_output["waterstand"],
        label=f"Waterstand na {int(start_time/3600)} uur, ~{int(start_time/3600/24)} dagen ",
    )
    axs[1].set_xlabel("Afstand vanaf gemaal Ankersmit [kilometer] ")
    axs[1].set_ylabel("Waterstand [mNAP]")

    axs[1].legend()
    axs[1].grid()

    fig.show()
    fig.savefig(output_name, dpi=400)


if __name__ == "__main__":
    schematisation_path = r"C:\Users\chris.kerklaan\Documents\threedi_models\NBW Overijssels Kanaal - winter GHG - NVO\work in progress\schematisation/NBW Overijssels Kanaal - winter GHG - NVO.sqlite"
    result_folder = r"C:\Users\chris.kerklaan\Documents\threedi_models\NBW Overijssels Kanaal - winter GHG - NVO\work in progress\results\sim_86859_T100_winter_GHG_-_NVO_(n_=_0.3)"
    selection_shape = (
        r"C:/Users/chris.kerklaan/Documents/Projecten/overijssels_kanaal/data/ovk.shp"
    )
    output_name = r"C:\Users\chris.kerklaan\Documents\Projecten\overijssels_kanaal\processing/sideview_NVO_manning_03.png"
    create(
        "NVO (Manning n=0.3)",
        result_folder,
        schematisation_path,
        selection_shape,
        output_name,
        start_time=120 * 3600,
        start_node_id=74507,
    )
