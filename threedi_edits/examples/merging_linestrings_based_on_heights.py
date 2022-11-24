# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:15:34 2021

@author: chris.kerklaan

Example of sampling linestrings

"""

# imports
import os
from osgeo import ogr
import numpy as np
from threedi_edits import Vector, Raster, Progress, LineString, RasterGroup


def split(levees: Vector, distance=5):
    """splitting linestrings with a certain distance"""
    # first cut into pieces of 5 meter and sample some dem values
    split = levees.copy(shell=True)
    for levee in Progress(levees, "Splitting"):
        # create single levees due to multilinestrings
        for single_levee in levee.geometry.to_single():
            # split the levee in sections of 5 meter
            levee_lines = LineString(single_levee).split(distance)
            for levee_line in levee_lines:
                split.add(geometry=levee_line, items=levee.items)
    return split


def sample(levees: Vector, dem: Raster, distance=2, perpendicular_length=10):
    """
    Sample levees on a dem. Every distance a perpendicuar length is taken.
    The maxium on every perpendicular linestring is taken an averaged.

    params:
        levees: Vector
        dem: Raster or Rastergroup
        distance: distance between perpendicular lines
        perpendicular_lenhgth: Length of the perpendicular linestrings on
        which the sample is taken.

    """

    sampled = levees.copy(shell=True)
    sampled.add_field("height", float)
    for levee in Progress(levees, "Sampling"):
        lines = levee.geometry.perpendicular_lines(distance, perpendicular_length)
        heights = []
        for line in lines:
            line_data = dem.read(line)
            try:
                if len(line_data) > 0:
                    heights.append(np.nanmax(line_data))
            except Exception as e:
                print(e)
                pass
        mean = float(np.nanmean(heights))
        sampled.add(items=levee.items, height=mean, geometry=levee.geometry)

    return sampled


class LeveeAggregator:
    def __init__(
        self,
        levees: Vector,
        height_field: str = "height",
        mode: str = "min_max",
        threshold: float = 0.25,
        type_field: str = None,
        thresholds: dict = None,
        feature: bool = False,
        feature_dem: Raster = None,
        feature_split_distance: int = 10,
        feature_perpendicular_distance: int = 2,
        feature_perpendicular_length: int = 10,
    ):
        """
        Class to aggregate levees based on a field.
        Using 'type_field' it is possible to run the algorithm
        per unique field value. 'mode' sets the method in which the levees
        are aggregated.

        If the aggregator is run per feature, it has to split and sample the features
        on the go. So, the levee input does not have to be splitted and sampled which
        is the case if the aggregator is run differently.

        There are generally two ways in which  it is run:

           Without fields or worrying about features.
           ------
           aggregator = LeveeAggregator(levees, height_field='height', threshold=0.25)
           output = aggregator.run()
           output.write(output_path)

           With fields and not worrying about features.
           ------
           thresholds = {
               "kering-primair": 0.25,
               "kering-secundair": 0.25,
               "weg": 0.5,
               "spoor": 0.5,
           }
           aggregator = LeveeAggregator(levees, type_field='type', height_field='height', thresholds=thresholds)
           output = aggregator.run()
           output.write(output_path)

           With fields and using it per feature.
           ------
           thresholds = {
               "kering-primair": 0.25,
               "kering-secundair": 0.25,
               "weg": 0.5,
               "spoor": 0.5,
           }
           aggregator = LeveeAggregator(levees, type_field='type', height_field='height', thresholds=thresholds,
                                        feature=True,
                                        feature_dem=dem,
                                        feature_split_distance,
                                        feature_perpendicular_distance,
                                        feature_perpendicular_length)
           output = aggregator.run()
           output.write(output_path)




        params:
            levees: Vector with LineStrings, with or without height field
            height_field: Field in which the height is stored
            type_field: Field in which the vector gets separated.
            mode: Can either be 'min_max' or median.

            min_max:
            1. The adjacent levee height is not lower than the threshold minus the
            minium of the current heights.
            2. The adjecent levee height is not lower than the threshold plus the
            maximum of the current heights

            median:
            1. The difference between the median of the current heights
            and the next levee height is not larger than the threshold.

            feature: levees are converged per feature.
            feature_dem: Raster on which the levees are sampled,
            feature_split_distance: distance on which the levees are split,
            feature_perpendicular_distance: distance on wici perpedicular lines are created,
            feature_perpendicular_length: length of the perpendicular lines, we sample on this.

        """

        # making sure there are only np.nan's
        if not feature:
            data = levees.copy(shell=True)
            for levee in levees:
                items = levee.items
                if items[height_field] is None:
                    items[height_field] = np.nan
                data.add(items=items, geometry=levee.geometry)
        else:
            data = levees.copy()

        # separate the data if type_field is used
        if type_field:
            field_data = {}
            for value in list(set(data.table[type_field])):
                field_data[value] = data.filter(
                    return_vector=True, **{type_field: value}
                )
            self.field_data = field_data

        # adding to self to use in other functions
        self.data = data
        self.mode = mode
        self.threshold = threshold
        self.thresholds = thresholds
        self.feature = feature
        self.type_field = type_field
        self.dem = feature_dem
        self.feature_split_distance = (feature_split_distance,)
        self.feature_perpendicular_distance = feature_perpendicular_distance
        self.feature_perpendicular_length = feature_perpendicular_length

    def add_fields(self, output):
        output.add_field("m_count", int)
        output.add_field("m_aheight", float)
        output.add_field("m_max", float)
        output.add_field("m_min", float)
        output.add_field("height", float)
        return output

    def run(self):
        """run with standard settings"""
        if self.type_field:
            return self._run_fields(
                threshold=self.threshold, thresholds=self.thresholds
            )
        return self._run(self.data, self.threshold)

    def _run_fields(self, threshold=0.25, thresholds={}):
        field_output = self.data.copy(shell=True)
        field_output = self.add_fields(field_output)
        for field, data in self.field_data.items():
            if thresholds:
                threshold = thresholds[field]

            print("Running:", field, "Threshold:", threshold, "Input count", data.count)
            if self.feature:
                output = self._run_feature(
                    data,
                    self.dem,
                    threshold,
                    self.feature_split_distance,
                    self.feature_perpendicular_distance,
                    self.feature_perpendicular_length,
                )
            else:
                output = self._run(data, threshold)

            print("Output count:", output.count)
            for feature in output:
                field_output.add(feature)
        return field_output

    def _run(self, data, threshold=0.25):
        output = data.copy(shell=True)
        output = self.add_fields(output)

        data = data.copy()
        data.reset()

        levee = next(data)
        levee_geometry = levee.geometry
        levee_items = levee.items
        heights = [levee["height"]]

        data.delete(levee)

        while len(data) != 0:
            print(len(data), "levee_id:", levee.id)
            data.reset()
            levee_geometry, added_levees = add_adjacent_levees(
                heights, levee_geometry, data, threshold, self.mode
            )

            # add to output if nothing has to be added anymore
            if len(added_levees) == 0:

                output.add(
                    geometry=levee_geometry,
                    items=levee_items,
                    m_count=len(heights),
                    m_aheight=np.nanmean(heights),
                    m_max=np.nanmax(heights),
                    m_min=np.nanmin(heights),
                )

                levee = next(data)  # retrieves a new feature even if some are deleted
                levee_items = levee.items
                levee_geometry = levee.geometry
                heights = [levee["height"]]
                data.delete(levee)

            else:
                # delete from data and to heights
                for levee_id in added_levees:
                    added_levee = data[levee_id]
                    heights.append(added_levee["height"])
                    data.delete(added_levee)

        data.close()

        # finally
        output.add(
            geometry=levee_geometry,
            items=levee_items,
            m_count=len(heights),
            m_aheight=np.nanmean(heights),
            m_max=np.nanmax(heights),
            m_min=np.nanmin(heights),
        )

        return output

    def _run_feature(
        self,
        levees,
        dem,
        threshold,
        split_distance=10,
        perpendicular_distance=2,
        perpendicular_length=10,
    ):
        """the same thing as above but then per feature and splitting and sampling at the same time"""

        output = self.data.copy(shell=True)
        output = self.add_fields(output)

        # copy the original
        # get the first params and delete the first levee:
        #   levee_geometry
        #   levee_items
        #   heights

        # do it per feature
        for levee in levees:

            print(levee.id)
            # split the levee into sections
            print("Splitting")
            split_lines = []
            for single_levee in levee.geometry.to_single():
                split_lines.extend(LineString(single_levee).split(split_distance))

            # sample the data on the dem
            print("Sampling")
            sampled_lines = []
            for line in Progress(split_lines, "Sampling"):

                heights = []
                lines = line.perpendicular_lines(
                    perpendicular_distance, perpendicular_length
                )
                for line in lines:
                    line_data = dem.read(line)
                    try:
                        if len(line_data) > 0:
                            heights.append(np.nanmax(line_data))
                    except Exception as e:
                        print(e)
                        pass

                sampled_lines.append(float(np.nanmean(heights)))

            # merge into a vector
            print("Merging")
            vector = output.copy(shell=True)
            for geometry, height in zip(split_lines, sampled_lines):
                vector.add(items=levee.items, geometry=geometry, height=height)

            # run the convergence
            print("Converging")
            feature_output = self._run(vector, threshold)
            for feature in feature_output:
                output.add(feature)

        return output


def add_adjacent_levees(
    heights: list,
    levee_geometry: ogr.Geometry,
    levee_vector: Vector,
    threshold=0.25,
    mode="min_max",
):
    """
    adds adjacent levees if based on a list of previously generated
    heights and the height of the adjecent levee.
    The adjacent levee and its height are added in the following cases:

        min_max:
        1. The adjacent levee height is not lower than the threshold minus the
        minium of the current heights.
        2. The adjecent levee height is not lower than the threshold plus the
        maximum of the current heights

        median:
        1. The difference between the median of the current heights
        and the next levee height is not larger than the threshold.

    other params:
        heights: list of heights
        levee_geometry: geometry of merged levees
        levee_vector: levee vector


    """

    output_geometry = levee_geometry.copy()
    added_levees = []
    for levee in levee_vector.spatial_filter(levee_geometry.buffer(0.0001)):
        next_levee_height = levee["height"]

        is_nan = np.isnan(next_levee_height)
        if not is_nan:
            if mode == "min_max":
                minimum_difference = abs(next_levee_height - min(heights)) < threshold
                maximum_difference = abs(next_levee_height - max(heights)) < threshold
                add = minimum_difference and maximum_difference
            elif mode == "median":
                median = np.percentile(heights, 50)
                add = abs(median - next_levee_height) < threshold
        else:
            add = False

        if add or is_nan:
            added_levees.append(levee.id)
            output_geometry = output_geometry.union(levee.geometry)
            output_geometry = output_geometry.dissolve()

    return output_geometry, added_levees


if __name__ == "__main__":
    # from here the code is run
    os.chdir(
        r"C:\Users\chris.kerklaan\Documents\Projecten\basis kaarten hhnk\processing"
    )

    # load the levees
    levees = Vector("levees/all_levees.gpkg")

    # load the dem as a rastergroup, because they are in tiles
    dem = RasterGroup([Raster(f"data/dtm/{i}") for i in os.listdir("data/dtm")])

    # step 1 - split the levees
    split_levees = split(levees, distance=10)
    split_levees.write("levees/split_10m.gpkg")

    # step 2 - sample the levees on the dem
    sampled_levees = sample(split_levees, dem, distance=2, perpendicular_length=10)
    sampled_levees.write("levees/sampled_10m.gpkg")

    # step 3 - run the aggregator
    levees.add_field("height", float)
    aggregator = LeveeAggregator(levees, height_field="height", mode="min_max")
    levees_aggregated = aggregator.run(aggregator.data, 0.25)
    levees_aggregated.write("levees/aggregated.gpkg")

    # or everything at the same time
    thresholds = {
        "kering-primair": 0.25,
        "kering-secundair": 0.25,
        "weg": 0.5,
        "spoor": 0.5,
    }

    aggregator = LeveeAggregator(
        levees,
        type_field="type",
        mode="median",
        thresholds=thresholds,
        feature=True,
        feature_dem=dem,
    )
    levees_aggregated = aggregator.run()
    levees_aggregated.write("levees/aggregated_per_feature2.gpkg")
