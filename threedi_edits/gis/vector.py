"""
Created on Fri Apr 26 16:39:50 2019

@author: chris.kerklaan

Vector wrapper used for GIS programming

Goals:
    1. Provide a better interface for ogr
    2. Provides a table and spatial index which grows and deletes with adding
    and deleting features
    3. Provides basic geoprocessing options: Clipping, dissolving etc.
    4. Provides appropriate fixing of geometries

TODO for rewriting in shapely tre version 1.0:
    1. Make shapely available in tre's geometry functions.
        1.1. Add svg in geometries. #Done
        1.2. Convert easily between shapely and ogr. #Done
    2. Rewrite functions below vector.py so that they use the Vector class.
        2.0. Rewrite vectorbase so that it is as fast as possible.
        2.1. Check which functions are used within other functions.
        2.2 Voronoi can be rewritten with the multipoint voronoi.
        
    3. Rewrite the table and table features with numpy.
    4. Add @cached_property
    5. Accelerate table generation.
    6. We can remove rtree indexing as it is already present within shapely/pygeos.
    6. Finally do the testing.
    
Bugs:
    1. Multipart2singlepart does not write as single part file shape

Note:
    21. Ogr indices are created when copying the vector

Notes on adding shapely:
    There are multiple reasons for adding shapely:
        - Speed-ups for multi-geometry intersections.
        - Extra geometry operations.
    We still need to have a fast way to use the shapely functions.
    
        

Notes on pointers:
    Datasource is your pointer, pointer always needs to be present when
    data is edited in layer
    When a feature is created the feature is the pointer.

    Due to this pointers, gdal ds, layers and features are kept in place,
    hence the classs are not an extension of the object but built on top of them.
    If working with pure ogr, make sure that pointer are available:

        lines = Vector(TEST_DIRECTORY+"lines.shp")
        works:
            layer = lines.layer
            feature=layer[0]
            feature.geometry().Clone()
        does not work:
            layer = lines.layer
            feature.geometry().Clone()

    Also:
        works:
            for geometry in multipartgeometry:
                geometry.GetArea()
        does not work:
            for geometry in multipartgeometry:
                function(geometry)

            def function(geometry):
                geometry.GetArea()

    Functions outside the class can only use ogr-functions
    
Notes on gdal speed-ups
    - Function with 'Get' are the fastest gdal functions
        (GetFeature, GetGeometry) instead of indexing '[]' and .geometry() and .items()
        Use field_indices and feature_items to bypass and increase speed

"""
# First-party imports
import logging
import importlib
import numpy as np

# from multiprocessing import Pool

# Third party imports
from osgeo import ogr, gdal

HAS_SCIPY = importlib.util.find_spec("scipy") is not None
if HAS_SCIPY:
    from scipy import spatial

# Local imports

from threedi_edits.gis.geometry import (
    SINGLE_TYPES,
    SINGLE_TO_MULTIPLE,
    MULTI_TYPES,
    clip as clip_geometry,
    difference as difference_geometry,
    dissolve as dissolve_geometry,
)

from threedi_edits.gis.polygon import centerline as centerline_geometry
from threedi_edits.gis.polygon import Polygon
from threedi_edits.gis.vectorbase import VectorBase, mem_path
from threedi_edits.utils.project import Classes, Functions
from threedi_edits.utils.progress import Progress

# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)


# Global DRIVERS
DRIVER_GDAL_MEM = gdal.GetDriverByName("MEM")
DRIVER_OGR_MEM = ogr.GetDriverByName("Memory")
DRIVER_OGR_GPKG = ogr.GetDriverByName("GPKG")
DRIVER_OGR_SHAPEFILE = ogr.GetDriverByName("ESRI Shapefile")
DRIVERS_SQL = ["PostgreSQL", "GPKG"]

# Memory counter
_mem_num = 0

# Logger
logger = logging.getLogger(__name__)


class Vector(VectorBase):
    """wrapper of ogr layer, input is ogr datasource, only one layer
    is expected
    """

    instance = "gis.vector.Vector"

    def __init__(self, path, **kwargs):
        super().__init__(path, **kwargs)

    def copy(self, geometry_type=None, shell=False, fields=True, driver="GPKG"):
        """
        Copies a vector.
        params:
            geometry_type: change the geometry type
            shell: Does not add the features
            fields: Set to True is fields have to be included
        """
        return Vector(
            path=None, datasource=self._copy(geometry_type, shell, fields, driver)[0]
        )

    def subset(self, fids):
        """Returns a subset of features based on a list of fids"""
        return Vector(self._subset(*fids))

    def spatial_subset(
        self,
        geometry: ogr.Geometry,
        method: str = "intersects",
    ):
        """
        Spatially filtered geometries using the shapely predicates.
        params:
            geometry: ogr.Geometry, shapely geometry, te.Geometry
            method: Use the shapely predicates.

        """
        if hasattr(geometry, "ds"):
            return self.full_spatial_filter(geometry, method)

        return Vector(self._subset(*self._sfilter(geometry, method)))

    def filter(self, **filtering):
        """
        Filters the dataset on fields and values of those fields.
        E.g., percentile=1 or class=1 or both.
        """
        return Vector(self._filter(**filtering))

    def spatial_filter(
        self,
        geometry: ogr.Geometry,
        method: str = "intersects",
    ):
        """
        Spatially filtered geometries using the shapely predicates.
        params:
            geometry: ogr.Geometry, shapely geometry, te.Geometry
            method: Use the shapely predicates.
        """
        print("deprecated, use spatial_subset")

        if hasattr(geometry, "ds"):
            return self.full_spatial_filter(geometry, method)

        return Vector(self._subset(*self._sfilter(geometry, method)))

    def add(self, feature=None, geometry=None, fid=-1, **items):
        """adds a feature"""
        return self._add(feature, geometry, fid, **items)

    def full_spatial_subset(self, vector, method="intersects", quiet=True):
        """
        Spatially filteres with a whole vector.
        """
        fids = []
        for geometry in Progress(vector.geometries, f"{method} filter", quiet):
            fids.extend(self._sfilter(geometry, method))

        return Vector(self._subset(*fids))

    def buffer(self, buffer_size: float, quiet: bool = True):
        """Buffers a layer, output is always a polygon"""
        target = self.copy(geometry_type=ogr.wkbPolygon, shell=True)
        for feature in Progress(self, "Buffer", quiet):
            buffered = feature.geometry.Buffer(buffer_size)
            target.add(geometry=buffered, items=feature.items)
        return target

    def simplify(self, simple_factor: float, quiet: bool = True):
        """Simplifies all geometries."""
        target = self.copy(geometry_type=ogr.wkbPolygon, shell=True)
        for feature in Progress(self, "Simplify", quiet):
            simple = feature.geometry.Simplify(simple_factor)
            target.add(geometry=simple, items=feature.items)
        return target

    def centroid(self, quiet: bool = True):
        """Creates a centroid of all geometries."""
        target = self.copy(geometry_type=ogr.wkbPoint, shell=True)
        for feature in Progress(self, "Centroid", quiet):
            centroid = feature.geometry.Centroid()
            target.add(geometry=centroid, items=feature.items)
        return target

    def to_single(self, quiet: bool = True):
        """Converts multiparts to single parts"""
        target = self.copy(geometry_type=SINGLE_TYPES[self.geometry_type], shell=True)
        for feature in Progress(self, "To singles", quiet):
            geometry_type = feature.geometry.type
            items = feature.items
            if geometry_type > 3:
                for part in feature.geometry:
                    target.add(items=items, geometry=part)
            else:
                target.add(items=items, geometry=feature.geometry)
        return target

    def reproject(self, epsg: int, quiet=True):
        """Takes an ogr layer, reprojects and returns an ogr datasource"""
        target = Vector.from_scratch("reprojected", self.geometry_type, epsg)
        for field_name, field_type in self.fields.items.items():
            target.add_field(field_name, field_type[0])

        for feature in Progress(self, "Reproject", quiet):
            geometry = feature.geometry.reproject(epsg)
            target.add(geometry=geometry, items=feature.items)

        return target

    def centerlines(self, resolution, quiet: bool = True):
        """creates centerlines from polygons"""
        target = self.copy(geometry_type=ogr.wkbLineString, shell=True)
        for feature in Progress(self, "Centerline", quiet):
            centerline = centerline_geometry(feature.geometry, resolution)
            target.add(geometry=centerline, items=feature.items)
        return target

    def polygon_to_lines(self, quiet=True):
        """Takes a polygon and returns a linestring"""
        target = self.copy(geometry_type=ogr.wkbLineString, shell=True)
        for feature in Progress(self, "Centerline", quiet):
            target.add(geometry=feature.geometry.Boundary(), items=feature.items)
        return target

    def fix(self, geometry_type=None, quiet=True):
        """
        Fixes the geometries of a layer.
        Invalid geometries but also 3D features will be
        set to the layers' geometry type.
        params:
            geometry_type: Sets the target geometry type
                1 = point
                2 = linestring
                3 = polygon
                4 = multipoint
                5 = multilinestring
                6 = multipolygon
            quiet: Shows the progress

        """
        if not geometry_type:
            geometry_type = self.geometry_type

        target = self.copy(geometry_type=geometry_type, shell=True)
        for feature in Progress(self, "geometry fixes", quiet):
            target.add(
                geometry=feature.geometry.fix(geometry_type), items=feature.items
            )
        return target

    def clip(
        self,
        mask=None,
        geometry: ogr.Geometry = None,
        quiet: bool = True,
    ):
        """
        Clip off areas that are not covered by the input vector.
        If the feature falls outside of the masks, the will not be included.
        We are skipping features without geometries.


        params:
            mask: Vector on which it is clipped.
            quiet: Shows the progress
            use_ogr: uses the ogr variant of this function
        """
        if geometry:
            mask = Vector.from_scratch("", 3, geometry.epsg)
            mask.add(geometry=geometry)

        target = self.copy(shell=True)
        type_geometry = target.geometry_type
        for feature in Progress(self, "Clipping", quiet):
            geometry = feature.geometry
            if geometry is None:
                continue
            for mask_feature in mask.spatial_subset(geometry):
                clipped = clip_geometry(geometry, mask_feature.geometry, type_geometry)
                if clipped is None:
                    continue
                for clip in clipped:
                    target.add(geometry=clip, items=feature.items)

        return target

    def difference(self, mask, quiet: bool = True, ignore_errors: bool = False):
        """
        This function takes a difference between vector layer and difference layer.
        - Takes into account multiparts and single parts.
        - It also leaves geometries which are not valid.
        params:
            mask: Vector in which the difference is taken.
            quiet: Shows progress or not.
            ignore_errors: Ignores errors if True.

        """

        geometry_types = SINGLE_TO_MULTIPLE[self.geometry_type]
        target = self.copy(shell=True)
        for feature in Progress(self, "Difference", quiet):
            if feature is None:
                continue

            difference = feature.geometry
            masks = mask.spatial_subset(difference)

            for loop, mask_feature in enumerate(masks):
                mask_geometry = mask_feature.geometry
                difference = difference_geometry(
                    difference, mask_geometry, ignore_errors
                )

                # print if multiplication of 1000 features
                if loop % 5000 == 0 and loop > 0:
                    print(f"Doing difference at {loop} mask features")

            diff_part_type = difference.GetGeometryType()
            if diff_part_type not in geometry_types:
                if diff_part_type == ogr.wkbGeometryCollection:
                    for part in difference:
                        if part.GetGeometryType() == self.geometry_type:
                            target.add(geometry=part, items=feature.items)
            else:
                target.add(geometry=difference, items=feature.items)
        return target

    def spatial_clusters(self, quiet=True) -> (list, list):
        """
        Returns a list of ids of geometries that intersect each other.
        Returns a list of multipolygons of clustered geometries.
        """
        # TODO does not work yet.
        work = self.copy(shell=False)
        target = self.copy(shell=True)
        target.add_field("dids", str)

        total = len(work)
        work.reset()

        progress = Progress(message="Spatial clustering", total=total, quiet=quiet)

        # keep repeating until the work layer is empty
        while len(work) != 0:

            # new cluster id list
            cluster_ids = []
            cluster_geometry = ogr.Geometry(ogr.wkbMultiPolygon)

            # get first feature, geometry, first spatial filter, delete from work
            feature = next(work)
            if feature is None:
                continue

            cluster_ids.append(feature.id)
            work.delete(feature.id)
            cluster_geometry.AddGeometry(feature.geometry)
            cluster_geometry = ogr.ForceToMultiPolygon(cluster_geometry.UnionCascaded())

            # keep repeating until no more growth of the cluster
            while True:
                growth = False  # we are not growing
                for intersect in work.spatial_subset(cluster_geometry):
                    intersect_geometry = intersect.geometry
                    # we are growing
                    growth = True
                    # add this to the cluster
                    cluster_ids.append(intersect.id)
                    cluster_geometry.AddGeometry(intersect_geometry)
                    cluster_geometry = ogr.ForceToMultiPolygon(
                        cluster_geometry.UnionCascaded()
                    )

                    work.delete(intersect.id)

                if not growth:  # no more growth
                    break

            # add the cluster to the cluster_list

            dids = "_".join([str(c) for c in cluster_ids])
            target.add(geometry=cluster_geometry, dids=dids)
            progress.show(custom_add=len(cluster_ids))

        return target

    def dissolve(self, field: str = None, quiet: bool = True):
        """
        Dissolves the vector into a single multipolygon feature.
        params:
            field: The vector will be dissolved on the fields' unique values.
            quiet: verbosity.
        """

        multi_type = MULTI_TYPES[self.geometry_type]
        target = self.copy(shell=True, geometry_type=multi_type)
        if field is not None:  # finds unique fields, then dissolves

            unique = {f: [] for f in set(self.table[field])}

            for feature in self:
                unique[feature[field]].append(feature.id)

            for field_value, fid_list in Progress(
                unique.items(), f"Dissolve field {field}", quiet
            ):
                multi = ogr.Geometry(multi_type)
                for fid in fid_list:
                    multi.AddGeometry(self[fid].geometry)

                target.add(geometry=dissolve_geometry(multi))

        else:
            # here we are just dissolving all
            singles = self.to_single()
            multi = ogr.Geometry(multi_type)
            for feature in singles:
                multi.AddGeometry(feature.geometry)

            logger.info("Dissolving all")
            target.add(geometry=dissolve_geometry(multi))

        return target

    def voronoi(self, buffer=1, quiet=True):
        """
        Create vonoroi polygons based on a point layer.
        params:
            Buffer: used to include outer points.
            quiet: verbosity.

        """
        target = self.copy(shell=True, geometry_type=ogr.wkbPolygon)
        extent_geometry = self.extent_geometry

        point_list = []
        point_items = {}
        for fid, feature in Progress(
            self, "Vonoroi - Setting input points", enum=True, quiet=quiet
        ):
            point_list.append([*feature.geometry.points[0]])
            point_items[fid] = feature.items

        # add extent to properly add edge points
        for line in extent_geometry:
            for point in line.GetPoints():
                point_list.append([point[0], point[1]])

        vor = spatial.Voronoi(point_list)
        for fid, region_index in Progress(
            vor.point_region,
            "Vonoroi - Processing output polygons",
            enum=True,
            quiet=quiet,
        ):
            rpoints = []
            for vertex_index in vor.regions[region_index]:
                if vertex_index != -1:  # the library uses this for infinity
                    rpoints.append(list(vor.vertices[vertex_index]))
            rpoints.append(rpoints[0])
            if len(rpoints) < 4 or fid > len(point_items) - 1:  # anders geen polygon
                continue

            polygon = Polygon.from_points(rpoints)
            if not polygon.Within(extent_geometry):
                polygon = polygon.Intersection(extent_geometry)

            target.add(geometry=polygon, items=point_items[fid])
        return target

    def rasterize(
        self,
        rows=None,
        columns=None,
        geotransform=None,
        resolution=None,
        nodata=-9999,
        field=None,
        extent=None,
        all_touches=False,
        options=None,
        return_ds=True,
        data_type=gdal.GDT_Float32,
    ):
        """Rasterizes the vector as a boolean,
        If field or all touches is given that is used
        """

        if not columns and not rows:
            if not resolution:
                logger.error("please provide resolution")
                return
            if not extent:
                extent = self.extent
            x1, x2, y1, y2 = extent
            columns = int((x2 - x1) / resolution)
            rows = int((y2 - y1) / resolution)
            geotransform = (x1, resolution, 0, y2, 0, -resolution)

        if data_type == float:
            data_type = gdal.GDT_Float32
        elif data_type == int:
            data_type = gdal.GDT_Int32

        return rasterize(
            self.layer,
            rows,
            columns,
            geotransform,
            self.spatial_reference.wkt,
            nodata,
            field,
            all_touches,
            options,
            return_ds,
            data_type,
        )

    def interpolate(
        self,
        field,
        resolution,
        algorithm="invdist",
        nodata_value=-9999,
        quiet=True,
        **arguments,
    ):
        """
        Returns a gdal datasource.  Can be loaded with Raster.

        Interpolates vector points using gdal grid
        https://gdal.org/programs/gdal_grid.html#interpolation-algorithms


        params:
            field: Fields of which is used to interpolate
            resolution: Resulution of the final raster
            algorithm: See the above link to see how it can be filled.
            argument: See the above link to how it can be filled.

        """
        logger.debug("vector - interpolate - Start")
        progress = Progress(total=100, message="GDALGrid Interpolate", quiet=quiet)
        argument = ""
        for k, v in arguments.items():
            argument = argument + f"{k}={v}:"

        if len(argument) == 0:
            raise ValueError(
                "Please fill in an arguments, look at (https://gdal.org/programs/gdal_grid.html#interpolation-algorithms)"
            )

        assert field in self.fields

        extent = self.extent
        width = int((extent[1] - extent[0]) / resolution)
        height = int((extent[3] - extent[2]) / resolution)
        grid_options = gdal.GridOptions(
            zfield=field,
            algorithm=f"{algorithm}:{argument}",
            width=width,
            height=height,
            callback=progress.gdal,
        )

        path = mem_path()
        ds = gdal.Grid(
            destName=path,
            srcDS=self.as_gdal_dataset,
            options=grid_options,
        )

        # The raster which is extracted from the gdalGrid is acutally upside
        # down =with epsg 28992.
        # So we have to flip the array and change the geotransform.

        # We will now reverse it to the left upper corner.
        transform = list(ds.GetGeoTransform())
        transform[3] = transform[3] + (height * resolution)
        transform[-1] = -transform[-1]
        ds.SetGeoTransform(tuple(transform))

        band = ds.GetRasterBand(1)
        array = band.ReadAsArray()
        band.WriteArray(np.flip(array, 0), 0, 0)

        # set nodata value
        band.SetNoDataValue(-9999)

        band = None
        ds = None

        return gdal.Open(path)


def rasterize(
    layer,
    rows,
    columns,
    geotransform,
    spatial_reference_wkt,
    nodata=-9999,
    field=None,
    all_touches=False,
    options=None,
    return_ds=False,
    data_type=gdal.GDT_Float32,
):
    target_ds = DRIVER_GDAL_MEM.Create("", columns, rows, 1, data_type)

    # set nodata
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.Fill(nodata)
    band.FlushCache()

    # set metadata
    target_ds.SetProjection(spatial_reference_wkt)
    target_ds.SetGeoTransform(geotransform)

    # set options
    gdal_options = []

    if field:
        logger.info(field)
        gdal_options.append(f"ATTRIBUTE={field}")

    if all_touches:
        gdal_options.append("ALL_TOUCHES=TRUE")

    if options:
        gdal_options.extend(options)

    if len(gdal_options) == 0:
        gdal.RasterizeLayer(target_ds, (1,), layer, burn_values=(1,))
    else:
        gdal.RasterizeLayer(target_ds, [1], layer, options=gdal_options)

    if return_ds:
        return target_ds
    else:
        array = target_ds.ReadAsArray()
        target_ds = None
        return array
