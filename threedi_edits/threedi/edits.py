# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:43:02 2020

@author: chris.kerklaan

ThreediEdits is a class based on vector.py and raster.py to edit threedi models
The class always edits the model in memory and you'll have to write the model
to show the results
"""

# System imports
import os
import logging
import pathlib

# Third-party imports
from osgeo import ogr
from cached_property import cached_property

# Local imports
from threedi_edits.gis.vector import Vector
from threedi_edits.threedi.grid import make_grid
from threedi_edits.threedi.rastergroup import ThreediRasterGroup
from threedi_edits.threedi.vectorgroup import ThreediVectorGroup
from threedi_edits.utils.dependencies import DEPENDENCIES
from threedi_edits.threedi.constants.constants import TABLES, raster_fields, TRANSLATE
from threedi_edits.utils.project import Classes, Functions


# structure
classes = Classes(__name__, local_only=True)
functions = Functions(__name__)

# logger
logger = logging.getLogger(__name__)


# Globals
OGR_SQLITE_DRIVER = ogr.GetDriverByName("SQLite")
OGR_MEM_DRIVER = ogr.GetDriverByName("Memory")
SUPPORTED_VERSIONS = ["0209"]


class ThreediEdits(ThreediVectorGroup):
    """
    An object for editing a threedi model,
    Can be openen from scratch, from any ogr format (sqlite, postgres, geopackage).

        mode:
            'write': Used to write a file
            'empty': Returns an empty in-mem ogr threedi model
            'memory': Used when a full threedimodel is presented in memory
            'read': Used for reading a model only

    Example code:
        sqlite_path = "path/to/schema.sqlite"
        schema = ThreediEdits(sqlite_path)

        # retrieving tables
        nodes = schema['v2_connection_nodes']
        nodes = schema.nodes

        # resetting tables
        new_nodes = nodes.copy()
        schema['nodes'] = new_nodes
        schema.nodes = new_nodes

        # retrieving rasters
        all_rasters = schema.rasters
        dem = schema.dem

        # resetting rasters
        new_dem = dem.copy()
        schema.dem = dem

        # node info
        schema.nodes_view
        schema.nodes_height


    """

    instance = "threedi.edits.ThreediEdits"

    def __init__(
        self,
        sqlite_path: str = None,
        mode: str = "read",
        datasource: ogr.DataSource = None,
    ):

        if sqlite_path is not None and not os.path.exists(str(sqlite_path)):
            raise FileNotFoundError("Path does not exist")

        if sqlite_path:
            datasource = ogr.Open(str(sqlite_path), 0)
            self.model_dir = os.path.dirname(sqlite_path) + "/"
            self.name = pathlib.Path(sqlite_path).stem
            self.path = os.path.join(os.getcwd(), sqlite_path)

        if mode == "empty":
            datasource = None

        super().__init__(datasource, mode)

        logger.debug("Initializing tables")
        self.tables
        self.supported(self.version)

    @classmethod
    def from_pg(
        cls,
        dbname,
        host="nens-3di-db-03.nens.local",
        port="5432",
        user="threedi",
        password="1S418lTYWLFsYxud4don",
    ):

        pg_str = ("PG:host={} port={} user='{}'" "password='{}' dbname='{}'").format(
            host, port, user, password, dbname
        )
        return cls(pg_str=ogr.Open(pg_str, 0))

    @classmethod
    def from_scratch(cls):
        return cls(mode="empty")

    def copy(self):
        return ThreediEdits(datasource=self._copy(), mode="read")

    def grid(self):
        if DEPENDENCIES.threedigrid_builder.installed:
            return make_grid(self)
        else:
            raise ImportError("Please install threedigrid builder!")

    @property
    def empty(self) -> bool:
        """
        Returns true or false based on the first global setting row
        """
        layer = self.ds.GetLayerByName("v2_global_settings")
        return layer.GetFeatureCount() == 0

    @cached_property
    def rasters(self):
        rasters = ThreediRasterGroup(**self._existing_raster_paths)
        for raster in rasters:
            setattr(self, raster.name, raster)
        return rasters

    @cached_property
    def tables(self):
        """
        The 3Di tables are cached so they can be edited in place.
        The reference stays the same, no matter the place.
        We should also be able to replace the cache.

        """

        tables = self._tables()
        for table in TABLES["all"]:
            if table in TRANSLATE:
                attribute = getattr(tables, TRANSLATE[table])
                setattr(self, TRANSLATE[table], attribute)

        return tables

    @property
    def _get_raster_files(self) -> dict:
        """Retrieves all raster files from the global settings,
        Returns empty dict if no global settings.

        """

        if self.empty:
            return {}

        global_settings = self.ds.GetLayerByName("v2_global_settings")
        gs = next(iter(global_settings))
        paths = {field: gs[field] for field in raster_fields if field in gs.keys()}
        infiltration_settings = self.ds.GetLayerByName("v2_simple_infiltration")
        infiltration_id = gs["simple_infiltration_settings_id"]
        if type(infiltration_id) == int:
            infiltration = infiltration_settings.GetFeature(infiltration_id)
            if infiltration is not None:
                paths["infiltration_rate_file"] = infiltration["infiltration_rate_file"]
                paths["max_infiltration_capacity_file"] = infiltration[
                    "max_infiltration_capacity_file"
                ]
        return paths

    @property
    def _existing_raster_paths(self) -> dict:
        """
        Retrieves raster files if they exist.
        Return an empty dictionary otherwise.
        """
        if not hasattr(self, "model_dir"):
            return {}
        else:
            paths = {}
            for k, v in self._get_raster_files.items():
                if v is None:
                    continue
                if v == "":
                    continue

                if os.path.exists(self.model_dir + v):
                    paths[k] = self.model_dir + v
            return paths

    @property
    def files(self):
        return self._get_raster_files

    @property
    def placed_rasters(self):
        return self.rasters.count > 0

    @property
    def has_existing_rasters(self):
        return len(self._existing_raster_paths) > 0

    @property
    def version(self):
        versioning = self._table("schema_version")
        if len(versioning) > 0:
            return versioning[1]["version_num"]
        else:
            logger.info("Version not set!")
            return SUPPORTED_VERSIONS[0]

    @cached_property
    def extent(self):
        extent = Vector.from_scratch("", 3, self.rasters.dem.epsg)
        extent.add(geometry=self.rasters.dem.extent_geometry)
        return extent

    @cached_property
    def shape(self):
        raster = self.dem.copy()
        raster.array = self.dem.mask
        return Vector(raster.polygonize(quiet=False, mask_value=0))

    @cached_property
    def extent_indices(self):
        return self.rasters.dem.extent

    @cached_property
    def extent_geometry(self):
        return self.rasters.dem.extent_geometry

    def write(self, path, rasters=False, rasters_correct=True, quiet=False):
        path = str(path)

        if (self.has_existing_rasters or self.placed_rasters) and rasters:
            group = self.rasters
            if rasters_correct:
                group.correct()
            group.write(str(pathlib.Path(path).parent), self._get_raster_files)

        self._write(path, quiet)

    @property
    def epsg(self):
        return self.global_setting["epsg"]

    @cached_property
    def global_setting(self):
        return next(iter(self["v2_global_settings"]))

    # TODO what to do with the data.
    @cached_property
    def nodes_view(self, data=False):
        node_table = self.nodes.table
        new_dict = {fid: {} for fid in node_table["fid"]}

        # add connection node info
        for i in range(0, self.nodes.count):
            new_dict[node_table["fid"][i]]["v2_connection_nodes"] = {
                key: node_table[key][i] for key in node_table
            }

        for table in TABLES["single"] + TABLES["start_end"]:
            model_table = self[table]

            add = model_table.table
            for i in range(len(self[table])):

                if table in TABLES["single"]:
                    node_table = new_dict[add["connection_node_id"][i]]

                    if data:
                        node_table[table] = {key: add[key][i] for key in add}
                    else:
                        if table not in node_table:
                            node_table[table] = []

                        node_table[table].append(add["fid"][i])
                else:

                    node_table_start = new_dict[add["connection_node_start_id"][i]]
                    node_table_end = new_dict[add["connection_node_end_id"][i]]

                    if data:
                        if add["connection_node_start_id"][i]:
                            node_table_start[table + "_start"] = {
                                key: add[key][i] for key in add
                            }

                        if add["connection_node_end_id"][i]:
                            node_table_end[table + "_end"] = {
                                key: add[key][i] for key in add
                            }
                    else:

                        if table not in node_table_start:
                            node_table_start[table] = []

                        if table not in node_table_end:
                            node_table_end[table] = []

                        node_table_start[table].append(add["fid"][i])
                        node_table_end[table].append(add["fid"][i])

        return new_dict

    def set_rasters(self, value):
        self.__dict__["rasters"] = value

    def supported(self, version):
        if not version in SUPPORTED_VERSIONS:
            raise NotImplementedError(f"Version {version} not implemented")

    def nodes_height(self):
        """Samples all nodes using the dem."""
        dem = self.dem
        dem.return_raster = True
        project = dem.reproject(4326)
        return {node.fid: project.read(node.geometry) for node in self.nodes}

    def spatial_subset(self, geometry: ogr.Geometry):
        """
        Filter the model on connection nodes. Shapely or ogr geometry
        Also retrieves all data
        """
        nodes = self.nodes.spatial_subset(geometry)
        filtered = {}
        for table in self:
            filtered[table.name] = []
            for node in nodes:
                filtered[table.name].append(self.nodes_view[node.id])
        return filtered

    def delete_tables(self, deletes=None, quiet=True):
        """
        Deletes features in tabels in the correct order
        Parameters
        ----------
        db : ThreediEdits
        deletes : dict
            dict with the table name and fids to delete e.g., {'v2_pipe': [1,2,3,4]}.
        quiet : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """
        for table in TABLES["order"]:
            if deletes:
                if table not in deletes:
                    continue

                fids = deletes[table]
            else:
                fids = self[table].fids

            for fid in fids:
                self[table].delete(fid)

    def delete_node(self, node_id, clear=True):
        """
        Deletes a node and all its assciates.
        """

        view = self.nodes_view[node_id]
        deletes = {k: [] for k in TABLES["order"]}
        for delete_table in deletes:
            if delete_table in view:
                if type(view[delete_table]) == dict:
                    fid = view[delete_table]["fid"]
                    if not fid in deletes[delete_table]:
                        deletes[delete_table].append(fid)
                else:
                    for fid in view[delete_table]:
                        if not fid in deletes[delete_table]:
                            deletes[delete_table].append(fid)
        # associates = [(i, deletes[i]) for i in deletes if len(deletes[i]) > 0]
        logger.info("Deleting associates f{associates}")

        if clear:
            self.delete_tables(deletes, quiet=True)
        else:
            return deletes
