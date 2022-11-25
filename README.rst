threedi-edits
==========================================

Introduction

An experimental pythonic 3Di schematisation api. Using this api, we can access, alter en write a 3Di database within python. Within the package gis tools are provided as well.  
Currently version '0209' of the schematisation is supported. If your 3Di database is not up to date, please use: 'threedi_modelchecker (0.30) --sqlite "path/to/sqlite.sqlite" migrate' to migrate your database to version '0209'.

This package is the continuation of threedi-raster-edits (0.27)

Usage - 3Di schematisation edits
------------


Importing the module::

  >>> import threedi_edits as te
  
  
Creating a new schematisation::

	>>> schema = te.ThreediEdits.from_scratch()
	
Basic schematisation reading::

  >>> schema = te.ThreediEdits(sqlite_path, mode="read")
  >>> schema = te.ThreediEdits(sqlite_path, mode="memory")
  
Copying a schematisation::

	>>> copied = schema.copy()
	
Writing a schematisation::

	>>> schema.write(sqlite_output_path, rasters=True)


Investigate the data::

  >>> cross_section_locations = schema.cross_section_locations
  >>> channels = schema.channels
  >>> print(cross_section_locations.first()) # show data of the first cross section location
  >>> print(schema.tables) # show all tables
  >>> shapely_geom = schema.channels.first().geometry.shape # retrieve the shapely geometry
 
Creating the grid::

  >>> grid = schema.grid()
  >>> cells = grid['cells']
  >>> nodes = grid['nodes']
  >>> lines = grid['lines']
  >>> cells.write("path_to_lines.gpkg") # write cells only
  >>> grid.write("path_to_grid.gpkg") # Write lines, cells and nodes
 
 
Adding a feature using a template::

  >>> from te.globals import SUPPORTED_THREEDI_VERSIONS
  >>> from te.threedi.constants import get_version
  >>> templates = get_version(SUPPORTED_THREEDI_VERSIONS[0]).Templates()
  >>> node = templates.node
  >>> print(node) # show what you should fill
  >>> node["initial_waterlevel"] = 1
  >>> node["the_geom"] = Point.from_point((1, 1))
  >>> new_fid = schema.nodes.add(node)
  
Adding a feature using items and geometry::

  >>> items = {"initial_waterlevel": 0.1, "storage_are": 0.64, "code": "x"}
  >>> geometry = te.Point.from_point((1, 1))
  >>> new_fid = schema.nodes.add(items=items, geometry=geometry)

Adding a feature using another feature::

  >>> feature = other_schema.nodes.first()
  >>> new_fid = schema.add(feature)

Editing and writing rasters::

  >>> schema.rasters  # initialize rasters
  >>> dem = schema.dem
  >>> new_dem = dem.copy()
  >>> array = new_dem.array
  >>> array[0, 0] = 2
  >>> new_dem.array = array
  >>> schema.dem = new_dem
  
  
For more examples:
	te.get_examples("example_folder")