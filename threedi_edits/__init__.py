"""  
threedi_edits is a python toolbox for gis and threedi related processing.
Note that the only depedency for base processing should be gdal.

Release info for gdal:
https://github.com/OSGeo/gdal/blob/v3.4.1/gdal/NEWS.md

General imports from threedi_edits are listed below:
    
E.g., from threedi_edits import Raster

-- Use help(Raster) for further information
-- Use .modules for all possible imports of modules
-- Use .classes for all possible imports of classes
-- Use .functions for all possible imports of functions


gis-processing:
    - Raster
    - RasterGroup
    - Vector
    - VectorGroup
    - LineString
    - MultiLineString
    - Polygon
    - MultiPolygon
    - Point 
    - MultiPoint
    
threedi-processing
    - ThreediRasterGroup (used for the rasters of threedi)
    - ThreediEdits (used for the model database (.sqlite) of threedi)
            
Philosophy:
    We want to create a simple api on a threedi database.
    Using as few dependencies as possible.
    
About version 1.0:
    In version 1.0 we changed vector processing is changed by using shapely.
    
    
"""
from .utils.logging import show_console_logging

show_console_logging()

from .utils.dependencies import DEPENDENCIES

DEPENDENCIES.missing()


from osgeo import gdal, ogr

# options
gdal.SetConfigOption("SQLITE_LIST_ALL_TABLES", "YES")
gdal.SetConfigOption("PG_LIST_ALL_TABLES", "YES")

# Exceptions
gdal.SetConfigOption("CPL_LOG", "NUl")
ogr.UseExceptions()
gdal.UseExceptions()

# versioning
__version__ = "1.0"

# gis
import shapely  # noqa
from . import gis  # noqa
from .gis.raster import Raster  # noqa
from .gis.rastergroup import RasterGroup  # noqa
from .gis.vector import Vector  # noqa
from .gis.vectorbase import VectorBase  # noqa
from .gis.vectorgroup import VectorGroup  # noqa
from .gis.linestring import LineString, MultiLineString  # noqa
from .gis.polygon import Polygon, MultiPolygon  # noqa
from .gis.point import Point, MultiPoint  # noqa

# Threedi
from . import threedi  # noqa
from .threedi.rastergroup import ThreediRasterGroup  # noqa
from .threedi.rastergroup import retrieve_soil_conversion_table  # noqa
from .threedi.rastergroup import retrieve_landuse_conversion_table  # noqa
from .threedi.edits import ThreediEdits  # noqa
from .threedi.constants.templates import Templates  # noqa

# Logging
from . import utils  # noqa
from .utils.progress import Progress  # noqa
from .utils.project import log_time  # noqa
from .utils.project import Files  # noqa
from .utils.project import Functions  # noqa
from .utils.project import Classes  # noqa
from .utils.project import Modules  # noqa
from .utils.project import Logger  # noqa

# Examples
from .examples import get_examples  # noqa

# structure
files = Files(__file__)
functions = Functions(__name__)
classes = Classes(__name__, local_only=False)
modules = Modules(__name__)
