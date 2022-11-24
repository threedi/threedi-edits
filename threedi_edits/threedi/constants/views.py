# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 09:50:02 2021

@author: chris.kerklaan

Note ogr uses spatialite 4.X instead 3, 
which requires a read-only (read_only ) parameter and changing values to lowercase

"""


def create_views(ogr_ds):
    """
    Change view geometry table to make read-only have a default
    And add views per layer
    """

    # sql = """DROP TABLE IF EXISTS views_geometry_columns;"""
    # ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    # sql = """CREATE TABLE views_geometry_columns (
    #         view_name TEXT NOT NULL,
    #         view_geometry TEXT NOT NULL,
    #         view_rowid TEXT NOT NULL,
    #         f_table_name TEXT NOT NULL,
    #         f_geometry_column TEXT NOT NULL,
    #         read_only INTEGER DEFAULT 1 NOT NULL,
    #         CONSTRAINT pk_geom_cols_views PRIMARY KEY (view_name, view_geometry),
    #         CONSTRAINT fk_views_geom_cols FOREIGN KEY (f_table_name, f_geometry_column) REFERENCES geometry_columns (f_table_name, f_geometry_column) ON DELETE CASCADE,
    #         CONSTRAINT ck_vw_rdonly CHECK (read_only IN (0,1)))"""

    # ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        CREATE VIEW IF NOT EXISTS v2_manhole_view
        AS SELECT manh.rowid AS ROWID, node.id AS node_id, manh.bottom_level
        AS manh_bottom_level, manh.surface_level AS manh_surface_level,
        manh.display_name AS manh_display_name, manh.shape AS manh_shape,
        manh.width AS manh_width, manh.length AS manh_length,
        manh.manhole_indicator AS manh_manhole_indicator, manh.calculation_type
        AS manh_calculation_type, manh.drain_level AS manh_drain_level,
        manh.zoom_category AS manh_zoom_category, node.initial_waterlevel AS
        node_initial_waterlevel, manh.id AS manh_id, manh.connection_node_id  AS
        manh_connection_node_id, node.storage_area AS node_storage_area,
        manh.code AS manh_code, node.code AS node_code, node.the_geom,
        node.the_geom_linestring AS node_the_geom_linestring,
        manh.sediment_level AS manh_sediment_level
        FROM v2_manhole manh, v2_connection_nodes node
        WHERE manh.connection_node_id = node.id;
        """
    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        DELETE FROM views_geometry_columns
        WHERE view_name = 'v2_manhole_view';"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        INSERT INTO views_geometry_columns (view_name, view_geometry,
        view_rowid, f_table_name, f_geometry_column)
        VALUES('v2_manhole_view', 'the_geom', 'rowid',
        'v2_connection_nodes', 'the_geom');"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        CREATE VIEW IF NOT EXISTS v2_pumpstation_point_view
        AS SELECT a.ROWID AS ROWID, a.id AS pump_id, a.display_name, a.code,
        a.classification, a.sewerage, a.start_level, a.lower_stop_level,
        a.upper_stop_level, a.capacity, a.zoom_category,
        a.connection_node_start_id, a.connection_node_end_id, a.type,
        b.id AS connection_node_id, b.storage_area, b.the_geom
        FROM v2_pumpstation a
        JOIN v2_connection_nodes b
        ON a.connection_node_start_id = b.id;"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        DELETE FROM views_geometry_columns
        WHERE view_name = 'v2_pumpstation_point_view';"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        INSERT INTO views_geometry_columns (view_name, view_geometry,
        view_rowid, f_table_name, f_geometry_column)
        VALUES('v2_pumpstation_point_view', 'the_geom',
        'connection_node_start_id', 'v2_connection_nodes', 'the_geom');"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        CREATE VIEW IF NOT EXISTS v2_1d_lateral_view
        AS SELECT a.ROWID AS ROWID, a.id AS id,
        a.connection_node_id AS connection_node_id,
        a.timeseries AS timeseries, b.the_geom
        FROM v2_1d_lateral a
        JOIN v2_connection_nodes b ON a.connection_node_id = b.id;"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        DELETE FROM views_geometry_columns
        WHERE view_name = 'v2_1d_lateral_view';"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        INSERT INTO views_geometry_columns (view_name, view_geometry,
        view_rowid, f_table_name, f_geometry_column)
        VALUES('v2_1d_lateral_view', 'the_geom', 'connection_node_id',
        'v2_connection_nodes', 'the_geom');"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        CREATE VIEW IF NOT EXISTS v2_1d_boundary_conditions_view
        AS SELECT a.ROWID AS ROWID, a.id AS id,
        a.connection_node_id AS connection_node_id,
        a.boundary_type AS boundary_type, a.timeseries AS timeseries,
        b.the_geom
        FROM v2_1d_boundary_conditions a
        JOIN v2_connection_nodes b ON a.connection_node_id = b.id;"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        DELETE FROM views_geometry_columns
        WHERE view_name = 'v2_1d_boundary_conditions_view';"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        INSERT INTO views_geometry_columns (view_name, view_geometry,
        view_rowid, f_table_name, f_geometry_column)
        VALUES('v2_1d_boundary_conditions_view', 'the_geom',
        'connection_node_id', 'v2_connection_nodes', 'the_geom');"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        CREATE VIEW IF NOT EXISTS v2_cross_section_location_view
        AS SELECT loc.ROWID as ROWID, loc.id as loc_id, loc.code as loc_code,
        loc.reference_level as loc_reference_level,
        loc.bank_level as loc_bank_level, loc.friction_type as
        loc_friction_type, loc.friction_value as loc_friction_value,
        loc.definition_id as loc_definition_id, loc.channel_id as
        loc_channel_id, loc.the_geom as the_geom, def.id as def_id,
        def.shape as def_shape, def.width as def_width, def.code as
        def_code, def.height as def_height
        FROM v2_cross_section_location loc, v2_cross_section_definition def
        WHERE loc.definition_id = def.id;"""
    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        DELETE FROM views_geometry_columns
        WHERE view_name = 'v2_cross_section_location_view';"""
    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        INSERT INTO views_geometry_columns (view_name, view_geometry,
        view_rowid, f_table_name, f_geometry_column)
        VALUES('v2_cross_section_location_view', 'the_geom',
        'rowid', 'v2_cross_section_location', 'the_geom');"""
    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """CREATE VIEW IF NOT EXISTS v2_pipe_view AS
            SELECT pipe.rowid AS rowid,
            pipe.id AS pipe_id,
            pipe.display_name AS pipe_display_name,
            pipe.code AS pipe_code,
            pipe.profile_num AS pipe_profile_num,
            pipe.sewerage_type AS pipe_sewerage_type,
            pipe.calculation_type AS pipe_calculation_type,
            pipe.invert_level_start_point AS pipe_invert_level_start_point,
            pipe.invert_level_end_point AS pipe_invert_level_end_point,
            pipe.cross_section_definition_id AS pipe_cross_section_definition_id,
            pipe.friction_value AS pipe_friction_value,
            pipe.friction_type AS pipe_friction_type,
            pipe.dist_calc_points AS pipe_dist_calc_points,
            pipe.material AS pipe_material,
            pipe.pipe_quality AS pipe_pipe_quality,
            pipe.original_length AS pipe_original_length,
            pipe.zoom_category AS pipe_zoom_category,
            pipe.connection_node_start_id AS pipe_connection_node_start_id,
            pipe.connection_node_end_id AS pipe_connection_node_end_id,
            def.id AS def_id,
            def.shape AS def_shape,
            def.width AS def_width,
            def.height AS def_height,
            def.code AS def_code,
            MakeLine(start_node.the_geom, end_node.the_geom) AS the_geom
            FROM v2_pipe AS pipe
            , v2_cross_section_definition AS def
            , v2_connection_nodes AS start_node
            , v2_connection_nodes AS end_node
            WHERE pipe.connection_node_start_id = start_node.id
            AND pipe.connection_node_end_id = end_node.id
            AND pipe.cross_section_definition_id = def.id"""
    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        DELETE FROM views_geometry_columns
        WHERE view_name = 'v2_pipe_view';"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    sql = """
        INSERT INTO views_geometry_columns (view_name, view_geometry,
        view_rowid, f_table_name, f_geometry_column)
        VALUES('v2_pipe_view', 'the_geom','rowid',
               'v2_connection_nodes', 'the_geom_linestring');"""

    ogr_ds.ExecuteSQL(sql, dialect="SQLite")

    return ogr_ds
