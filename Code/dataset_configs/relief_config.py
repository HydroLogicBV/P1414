class FixedWeirs:
    ## PATHS
    dm_attribute = "overiglijnelement"
    p_folder = r"D:\Work\Project\P1414\GIS"
    flood_defences_path = p_folder + r"\Keringen_met_hoogte\relief_v2.shp"
    fixed_weir_index_mapping = dict(
        [
            ("code", "lokaalid"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )
