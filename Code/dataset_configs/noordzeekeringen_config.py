class FixedWeirs:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    flood_defences_path = p_folder + r"\Keringen_met_hoogte\Noordzeekeringen.shp"
    fixed_weir_index_mapping = dict(
        [
            ("code", "id"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )
