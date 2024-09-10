# Specify here the folder (location) in which the GIS folder is located
folder_path_GIS = r"D:\Work\Project\P1414"

class FixedWeirs:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    flood_defences_path = p_folder + r"\Keringen_met_hoogte\Noordzeekeringen.shp"
    fixed_weir_index_mapping = dict(
        [
            ("code", "id"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )
