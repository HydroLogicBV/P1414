import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class FixedWeirs:
    ## PATHS
    dm_attribute = "overiglijnelement"
    p_folder = folder_path_GIS + r"\GIS"
    flood_defences_path = p_folder + r"\Keringen_met_hoogte\hoofd_en_spoorwegen_v2.shp"
    fixed_weir_index_mapping = dict(
        [
            ("code", "lokaalid"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )
