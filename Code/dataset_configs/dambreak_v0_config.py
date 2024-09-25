import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Dambreak:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    dambreak_path = p_folder + r"\Dijkdoorbraken\Dijkdoorbraken_v0_lijn.shp"
    dambreak_index_mapping = dict(
        [
            ("algorithm", None),
            ("breachwidthini", "ini_breedt"),
            ("code", "id"),
            ("crestlevelini", "ini_diepte"),
            ("crestlevelmin", "min_hoogte"),
            ("f1", None),
            ("f2", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("t0", "t0"),
            ("timetobreachtomaximumdepth", None),
            ("ucrit", None),
        ]
    )
