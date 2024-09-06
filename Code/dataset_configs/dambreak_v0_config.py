# Specify here the folder (location) in which the GIS folder is located
folder_path_GIS = r"D:\Work\Project\P1414"

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
