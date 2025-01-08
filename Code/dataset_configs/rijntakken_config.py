import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Rijntakken'

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 0.1
            node_distance = 100

        class two_d:
            coupling_type = "1Dto2D"
            dx = 500
            dy = 500
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            two_d_buffer = 100


class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    branches_path = p_folder + r"\Rijntakken\RTK_BranchesV3.shp"
    river_profile_path = p_folder + r"\Rijntakken\ZW_cross_v3.csv"
    river_roughness_path = p_folder + r"\Rijntakken\roughness_v2.csv"
    weir_path = p_folder + r"\Rijntakken\RTK_Weirs.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    class Peil:
        default_peil = 4

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "Name"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("is_duiker", None),
            ("typeruwheid", None),
        ]
    )

    ## Weirs
    weir_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "Name"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "CrestWidth"),
            ("hoogstedoorstroomhoogte", "CrestLevel"),
            ("laagstedoorstroombreedte", "CrestWidth"),
            ("laagstedoorstroomhoogte", "CrestLevel"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
            ("flowdir", None),
        ]
    )
