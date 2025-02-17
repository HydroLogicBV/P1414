import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Onderdoorgangen'

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 0.1
            node_distance = 50

        class two_d:
            coupling_type = None
            dx = 500
            dy = 500
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            initial_peil_raster_path = folder_path_GIS + r"\GIS\peilen\peilen_jp_25m_full.tif"
            two_d_buffer = 100


class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    branches_path = dict(
        [
            ("base", p_folder + r"\Tunnels en onderdoorgangen\Underpasses_unfiltered_500m_len.shp"),
            ("concat", p_folder + r"\WAGV\AGV_Onderdoorgangen_Extra\AGV_Onderdoorgangen_Extra.shp"),
        ]
    )

    norm_profile_path = branches_path

    ahn_path = folder_path_GIS + r"\GIS\AHN\AHN4_WSS_filled.tif"
    keringen_path = folder_path_GIS + r"\GIS\HYDAMO\Combined_keringen.gpkg"
    tunnel_path = folder_path_GIS + r"\GIS\Keringen_met_hoogte\tunnel.shp"

    class Peil:
        default_peil = -9.75

    

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", ["lokaalid","id_agv"]),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", True),
            ("typeruwheid", None),
            ("is_duiker", None)
        ]
    )

    ## Norm Profiles
    np_index_mapping = dict(
        [
            ("bodembreedte", ["width", "Width"]),
            ("bodemhoogte benedenstrooms", ['z_min', 'z-height']),
            ("bodemhoogte bovenstrooms", ['z_min', 'z-height']),
            ("code", ["lokaalid", "id_agv"]),
            ("diepte", ["height", "Height"]),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("taludhelling linkerzijde", None),
            ("taludhelling rechterzijde", None),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", None),
        ]
    )
