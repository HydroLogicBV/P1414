import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Tunnel'
    
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
            coupling_type = "1Dto2D"
            dx = 500
            dy = 500
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            initial_peil_raster_path = folder_path_GIS + r"\GIS\peilen\peilen_jp_25m_full.tif"
            two_d_buffer = 100


class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    branches_path = p_folder + r"\Tunnels en onderdoorgangen\tunnel_verdiepteA4V2.shp"
    culvert_path = p_folder + r"\Tunnels en onderdoorgangen\tunnel.shp"
    norm_profile_path = p_folder + r"\Tunnels en onderdoorgangen\tunnel_verdiepteA4.shp"

    class Peil:
        default_peil = -9.75

    

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "lokaalid"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", True),
            ("is_duiker", None),
            ("typeruwheid", None),
        ]
    )

    culvert_index_mapping = dict(
        [
            ("breedteopening", "width"),
            ("code", "lokaalid"),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "z_min"),
            ("hoogtebinnenonderkantbov", "z_min"),
            ("hoogteopening", "height"),
            ("intreeverlies", None),
            ("lengte", None),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("vormkoker", None),
        ]
    )

    ## Norm Profiles
    np_index_mapping = dict(
        [
            ("bodembreedte", "width"),
            ("bodemhoogte benedenstrooms", "z_min"),
            ("bodemhoogte bovenstrooms", "z_min"),
            ("code", "lokaalid"),
            ("diepte", "height"),
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
