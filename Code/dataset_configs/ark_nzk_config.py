import os
# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Amsterdam-Rijnkanaal & Noordzeekanaal'
    
class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 5
            node_distance = 500

        class two_d:
            coupling_type = "1Dto2D"
            dx = 500
            dy = 500
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            two_d_buffer = 100

### PEIL ++
class FixedWeirs:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    

    flood_defences_path = dict(
                [
                     ("base", p_folder + r"\Keringen_met_hoogte\ARKV2.shp"),
                     ("concat", p_folder + r"\Keringen_met_hoogte\ExtraKeringenNZK.shp"),
                ]
            )

    fixed_weir_index_mapping = dict(
        [
            ("code", "objectid"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )


class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    
    branches_path = p_folder + r"\ARKNZK\netwerk_selectie_v7.shp"
    river_profile_path = p_folder + r"\ARKNZK\ZW_cross_v3.csv"
    sluice_path = p_folder + r"\ARKNZK\RWS_sluizen_clipped.shp"
    weir_path = p_folder + r"\ARKNZK\netwerk_watergangen_boezemmodel_WeirsV2.shp"

    class Peil:
        default_peil = -0.4

    
    weir_selection = dict([("column", "OPNEMEN"), ("value", "JA")])

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "Name"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("is_duiker",None),
            ("typeruwheid", None),
        ]
    )
    ## Sluice
    sluice_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "objectnaam"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", None),
            ("hoogstedoorstroomhoogte", None),
            ("laagstedoorstroombreedte", None),
            ("laagstedoorstroomhoogte", None),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
            ("flowdir", None),
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
    