# Specify here the folder (location) in which the GIS folder is located
folder_path_GIS = r"D:\Work\Project\P1414"

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 50


class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    # p_folder = r"D:\work\P1414_ROI\Boezemmodel_Waternet_dimr"
    branches_path = p_folder + r"\Noordzee\Noordzee.shp"
    norm_profile_path = p_folder + r"\Noordzee\Noordzee.shp"
    weir_path = p_folder + r"\Noordzee\stuw.shp"

    class Peil:
        default_peil = 4

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "Name"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("typeruwheid", None),
        ]
    )

    ## Branches
    np_index_mapping = dict(
        [
            ("bodembreedte", "Breedte"),
            ("bodemhoogte benedenstrooms", "H_NAP"),
            ("bodemhoogte bovenstrooms", "H_NAP"),
            ("code", "Name"),
            ("diepte", "Diepte"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("taludhelling linkerzijde", None),
            ("taludhelling rechterzijde", None),
            ("thalweg offset", "thal_offse"),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", None),
        ]
    )

    ## Weirs
    weir_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "id"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "BREEDTE"),
            ("hoogstedoorstroomhoogte", "HOOGTE"),
            ("laagstedoorstroombreedte", "BREEDTE"),
            ("laagstedoorstroomhoogte", "HOOGTE"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
        ]
    )
