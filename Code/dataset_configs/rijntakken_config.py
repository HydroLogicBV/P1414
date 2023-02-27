class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 0.1
            node_distance = 100

        class two_d:
            coupling_type = "1Dto2D"
            dx = 500
            dy = 500
            elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            two_d_buffer = 100


class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    branches_path = p_folder + r"\Rijntakken\without_waal\RTK_Branches.shp"
    river_profile_path = p_folder + r"\Rijntakken\without_waal\ZW_cross.csv"
    river_roughness_path = p_folder + r"\Rijntakken\without_waal\roughness.csv"
    weir_path = p_folder + r"\Rijntakken\without_waal\RTK_Weirs.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    class Peil:
        default_peil = 4

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", None),
            ("bodemhoogte benedenstrooms", None),
            ("bodemhoogte bovenstrooms", None),
            ("code", "Name"),
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
        ]
    )
