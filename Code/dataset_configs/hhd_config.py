class Models:
    class FM:
        dx = 100
        dy = 100
        elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
        max_snap_dist = 1
        one_d = True
        start_time = 20160601
        stop_time = 2 * 86400
        two_d = True
        two_d_buffer = 100


class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    branches_path = p_folder + r"\Uitgesneden watergangen\HHD_v4_test.shp"
    # bridges_path = p_folder + r"\HDSR\Legger\Bruggen\Bruggen.shp"
    # culvert_path = p_folder + r"\HDSR\Legger\Kokers_Lijnen\Kokers_Lijnen_edited.shp"
    pump_path = p_folder + r"\HHDelfland\Niet legger\Gemaal_peil.shp"
    weir_path = p_folder + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Stuw.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", None),
            ("bodemhoogte benedenstrooms", None),
            ("bodemhoogte bovenstrooms", None),
            ("code", "CODE"),
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

    # ## Bridges
    # bridge_index_mapping = dict(
    #     [
    #         ("code", "CODE"),
    #         ("geometry", "geometry"),
    #         ("globalid", "globalid"),
    #         ("intreeverlies", None),
    #         ("typeruwheid", None),
    #         ("ruwheid", None),
    #         ("uittreeverlies", None),
    #         ("lengte", "WS_DOORVAA"),
    #     ]
    # )

    ## Culverts
    # culvert_index_mapping = dict(
    #     [
    #         ("breedteopening", "BREEDTEOPE"),
    #         ("code", "CODE"),
    #         ("geometry", "geometry"),
    #         ("gesloten", None),
    #         ("globalid", "globalid"),
    #         ("hoogtebinnenonderkantbene", "HOOGTEBOKB"),
    #         ("hoogtebinnenonderkantbov", "HOOGTEBO_1"),
    #         ("hoogteopening", "HOOGTEOPEN"),
    #         ("intreeverlies", None),
    #         ("lengte", "LENGTE"),
    #         ("typeruwheid", None),
    #         ("ruwheid", None),
    #         ("uittreeverlies", None),
    #         ("vormkoker", "VORMKOKER"),
    #     ]
    # )

    ## Pumps
    pump_index_mapping = dict(
        [
            ("code", "CODE"),
            ("doelvariabele", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("maximalecapaciteit", "MAXCAPACIT"),
            ("streefwaarde", "streefpeil"),
            ("peil_marge", None),
        ]
    )

    ## Weirs
    weir_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "DOORSTROOM"),
            ("hoogstedoorstroomhoogte", "HOOGSTEDOO"),
            ("laagstedoorstroombreedte", "DOORSTROOM"),
            ("laagstedoorstroomhoogte", "LAAGSTEDOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", "REGELBAARH"),
            ("soortstuw", None),
            ("vormopening", None),
        ]
    )
