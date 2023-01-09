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
    branches_path = p_folder + r"\Uitgesneden watergangen\HDSR_v4_test.shp"
    bridges_path = p_folder + r"\HDSR\Legger\Bruggen\Bruggen.shp"
    culvert_path = p_folder + r"\HDSR\Legger\Kokers_Lijnen\Kokers_Lijnen_edited.shp"
    pump_path = p_folder + r"\HDSR\Legger\Gemalen\Gemalen_peil.shp"
    weir_path = p_folder + r"\HDSR\Legger\Stuwen\BR_Stuwen.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("hoogte insteek linkerzijde", "IWS_W_INST"),
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
            ("taludhelling linkerzijde", "IWS_W_TALU"),
            ("taludhelling rechterzijde", "IWS_W_TA_1"),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", "IWS_W_WATB"),
        ]
    )

    ## Bridges
    bridge_index_mapping = dict(
        [
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("intreeverlies", None),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("lengte", "WS_DOORVAA"),
        ]
    )

    ## Culverts
    culvert_index_mapping = dict(
        [
            ("breedteopening", "BREEDTEOPE"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "HOOGTEBOKB"),
            ("hoogtebinnenonderkantbov", "HOOGTEBO_1"),
            ("hoogteopening", "HOOGTEOPEN"),
            ("intreeverlies", None),
            ("lengte", "LENGTE"),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("vormkoker", "VORMKOKER"),
        ]
    )

    ## Pumps
    pump_index_mapping = dict(
        [
            ("code", "GEMAALID"),
            ("doelvariabele", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("maximalecapaciteit", "MAXIMALECA"),
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
            ("hoogstedoorstroombreedte", None),
            ("hoogstedoorstroomhoogte", "HOOGSTEDOO"),
            ("laagstedoorstroombreedte", "DOORSTROOM"),
            ("laagstedoorstroomhoogte", "LAAGSTEDOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", "SOORTREGEL"),
            ("soortstuw", "SOORTSTUW"),
            ("vormopening", None),
        ]
    )
