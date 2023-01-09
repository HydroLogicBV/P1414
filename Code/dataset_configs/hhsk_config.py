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
    branches_path = p_folder + r"\Uitgesneden watergangen\HHSK_v4_test.shp"
    culvert_path = p_folder + r"\HHSK\Legger\Duiker.shp"
    # pump_path = p_folder + r"\HHRijnland\Legger\Gemalen_peil\gemaal_peil.shp"
    weir_path = p_folder + r"\HHSK\Legger\Stuw.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", "BODEMBREED"),
            ("bodemhoogte benedenstrooms", None),
            ("bodemhoogte bovenstrooms", None),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("taludhelling linkerzijde", "TALUDLINKS"),
            ("taludhelling rechterzijde", "TALUDRECHT"),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", "WATERBREED"),
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
    #         ("lengte", "DOORSTRO_1"),
    #     ]
    # )

    ## Culverts
    culvert_index_mapping = dict(
        [
            ("breedteopening", "BREEDTE"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "HOOGTEBOK"),
            ("hoogtebinnenonderkantbov", "HOOGTEBOK"),
            ("hoogteopening", "HOOGTE"),
            ("intreeverlies", None),
            ("lengte", "LENGTE"),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("vormkoker", "VORM"),
        ]
    )

    ## Pumps
    pump_index_mapping = dict(
        [
            ("code", "CODE"),
            ("doelvariabele", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("maximalecapaciteit", "CAPACITEIT"),
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
            ("hoogstedoorstroombreedte", "KRUINBREED"),
            ("hoogstedoorstroomhoogte", "KRUINHOOGT"),
            ("laagstedoorstroombreedte", "KRUINBREED"),
            ("laagstedoorstroomhoogte", "KRUINHOOGT"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
        ]
    )
