class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20160601
        stop_time = 2 * 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 0.1
            node_distance = 100

        class two_d:
            coupling_type = "2Dto1D"
            dx = 100
            dy = 100
            elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            two_d_buffer = 100


class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    branches_path = p_folder + r"\WAGV\Niet legger\hydrovak_combined.gpkg"
    bridges_path = p_folder + r"\WAGV\brug_v13\brug_v13_clipped.shp"
    culvert_path = p_folder + r"\WAGV\duikersifonhevel_v13\duikersifonhevel_v13_clipped.shp"
    measured_profile_path = (
        p_folder + r"\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped_rm.shp"
    )
    pump_path = p_folder + r"\WAGV\Niet legger\pomp_gemaal_v13_clipped_streefpeil.shp"
    weir_path = p_folder + r"\WAGV\Niet legger\stuw_v13_clipped_with_do.shp"

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", "AVVBODDR"),
            ("bodemhoogte benedenstrooms", "AVVBODH"),
            ("bodemhoogte bovenstrooms", "AVVBODH"),
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", "IWS_W_WATP"),
            ("hoogte insteek rechterzijde", "IWS_W_WATP"),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("taludhelling linkerzijde", "AVVTALUL"),
            ("taludhelling rechterzijde", "AVVTALUR"),
            ("typeruwheid", "typeruwheid"),  # changed in hydrovak_combined
            ("water_width_index", "IWS_W_WATB"),
        ]
    )

    ## Bridges
    bridge_index_mapping = dict(
        [
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("intreeverlies", "intreeverl"),
            ("ruwheid", "ruwheid"),
            ("typeruwheid", "ruwheidsty"),
            ("uittreeverlies", "uittreever"),
            ("lengte", "lengte"),
        ]
    )

    ## Culverts
    culvert_index_mapping = dict(
        [
            ("breedteopening", "breedteope"),
            ("code", "code"),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "hoogtebinn"),
            ("hoogtebinnenonderkantbov", "hoogtebin0"),
            ("hoogteopening", "hoogteopen"),
            ("intreeverlies", "intreeverl"),
            ("lengte", "lengte"),
            ("ruwheid", "ruwheid"),
            ("typeruwheid", "ruwheidsty"),
            ("uittreeverlies", "uittreever"),
            ("vormkoker", "vormkokeri"),
        ]
    )

    ## Measured profiles
    measured_profile_index_mapping = dict(
        [
            ("code", "code"),
            ("codevolgnummer", "codevolgnu"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("profiel nummer", "metingprof"),
            ("ruwheidhoog", "ruwheidswa"),
            ("ruwheidlaag", "ruwheidsw0"),
            ("type meting", "typebodemi"),
            ("typeruwheid", "ruwheidsty"),
        ]
    )

    ## Pumps
    pump_index_mapping = dict(
        [
            ("code", "code"),
            ("doelvariabele", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("maximalecapaciteit", "maximaleca"),
            ("peil_marge", None),
            ("streefwaarde", "peil1"),
        ]
    )

    ## Weirs
    weir_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", "afvoercoef"),
            ("afvoercoefficient_opening", "afvoerco_1"),
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", None),
            ("hoogstedoorstroomhoogte", "hoogstedo0"),
            ("laagstedoorstroombreedte", "laagstedoo"),
            ("laagstedoorstroomhoogte", "laagstedo0"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", "soortregel"),
            ("soortstuw", "soortstuwi"),
            ("vormopening", None),
        ]
    )
