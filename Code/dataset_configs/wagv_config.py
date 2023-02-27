class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "2Dto1D"  # "1Dto2D"
            dx = 500
            dy = 500
            # elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            elevation_raster_path = "D:\work\P1414_ROI\GIS\AHN\AHN_merged.TIF"
            two_d_buffer = 100

        class hydrolib_core_options:
            class geometry:
                dxmin1d = 500
                usecaching = 1

            class numerics:
                cflmax = 0.7

            class output:
                hisinterval = [0]

            class time:
                dtmax = 60


class FixedWeirs:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    flood_defences_path = p_folder + r"\Keringen_met_hoogte\wagv.shp"
    fixed_weir_index_mapping = dict(
        [
            ("code", "DWKNAAM"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )


class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    # p_folder = r"D:\work\P1414_ROI\GIS"
    # branches_path = p_folder + r"\WAGV\Niet legger\hydrovak_combined_v10.shp"
    # branches_path = [
    #     p_folder + r"\Uitgesneden watergangen\AGV_v00_test.shp",
    #     p_folder + r"\WAGV\hydrovak\hydrovak.shp",
    # ]
    branches_path = [
        p_folder + r"\Uitgesneden watergangen\AGV_v2.2_test.shp",
        p_folder + r"\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp",
    ]

    bridges_path = p_folder + r"\WAGV\brug_v13\brug_v13_clipped.shp"
    culvert_path = p_folder + r"\WAGV\duikersifonhevel_v13\duikersifonhevel_v13_clipped.shp"
    measured_profile_path = (
        p_folder + r"\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped_rm.shp"
    )
    norm_profile_path = [
        p_folder + r"\WAGV\hydrovak\hydrovak.shp",
        p_folder + r"\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp",
    ]
    peil_gebieden_path = p_folder + r"\WAGV\vigerende_peilgebieden\peilgebieden.shp"
    pump_path = p_folder + r"\WAGV\Niet legger\pomp_gemaal_v13_clipped_streefpeil.shp"
    sluice_path = p_folder + r"\WAGV\sluis\sluis.shp"
    weir_path = p_folder + r"\WAGV\Niet legger\stuw_v13_clipped_with_do.shp"

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", None),
            ("bodemhoogte benedenstrooms", None),
            ("bodemhoogte bovenstrooms", None),
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("taludhelling linkerzijde", None),
            ("taludhelling rechterzijde", None),
            ("typeruwheid", "ruwheidsty"),  # changed in hydrovak_combined
            ("water_width_index", None),
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

    ## Norm Profiles
    np_index_mapping = dict(
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
            ("typeruwheid", "ruwheidsty"),  # changed in hydrovak_combined
            ("water_width_index", "IWS_W_WATB"),
        ]
    )
    ## Peil gebied
    peil_index_mapping = dict(
        [
            ("boven peil", ["FLEXIBEL_1", "ZOMERPEIL"]),
            ("geometry", "geometry"),
            ("onder peil", ["FLEXIBEL_W", "WINTERPEIL"]),
            ("vast peil", "VAST_PEIL"),
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
    ## Sluice
    sluice_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "KSLIDENT"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", None),
            ("hoogstedoorstroomhoogte", "KWKKERHG"),
            ("laagstedoorstroombreedte", "KSLDVBRE"),
            ("laagstedoorstroomhoogte", "KWKKERHG"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
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
