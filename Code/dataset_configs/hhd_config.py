class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 5
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
    flood_defences_path = dict(
        [
            ("base", p_folder + r"\Keringen_met_hoogte\hhd_zeewering.shp"),
            ("concat_1", p_folder + r"\Keringen_met_hoogte\hhd_regionale_kering.shp"),
            ("concat_2", p_folder + r"\Keringen_met_hoogte\hhd_polderkade.shp"),
            ("concat_3", p_folder + r"\Keringen_met_hoogte\hhd_landscheiding.shp"),
            ("concat_4", p_folder + r"\Keringen_met_hoogte\hhd_delflandsedijk.shp"),
        ]
    )
    fixed_weir_index_mapping = dict(
        [
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )


class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    # p_folder = r"D:\work\P1414_ROI\GIS"
    branches_path = p_folder + r"\Uitgesneden watergangen\HHD_v3.shp"  # Corrected from V7
    # branches_path = (
    #     p_folder + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water.shp"
    # )
    # bridges_path = p_folder + r"\HDSR\Legger\Bruggen\Bruggen.shp"
    # culvert_path = p_folder + r"\HDSR\Legger\Kokers_Lijnen\Kokers_Lijnen_edited.shp"
    culvert_path = dict(
        [
            (
                "base",
                p_folder
                + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Open duiker.shp",
            ),
            (
                "concat_1",
                p_folder + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Sifon.shp",
            ),
            (
                "concat_2",
                p_folder
                + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Stuwende duiker.shp",
            ),
            (
                "concat_3",
                p_folder
                + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Vispassageduiker.shp",
            ),
            (
                "concat_4",
                p_folder
                + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Inlaatduiker.shp",
            ),
        ]
    )
    norm_profile_path = p_folder + r"\Uitgesneden watergangen\HHD_v3_ww.shp"
    # norm_profile_path = dict(
    #     [
    #         ("base", p_folder + r"\Uitgesneden watergangen\HHD_v3.shp"),
    #         (
    #             "sjoin",
    #             p_folder
    #             + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water_ww.shp",
    #         ),
    #     ]
    # )
    peil_gebieden_path = p_folder + r"\HHDelfland\Peilbesluiten.shp\PeilgebiedPraktijk.shp"
    pump_path = p_folder + r"\HHDelfland\Niet legger\Gemaal_peil.shp"
    sluice_path = p_folder + r"\HHDelfland\Legger_Delfland_shp\Waterkeringen\Sluis.shp"
    watervlak_path = (
        p_folder
        + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Watervoerend deel.shp"
    )
    weir_path = p_folder + r"\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Stuw.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    # branch_index_mapping = dict(
    #     [
    #         ("bodembreedte", None),
    #         ("bodemhoogte benedenstrooms", None),
    #         ("bodemhoogte bovenstrooms", None),
    #         ("code", "CODE"),
    #         ("diepte", "LEGDIEPNUM"),
    #         ("geometry", "geometry"),
    #         ("globalid", "globalid"),
    #         ("hoogte insteek linkerzijde", None),
    #         ("hoogte insteek rechterzijde", None),
    #         ("taludhelling linkerzijde", None),
    #         ("taludhelling rechterzijde", None),
    #         ("typeruwheid", None),
    #         ("ruwheidhoog", None),
    #         ("ruwheidlaag", None),
    #         ("water_width_index", None),
    #     ]
    # )
    branch_index_mapping = dict(
        [
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("typeruwheid", None),
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

    culvert_index_mapping = dict(
        [
            ("breedteopening", "BREEDTEOPE"),  # Check
            ("code", "CODE"),  # Check
            ("geometry", "geometry"),  # Check
            ("gesloten", None),  # Check
            ("globalid", "globalid"),  # Check
            ("hoogtebinnenonderkantbene", ["HOOGBOKBEN", "HOOGTEBI00"]),
            ("hoogtebinnenonderkantbov", ["HOOGBOKBOV", "HOOGTEBINN"]),
            ("hoogteopening", "HOOGTEOPEN"),
            ("intreeverlies", None),  # Check
            ("lengte", "LENGTE"),
            ("typeruwheid", None),  # Check
            ("ruwheid", None),  # Check
            ("uittreeverlies", None),  # Check
            ("vormkoker", ["VORMKOKER_", "VORM"]),  # Check
        ]
    )

    # ## Normprofielen
    np_index_mapping = dict(
        [
            ("bodembreedte", "bodembreed"),
            ("bodemhoogte benedenstrooms", "bodemhoogt"),
            ("bodemhoogte bovenstrooms", "bodemhoo_1"),
            ("code", "code"),
            ("diepte", "diepte"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", "hoogte ins"),
            ("hoogte insteek rechterzijde", "hoogte i_1"),
            ("taludhelling linkerzijde", "taludhelli"),
            ("taludhelling rechterzijde", "taludhel_1"),
            ("typeruwheid", "typeruwhei"),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("water_width_index", None),
        ]
    )

    # ## Normprofielen
    # np_index_mapping = dict(
    #     [
    #         ("bodembreedte", None),
    #         ("bodemhoogte benedenstrooms", None),
    #         ("bodemhoogte bovenstrooms", None),
    #         ("code", "CODE"),
    #         ("diepte", "LEGDIEPNUM"),
    #         ("geometry", "geometry"),
    #         ("globalid", "globalid"),
    #         ("hoogte insteek linkerzijde", None),
    #         ("hoogte insteek rechterzijde", None),
    #         ("taludhelling linkerzijde", None),
    #         ("taludhelling rechterzijde", None),
    #         ("typeruwheid", None),
    #         ("ruwheidhoog", None),
    #         ("ruwheidlaag", None),
    #         ("water_width_index", None),
    #     ]
    # )

    ## Peil gebied
    peil_index_mapping = dict(
        [
            ("boven peil", "WS_HOOGPEI"),
            ("geometry", "geometry"),
            ("onder peil", "WS_LAAGPEI"),
            ("vast peil", None),
        ]
    )

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

    ## Sluice
    sluice_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "BREEDTE"),
            ("hoogstedoorstroomhoogte", None),
            ("laagstedoorstroombreedte", "BREEDTE"),
            ("laagstedoorstroomhoogte", "KERENDEHOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
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
