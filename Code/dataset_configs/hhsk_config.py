class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "1Dto2D"  # "1Dto2D"
            dx = 500
            dy = 500
            # elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            elevation_raster_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS_filled.tif"
            initial_peil_raster_path = r"D:\Work\Project\P1414\GIS\peilen\peilen_jp_25m_full.tif"
            roughness_2d_raster_path = (
                r"D:\Work\Project\P1414\GIS\Landgebruik\randstad_nikuradse_roughness_10m.tif"
            )
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
            ("base", p_folder + r"\Keringen_met_hoogte\hhsk_primaire_kering.shp"),
            ("concat_1", p_folder + r"\Keringen_met_hoogte\hhsk_regionale_kering.shp"),
            ("concat_2", p_folder + r"\Keringen_met_hoogte\hhsk_overige_kering.shp"),
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
    branches_path = p_folder + r"\Uitgesneden watergangen\HHSK_v3.shp"  # From V7
    culvert_path = p_folder + r"\HHSK\Legger\Duiker.shp"
    norm_profile_path = p_folder + r"\HHSK\Legger\Hoofdwatergang.shp"
    peil_gebieden_path = p_folder + r"\HHSK\Legger\Peilvakken.shp"
    pump_path = p_folder + r"\HHSK\Niet legger\Gemaal_peil.shp"
    sluice_path = p_folder + r"\HHSK\Legger\Sluis.shp"
    weir_path = p_folder + r"\HHSK\Legger\Stuw.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
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

    ## Normprofielen
    np_index_mapping = dict(
        [
            ("bodembreedte", "BODEMBREED"),
            ("bodemhoogte benedenstrooms", None),
            ("bodemhoogte bovenstrooms", None),
            ("code", "CODE"),
            ("diepte", ["DIEPTE", "MAXIMALEWA"]),
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

    ## Peil gebied
    peil_index_mapping = dict(
        [
            ("boven peil", "BOVENPEIL"),
            ("geometry", "geometry"),
            ("onder peil", "ONDERPEIL"),
            ("vast peil", "VASTPEIL"),
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
            ("laagstedoorstroomhoogte", None),
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
            ("hoogstedoorstroombreedte", "KRUINBREED"),
            ("hoogstedoorstroomhoogte", ["KRUINHOOGT", "MINIMAALPE", "MAXIMAALPE"]),
            ("laagstedoorstroombreedte", "KRUINBREED"),
            ("laagstedoorstroomhoogte", ["KRUINHOOGT", "MINIMAALPE", "MAXIMAALPE"]),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
        ]
    )
