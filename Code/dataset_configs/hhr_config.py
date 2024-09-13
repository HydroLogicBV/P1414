import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

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
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
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
    p_folder = folder_path_GIS + r"\GIS"
    flood_defences_path = dict(
        [
            ("base", p_folder + r"\Keringen_met_hoogte\hhr_primaire_kering.shp"),
            ("concat_1", p_folder + r"\Keringen_met_hoogte\hhr_regionale_kering.shp"),
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
    p_folder = folder_path_GIS + r"\GIS"
    # p_folder = r"D:\work\P1414_ROI\GIS"
    branches_path = p_folder + r"\Uitgesneden watergangen\HHR_v3.shp"  # From V7
    bridges_path = p_folder + r"\HHRijnland\Niet legger\brug_edited.shp"
    culvert_path = p_folder + r"\HHRijnland\Legger\Duiker\duiker.shp"
    norm_profile_path = p_folder + r"\HHRijnland\Legger\Watergang\Watergang_as_primair.shp"
    peil_gebieden_path = p_folder + r"\HHRijnland\Legger\Peilvakken\gerealiseerde_peilvakken.shp"
    pump_path = p_folder + r"\HHRijnland\Niet legger\gemaal_peil.shp"
    sluice_path = dict(
        [
            ("base", p_folder + r"\HHRijnland\Legger\Sluis\sluis.shp"),
            ("concat", p_folder + r"\HHRijnland\Legger\Noodwaterkering\noodwaterkering.shp"),
        ]
    )
    weir_path = p_folder + r"\HHRijnland\Legger\Stuw\stuw.shp"

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    # branch_selection = dict([("column", "CATEGORIEO"), ("value", "primair")])
    branch_index_mapping = dict(
        [
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("typeruwheid", "TYPERUWHEI"),
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
            ("lengte", "DOORSTRO_1"),
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
            ("hoogtebinnenonderkantbene", "HOOGTEBINN"),
            ("hoogtebinnenonderkantbov", "HOOGTEBI_1"),
            ("hoogteopening", "HOOGTEOPEN"),
            ("intreeverlies", None),
            ("lengte", "LENGTE"),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("vormkoker", "VORMKOKER"),
        ]
    )

    ## Normprofielen
    # np_selection = dict([("column", "CATEGORIEO"), ("value", "primair")])
    np_index_mapping = dict(
        [
            ("bodembreedte", None),
            ("bodemhoogte benedenstrooms", None),
            ("bodemhoogte bovenstrooms", None),
            ("code", "CODE"),
            ("diepte", "WATERDIEPT"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("taludhelling linkerzijde", "TALUDHELLI"),
            ("taludhelling rechterzijde", "TALUDHEL_1"),
            ("typeruwheid", "TYPERUWHEI"),
            ("ruwheidhoog", "RUWHEIDSWA"),
            ("ruwheidlaag", "RUWHEIDSWA"),
            ("water_width_index", "BREEDTE"),
        ]
    )

    ## Peil gebied
    peil_index_mapping = dict(
        [
            ("boven peil", ["ZOMERPEIL", "FLEXZOMERP", "FLEXZOME_1"]),
            ("geometry", "geometry"),
            ("onder peil", ["WINTERPEIL", "FLEXWINTER", "FLEXWINT_1"]),
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
            ("maximalecapaciteit", "MAXIMALECA"),
            ("streefwaarde", "streefpeil"),
            ("peil_marge", None),
        ]
    )

    ## Weirs
    sluice_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "BREEDTE"),
            ("hoogstedoorstroomhoogte", "KERENDEHOO"),
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
            ("hoogstedoorstroomhoogte", "HOOGTECONS"),
            ("laagstedoorstroombreedte", "DOORSTROOM"),
            ("laagstedoorstroomhoogte", "LAAGSTEDOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", "SOORTREGEL"),
            ("soortstuw", "SOORTSTUW"),
            ("vormopening", "KRUINVORM"),
        ]
    )
