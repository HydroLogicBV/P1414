# Specify here the folder (location) in which the GIS folder is located
folder_path_GIS = r"D:\Work\Project\P1414"

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = 20160601
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 5
            node_distance = 50

        class two_d:
            coupling_type = "1Dto2D"
            dx = 50
            dy = 50
            # elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN4_WSS_filled.TIF"
            initial_peil_raster_path = folder_path_GIS + r"\GIS\peilen\peilen_jp_25m_full.tif"
            roughness_2d_raster_path = (
                folder_path_GIS + r"\GIS\Landgebruik\randstad_nikuradse_roughness_10m.tif"
            )
            two_d_buffer = 100

        class hydrolib_core_options:
            class geometry:
                usecaching = 1

            class numerics:
                cflmax = 0.7

            class output:
                hisinterval = [300]
                wrimap_waterlevel_s0 = False
                wrimap_evaporation = False
                wrimap_velocity_component_u0 = False
                wrimap_upward_velocity_component = False
                wrimap_density_rho = False
                wrimap_horizontal_viscosity_viu = False
                wrimap_horizontal_diffusivity_diu = False
                wrimap_spiral_flow = False
                wrimap_taucurrent = False
                wrimap_chezy = True
                wrimap_turbulence = False
                wrimap_wind = False

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
    sluice_path = _path = dict(
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
            ("boven peil", "ZOMERPEIL"),
            ("geometry", "geometry"),
            ("onder peil", "WINTERPEIL"),
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


class Dambreak:
    ## PATHS
    dambreak_path = folder_path_GIS + r"\GIS\Dijkdoorbraken\Dijkdoorbraak_HM.shp"
    dambreak_index_mapping = dict(
        [
            ("algorithm", None),
            ("breachwidthini", "ini_breedt"),
            ("code", "id"),
            ("crestlevelini", "ini_diepte"),
            ("crestlevelmin", "min_hoogte"),
            ("f1", None),
            ("f2", None),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("t0", "t0"),
            ("timetobreachtomaximumdepth", None),
            ("ucrit", None),
        ]
    )
