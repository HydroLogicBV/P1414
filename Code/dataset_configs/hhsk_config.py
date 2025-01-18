import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Hoogheemraadschap Schieland & Krimpenerwaard'

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "1Dto2D"  # "1Dto2D"
            dx = 500
            dy = 500
            # elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN4_WSS_filled.tif"
            extent_path = folder_path_GIS + r"\GIS\Randstad_shape\dijkringen_randstad_merged.shp"
            initial_peil_raster_path = folder_path_GIS + r"\GIS\peilen\peilen_jp_25m_full.tif"
            roughness_2d_raster_path = (
                folder_path_GIS + r"\GIS\Landgebruik\randstad_nikuradse_roughness_10m.tif"
            )
            two_d_buffer = 0

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
            ("base", p_folder + r"\HHSK\Keringen met hoogte\hhsk_primaire_kering.shp"),
            ("concat_1", p_folder + r"\HHSK\Keringen met hoogte\hhsk_regionale_kering.shp"),
            ("concat_2", p_folder + r"\HHSK\Keringen met hoogte\hhsk_overige_kering.shp"),
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
    #branches_path = p_folder + r"\HHSK\Hoofdwatergang\Hoofdwatergangen_model_V5_no_coupures.shp"  
    branches_path = p_folder + r"\HHSK\Hoofdwatergang\Hoofdwatergangen_model_V6.5.shp"        # From V7, renamed the coupures via QGIS such that the name is not duplicated when cut
    culvert_path = p_folder + r"\HHSK\Duiker\Duiker_model2.shp"
    norm_profile_path = p_folder + r"\HHSK\Hoofdwatergang\Hoofdwatergangen_model_V6.5.shp"
    peil_gebieden_path = p_folder + r"\HHSK\Peilvakken\Praktijkpeilgebieden.shp"
    pump_path = p_folder + r"\HHSK\Gemaal\GemaalV2.shp"
    sluice_path = p_folder + r"\HHSK\Sluis\SluisV2.shp"
    weir_path = dict(
        [
            ("base", p_folder + r"\HHSK\Stuw\Stuw.shp"),
            ("concat_1", p_folder + r"\HHSK\Stuw\Coupures_primaire_waterkeringen_V2.shp"),
        ]
    )
    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    # Selection criteria
    branch_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    sluice_selection = dict([("column", "TOEVOEGEN"), ("value", "JA")])
    weir_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    np_selection = dict([("column", "OPNEMEN"), ("value", "JA")])

    class Peil:
        default_peil = -6 # Nodig voor coupures

    
    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "GLOBALIDWA"),
            ("tunnel", False),
            ("is_duiker","IS_DUIKER"),
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
            ("breedteopening", "BREEDTEOPE"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "HOOGTEBINN"),
            ("hoogtebinnenonderkantbov", "HOOGTEBI_1"),
            ("hoogteopening", "HOOGTEOPEN"),   # temporary, must be HOOGTEOPEN, but has invalid values for now
            ("intreeverlies", None),
            ("lengte", "LENGTE"),
            ("typeruwheid", None),
            ("ruwheid", None),
            ("uittreeverlies", None),
            ("vormkoker", "VORMKOKER"),
        ]
    )

    ## Normprofielen
    np_index_mapping = dict(
        [
            ("bodembreedte", "BODEMBREED"),
            ("bodemhoogte benedenstrooms", "BODEMHGTE"),       # bodemhoogte alleen aanwezig voor coupures
            ("bodemhoogte bovenstrooms", "BODEMHGTE"),
            ("code", "CODE"),
            ("diepte", "DIEPTE"), # Was: ["DIEPTE", "MAXIMALEWA"]), -> Aandachtspunt! 
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("taludhelling linkerzijde", "TALUDLINKS"), # Was: "TALUDLINKS"),
            ("taludhelling rechterzijde", "TALUDRECHT"), #Was: "TALUDRECHT"),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", None), # Was: "WATERBREED"),
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
            ("streefwaarde", None), # was: "streefpeil"), # Nog toevoegen aan script
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
            ("hoogstedoorstroomhoogte", "HOOGTE"),
            ("laagstedoorstroombreedte", "BREEDTE"),
            ("laagstedoorstroomhoogte", "HOOGTE"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
            ("flowdir", None),
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
            ("hoogstedoorstroombreedte", ["KRUINBREED", "BREEDTEOPE"]),
            ("hoogstedoorstroomhoogte", ["KRUINHOOGT", "MINIMAALPE", "MAXIMAALPE", "KERENDEHOO"]),
            ("laagstedoorstroombreedte", ["BREEDTEOPE","KRUINBREED"]),
            ("laagstedoorstroomhoogte", ["KRUINHOOGT", "MINIMAALPE", "MAXIMAALPE", "KERENDEHOO"]),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
            ("flowdir", None),
        ]
    )
