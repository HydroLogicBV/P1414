import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Waterschap Amstel, Gooi en Vecht'

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "2Dto1D"  # "1Dto2D"
            dx = 500
            dy = 500
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
    p_folder = folder_path_GIS + r"\GIS"
    
    branches_path = p_folder + r"\WAGV\AGV-Update\hydrovakken\hydroobject_v13_clipped_V4.shp"
    #branches_path = dict(
    #    [
    #        ("base", p_folder + r"\WAGV\AGV-Update\hydrovakken\hydroobject_v13_clipped.shp")
    #        #("base", p_folder + r"\Uitgesneden watergangen\AGV_v3.shp"),                # deze moet worden verwijderd?
    #        #("sjoin", p_folder + r"\WAGV\AGV-Update\hydrovakken\hydroobject_v13_clipped.shp"),
    #    ]
    #)

    bridges_path = p_folder + r"\WAGV\brug_v13_clipped\brug_v13_clipped.shp"
    culvert_path = p_folder + r"\WAGV\AGV-Update\duikers\duikersifonhevel_v13_clipped.shp"
    measured_profile_path = (p_folder + r"\WAGV\metingprofielpunt_v13_clipped\metingprofielpunt_v13_clipped.shp")
    norm_profile_path = dict(
        [
            ("base", p_folder + r"\WAGV\hydrovak\hydrovakV2.shp"),
            ("sjoin", p_folder + r"\WAGV\AGV-Update\hydrovakken\hydroobject_v13_clipped_V4.shp"),
        ]
    )
    peil_gebieden_path = p_folder + r"\WAGV\vigerende_peilgebieden\peilgebieden.shp"
    pump_path = p_folder + r"\WAGV\AGV-Update\gemalen\pomp_gemaal_v13_clippedV2.shp"
    
    #pump_path = dict(
    #    [
    #        #("base", p_folder + r"\WAGV\Niet legger\pomp_gemaal_v13_clipped_streefpeil.shp"),       # Nog aanpassen?, nieuwe in \WAGV\AGV-Update\gemalen\pomp_gemaal_v13_clipped.shp
    #        #("concat", p_folder + r"\WAGV\AGV-Update\gemalen\pomp_gemaal_v13_clipped.shp"),
    #        #("concat", p_folder + r"\WAGV\gemaal_v13\gemaal_v13_clipped.shp"),                      # Nog aanpassen?, nieuwe in \WAGV\AGV-Update\gemalen\pomp_gemaal_v13_clipped.shp
    #    ]
    #)

    sluice_path = p_folder + r"\WAGV\AGV-Update\sluizen\sluisV2.shp"                  
    #weir_path = p_folder + r"\WAGV\Niet legger\stuw_v13_clipped_with_do.shp"        # Nog aanpassen?, nieuwe in \WAGV\AGV-Update\stuwen\stuw_v13_clipped.shp 
    weir_path = p_folder + r"\WAGV\AGV-Update\stuwen\stuw_v13_clipped.shp"

    ## Selection criteria
    branch_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    
    class Peil:
        default_peil = -0.4

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("typeruwheid", "ruwheidsty"),  # changed in hydrovak_combined
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
            ("breedteopening", "BREEDTE"),
            ("code", "code"),
            ("geometry", "geometry"),
            ("gesloten", None),
            ("globalid", "globalid"),
            ("hoogtebinnenonderkantbene", "HOOGTE_BE"),
            ("hoogtebinnenonderkantbov", "HOOGTE_BO"),
            ("hoogteopening", "HOOGTE"),
            ("intreeverlies", "intreeverl"),
            ("lengte", "lengte"),
            ("ruwheid", "ruwheid"),
            ("typeruwheid", "typeruwhei"),
            ("uittreeverlies", "uittreever"),
            ("vormkoker", "vormkoker"),
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
            ("streefwaarde", None),              # was "peil1"
        ]
    )
    ## Sluice
    sluice_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", None),
            ("hoogstedoorstroomhoogte", "HOOGTE"),
            ("laagstedoorstroombreedte", "BREEDTE"),
            ("laagstedoorstroomhoogte", "HOOGTE"),
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
            ("afvoercoefficient_opening", "afvoercoef"),    
            ("code", "code"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", None),             
            ("hoogstedoorstroomhoogte", "HOOGTE"),     
            ("laagstedoorstroombreedte", "BREEDTE"),     
            ("laagstedoorstroomhoogte", "HOOGTE"),      
            ("overlaatonderlaat", None),            
            ("soortregelbaarheid", "typeregelb"),         
            ("soortstuw", "typestuw"),                    
            ("vormopening", None),
        ]
    )
