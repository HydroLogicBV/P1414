from hydrolib.core.io.bc.models import (ForcingBase, ForcingModel,
                                        QuantityUnitPair)
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral
import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Hoogheemraadschap De Stichtse Rijnlanden'
    
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
            coupling_type = "2Dto1D"
            dx = 500
            dy = 500
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            initial_peil_raster_path = folder_path_GIS + r"\GIS\peilen\peilen_jp_25m_full.tif"
            two_d_buffer = 100

        class hydrolib_core_options:
            class external_forcing:
                __fb1 = ForcingBase(
                    name="109975.000000_446940.000000",
                    function="constant",
                    quantityunitpair=[QuantityUnitPair(quantity="waterlevelbnd", unit="m\n0")],
                )
                __boundaries = [
                    Boundary(
                        quantity="waterlevelbnd",
                        nodeid="109975.000000_446940.000000",
                        forcingfile=ForcingModel(forcing=[__fb1]),
                    ),
                ]
                __laterals = [
                    Lateral(
                        id="LateralSource_1D_1",
                        name="LateralSource_1D_1",
                        branchId="hdsr_wl_H012375",
                        chainage=30,
                        discharge=3000,
                    )
                ]
                extforcefilenew = ExtModel(boundary=__boundaries, lateral=__laterals)

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
    flood_defences_path = p_folder + r"\HDSR\Keringen met hoogte\hdsr.shp"

    fixed_weir_index_mapping = dict(
        [
            ("code", "OBJECTID"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
        ]
    )


class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    branches_path = p_folder + r"\HDSR\HydroObject\HydroObject_v2.shp"  # Corrected from V7
    #bridges_path = p_folder + r"\HDSR\Bruggen\Bruggen.shp"
    culvert_path = p_folder + r"\HDSR\Kokers_Lijnen\Kokers_Lijnen.shp"
    norm_profile_path = p_folder + r"\HDSR\HydroObject\HydroObject_v2.shp"
    # norm_profile_path = dict(
    #     [
    #         ("base", p_folder + r"\Uitgesneden watergangen\HDSR_v3.shp"),
    #         ("sjoin", p_folder + r"\HDSR\Legger\Hydro_Objecten(2)\HydroObject_primair.shp"),
    #     ]
    # )
    peil_gebieden_path = p_folder + r"\HDSR\BR_Peilgebieden\BR_Peilgebieden.shp"
    pump_path = p_folder + r"\HDSR\Gemalen\Gemalen.shp"
    sluice_path = p_folder + r"\HDSR\Sluizen_Lijnen\Sluizen_Lijnen.shp"
    weir_path = dict(
        [
            ("base", p_folder + r"\HDSR\BR_Stuwen\BR_Stuwen.shp"),
            ("concat", p_folder + r"\HDSR\Kokers_Lijnen_Kerend\Kokers_Lijnen_Kerend.shp"),
        ]
    )

    # Selection criteria
    branch_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    np_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    sluice_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    weir_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    pump_selection = dict([("column", "OPNEMEN"), ("value", "JA")])
    culvert_selection = dict([("column", "OPNEMEN"), ("value", "JA")])

    ## Branches
    branch_index_mapping = dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", "IWS_W_INST"),
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
            ("taludhelling linkerzijde", "IWS_W_TALU"),
            ("taludhelling rechterzijde", "IWS_W_TA_1"),
            ("tunnel", False),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", "IWS_W_WATB"),
        ]
    )

    ## Bridges
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

    ## normprofiles
    np_index_mapping = dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("code", "CODE"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
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
    ## Peil gebied
    peil_index_mapping = dict(
        [
            ("boven peil", ["ZOMERPEIL", "BOVENPEIL"]),
            ("geometry", "geometry"),
            ("onder peil", ["WINTERPEIL", "ONDERPEIL"]),
            ("vast peil", "VASTPEIL"),
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
            ("streefwaarde", None), # Was: "streefpeil"),
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
            ("hoogstedoorstroombreedte", "DOORVAARTB"),
            ("hoogstedoorstroomhoogte", None),
            ("laagstedoorstroombreedte", "DOORVAARTB"),
            ("laagstedoorstroomhoogte", "KERENDEHOO"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", "SOORTSLUIS"),
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
            ("soortregelbaarheid", "SOORTREGEL"),
            ("soortstuw", "SOORTSTUW"),
            ("vormopening", None),
        ]
    )
