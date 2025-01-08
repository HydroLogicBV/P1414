import os
from hydrolib.core.io.bc.models import (ForcingBase, ForcingModel, QHTable,
                                        QuantityUnitPair)
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)


class Name:
    name = 'Markermeer'
    
class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 0.1
            node_distance = 50

        class hydrolib_core_options:
            class external_forcing:             
                pass             
                """ *** Forcing is added in combined config script ***
                # Add a timeseries data on the Markermeer nodes
                # The relevant nodes are: "127216.410000_487088.580000", "130578.000000_483044.000000", "137695.530000_482384.620000"
        
                # Provide here the timeseries for the markermeer, per 1 hour
                _timeseries_markermeer = [[i,0] for i in range(24)]
                _markermeer_forcings = []
                _forcing_nodes = ["127216.410000_487088.580000", "130578.000000_483044.000000", "137695.530000_482384.620000"]
                
                for i in range(len(_forcing_nodes)):
                    forcing=ForcingBase(
                        name=_forcing_nodes[i],
                        function="timeseries",
                        quantityunitpair=[QuantityUnitPair(quantity="time", unit=f"hours since 2016-06-01 00:00:00"),
                                        QuantityUnitPair(quantity="waterlevelbnd", unit="m")],
                        datablock=_timeseries_markermeer,
                    )
                    _markermeer_forcings.append(forcing)

                __boundaries = [
                    Boundary(
                        quantity="waterlevelbnd",
                        nodeid="127216.410000_487088.580000",
                        forcingfile=ForcingModel(forcing=_markermeer_forcings),
                    )
                ]
                extforcefilenew = ExtModel(boundary=__boundaries)
                """

class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    # p_folder = r"D:\work\P1414_ROI\Boezemmodel_Waternet_dimr"
    branches_path = p_folder + r"\Markermeer\MarkermeerV4.shp"
    norm_profile_path = p_folder + r"\Markermeer\MarkermeerV4.shp"
    sluice_path = p_folder + r"\Markermeer\Sluis.shp"

    class Peil:
        default_peil = 0

    # output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

    ## Branches
    branch_index_mapping = dict(
        [
            ("code", "Name"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("tunnel", False),
            ("is_duiker", None),
            ("typeruwheid", None),
        ]
    )

    ## Branches
    np_index_mapping = dict(
        [
            ("bodembreedte", "Breedte"),
            ("bodemhoogte benedenstrooms", "H_NAP"),
            ("bodemhoogte bovenstrooms", "H_NAP"),
            ("code", "Name"),
            ("diepte", "Diepte"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogte insteek linkerzijde", None),
            ("hoogte insteek rechterzijde", None),
            ("taludhelling linkerzijde", None),
            ("taludhelling rechterzijde", None),
            ("thalweg offset", "thal_offse"),
            ("typeruwheid", None),
            ("ruwheidhoog", None),
            ("ruwheidlaag", None),
            ("water_width_index", None),
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
            ("flowdir", None),
        ]
    )
