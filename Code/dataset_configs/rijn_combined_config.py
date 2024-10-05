from hydrolib.core.io.bc.models import ForcingBase, ForcingModel, QuantityUnitPair
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral
import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = False
        start_time = 20000101
        stop_time = 86400 * 7

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 0.1
            node_distance = 500

        class two_d:
            coupling_type = "2Dto1D"
            dx = 500
            dy = 500
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN_merged.TIF"
            extent_path = folder_path_GIS + r"\GIS\Randstad_shape\dijkringen_randstad_merged.shp"
            two_d_buffer = 100

        class hydrolib_core_options:
            class external_forcing:
                __fb1 = ForcingBase(
                    name="172541.443336_513910.402941",
                    function="constant",
                    quantityunitpair=[QuantityUnitPair(quantity="waterlevelbnd", unit="m\n0")],
                )
                __fb2 = ForcingBase(
                    name="62730.000000_445440.000000",
                    function="constant",
                    quantityunitpair=[QuantityUnitPair(quantity="waterlevelbnd", unit="m\n0")],
                )
                __boundaries = [
                    Boundary(
                        quantity="waterlevelbnd",
                        nodeid="62730.000000_445440.000000",
                        forcingfile=ForcingModel(forcing=[__fb1, __fb2]),
                    ),
                ]
                __laterals = [
                    Lateral(
                        id="LateralSource_1D_1",
                        name="LateralSource_1D_1",
                        branchId="rijn_DuitseRijn",
                        chainage=30,
                        discharge=18000,
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
