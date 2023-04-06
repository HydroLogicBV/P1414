from hydrolib.core.io.bc.models import (ForcingBase, ForcingModel, QHTable,
                                        QuantityUnitPair)
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral


class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 10
            max_snap_dist = 5
            node_distance = 500

        class two_d:
            coupling_type = "1Dto2D"
            dx = 500
            dy = 500
            elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN4_WSS.TIF"
            extent_path = "D:\Work\Project\P1414\GIS\Randstad_shape\dijkringen_randstad_merged.shp"
            initial_peil_raster_path = r"D:\Work\Project\P1414\GIS\peilen\peilen_jp_25m_full.tif"
            roughness_2d_raster_path = (
                r"D:\Work\Project\P1414\GIS\Landgebruik\randstad_nikuradse_roughness_10m.tif"
            )
            two_d_buffer = 1000

        class hydrolib_core_options:
            class external_forcing:
                # __fb1 = ForcingBase(
                #     name="172541.443336_513910.402941",
                #     function="constant",
                #     quantityunitpair=[QuantityUnitPair(quantity="waterlevelbnd", unit="m")],
                #     datablock=[[0]],
                # )
                __fb1 = QHTable(
                    name="172541.443336_513910.402941",
                    quantityunitpair=[
                        QuantityUnitPair(quantity="qhbnd waterlevel", unit="m"),
                        QuantityUnitPair(quantity="qhbnd discharge", unit="m³/s"),
                    ],
                    datablock=[
                        [-0.369, 50],
                        [-0.356, 100],
                        [-0.342, 150],
                        [-0.325, 200],
                        [-0.306, 250],
                        [-0.286, 300],
                        [-0.264, 350],
                        [-0.242, 400],
                        [-0.219, 450],
                        [-0.196, 500],
                        [-0.173, 550],
                        [-0.15, 600],
                        [-0.127, 650],
                        [-0.103, 700],
                        [-0.082, 750],
                        [-0.061, 800],
                        [-0.041, 850],
                        [-0.022, 900],
                        [-0.004, 950],
                        [0.013, 1000],
                        [0.029, 1050],
                        [0.045, 1100],
                        [0.059, 1150],
                        [0.073, 1200],
                        [0.086, 1250],
                        [0.099, 1300],
                        [0.111, 1350],
                        [0.122, 1400],
                        [0.132, 1450],
                        [0.143, 1500],
                        [0.152, 1550],
                        [0.162, 1600],
                        [0.171, 1650],
                        [0.179, 1700],
                        [0.188, 1750],
                        [0.196, 1800],
                        [0.204, 1850],
                        [0.212, 1900],
                        [0.219, 1950],
                        [0.227, 2000],
                        [0.234, 2050],
                        [0.242, 2100],
                        [0.249, 2150],
                        [0.257, 2200],
                        [0.264, 2250],
                        [0.272, 2300],
                        [0.279, 2350],
                        [0.287, 2400],
                        [0.295, 2450],
                        [0.303, 2500],
                        [0.311, 2550],
                        [0.32, 2600],
                        [0.328, 2650],
                        [0.337, 2700],
                        [0.346, 2750],
                        [0.355, 2800],
                        [0.364, 2850],
                        [0.374, 2900],
                        [0.384, 2950],
                        [0.394, 3000],
                        [0.404, 3050],
                        [0.414, 3100],
                        [0.425, 3150],
                        [0.832, 5000],
                    ],
                )
                __boundaries = [
                    Boundary(
                        quantity="qhbnd",
                        nodeid="172541.443336_513910.402941",
                        forcingfile=ForcingModel(forcing=[__fb1]),
                    ),
                ]
                __laterals = [
                    Lateral(
                        id="LateralSource_1D_1",
                        name="LateralSource_1D_1",
                        branchId="rijn_wl_DuitseRijn",
                        chainage=30,
                        discharge=3000,
                    )
                ]
                extforcefilenew = ExtModel(boundary=__boundaries, lateral=__laterals)

            class geometry:
                usecaching = 1

            class numerics:
                cflmax = 0.7

            class output:
                hisinterval = [0]
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
