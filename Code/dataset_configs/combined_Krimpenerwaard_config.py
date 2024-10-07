from hydrolib.core.io.bc.models import ForcingBase, ForcingModel, QHTable, QuantityUnitPair
from hydrolib.core.io.ext.models import Boundary, ExtModel, Lateral
import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

STARTTIME = 20000101
STOPTIME = 86400*7

class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = STARTTIME
        stop_time = STOPTIME

        class one_d:
            max_dist_to_struct = 10
            max_snap_dist = 5
            node_distance = 200

        class two_d:
            coupling_type = "1Dto2D"
            dx = 200
            dy = 200
            elevation_raster_path = folder_path_GIS + r"\GIS\AHN\AHN4_WSS_filled.TIF"
            extent_path = (
                folder_path_GIS + r"\GIS\Randstad_shape\dijkringen_randstad_merged.shp"
            )
            initial_peil_raster_path = folder_path_GIS + r"\GIS\peilen\wd_0_v4.tif"
            roughness_2d_raster_path = (
                folder_path_GIS + r"\GIS\Landgebruik\randstad_nikuradse_roughness_10m.tif"
            )
            two_d_buffer = 1000
            clip_extent_path = folder_path_GIS+r"\GIS\peilen\boezem_v2.shp"       # waterways that must be clipped from 2D
            clip_buffer = 0.5*dx                                                     # int, buffer around the polygons that must be clipped in [meters], 0.5*grid size works often to get a sharp clip
            clip_selection = ['> 125 meter', '50 - 125 meter']                    # selection of the buffer layer that must be kept, based on the 'breedtekla' column in the shapefile

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
                __forcings = [__fb1]

                # Add boundary conditions for the Markermeer and Noordzee
                __starttime_str = str(STARTTIME)[0:4]+"-"+str(STARTTIME)[4:6]+"-"+str(STARTTIME)[6:8]
                __n_hours = int(STOPTIME/3600)
                __timeseries_noordzee = [[__i,0] for __i in range(__n_hours)]
                __forcing_nodes = ["100240.000000_497850.000000", "85410.000000_471795.000000", "58820.000000_446520.000000", "49883.000000_431591.000000","49885.000000_431680.000000"]
                
                for __i in range(len(__forcing_nodes)):
                    __forcing=ForcingBase(
                        name=__forcing_nodes[__i],
                        function="timeseries",
                        quantityunitpair=[QuantityUnitPair(quantity="time", unit=f"hours since {__starttime_str} 00:00:00"),
                                        QuantityUnitPair(quantity="waterlevelbnd", unit="m")],
                        datablock=__timeseries_noordzee,
                    )

                    __boundary = Boundary(
                        quantity="waterlevelbnd",
                        nodeid=__forcing_nodes[__i],
                        forcingfile=ForcingModel(forcing=__forcing),
                    )

                    __forcings.append(__forcing)
                    __boundaries.append(__boundary)

                __timeseries_markermeer = [[__i,0] for __i in range(__n_hours)]
                __forcing_nodes = ["127216.410000_487088.580000", "130578.000000_483044.000000", "137695.530000_482384.620000"]
                
                for __i in range(len(__forcing_nodes)):
                    __forcing=ForcingBase(
                        name=__forcing_nodes[__i],
                        function="timeseries",
                        quantityunitpair=[QuantityUnitPair(quantity="time", unit=f"hours since {__starttime_str} 00:00:00"),
                                        QuantityUnitPair(quantity="waterlevelbnd", unit="m")],
                        datablock=__timeseries_markermeer,
                    )
                    __boundary = Boundary(
                        quantity="waterlevelbnd",
                        nodeid=__forcing_nodes[__i],
                        forcingfile=ForcingModel(forcing=__forcing),
                    )

                    __forcings.append(__forcing)
                    __boundaries.append(__boundary)

                __laterals = [
                    Lateral(
                        id="LateralSource_1D_1",
                        name="LateralSource_1D_1",
                        branchId="rijn_wl_DuitseRijn",
                        chainage=30,
                        discharge=16000,
                    )
                ]
                
                # All boundary conditions must be passed as 1 boundary element with multiple forcings, otherwise the boundary conditions overwrite each other in the .bc file
                __final_boundary = Boundary(
                        quantity="waterlevelbnd",
                        nodeid='127216.410000_487088.580000',
                        forcingfile=ForcingModel(forcing=__forcings),
                    )
                extforcefilenew = ExtModel(boundary=__final_boundary, lateral=__laterals)

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


class Dambreak:
    ## PATHS
    dambreak_path = folder_path_GIS + r"\GIS\Dijkdoorbraken\Dijkdoorbraak_Krimpenerwaard.shp"
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
