import sys
import os

# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)
sys.path.append(r'C:\Work\HL-P24050\P1414\HYDROLIB_adapted\hydrolib')
sys.path.append(r'C:\Work\HL-P24050\P1414\HYDROLIB_adapted')
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\Code")
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\HYDROLIB_adapted")
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\HYDROLIB_adapted\hydrolib")

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database"
folder_path_output = r"P:\HL-P24050\05_Analysis\02_Model"

os.environ['GIS_folder_path'] = folder_path_GIS

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Buitenwater_11_februari.gpkg"

gpkgs_list = [
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\Buitenwater_11_februari.gpkg",              'snap_dist': 10},
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\HHSK_11_februari.gpkg",                     'snap_dist': 10},
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\HDSR_11_februari.gpkg",                     'snap_dist': 35},
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\HHD_15_januari.gpkg",                       'snap_dist': 10},
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\HHR_11_februari.gpkg",                      'snap_dist': 10},
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\WAGV_14_januari.gpkg",                      'snap_dist': 10},
    {'gpkg_file': folder_path_GIS + r"\GIS\HYDAMO\Tunnels_50m_10_januari.gpkg",               'snap_dist': 0},
]

# Model name
output_folder = folder_path_output + r"\Rivieren_V3.16"

# Existing meshes available
existing_meshes = {
    '500m': r"P:\HL-P24050\05_Analysis\02_Model\Combined_V3.16_50m\dflowfm\network.nc",
    '100m': r"P:\HL-P24050\05_Analysis\02_Model\Combined_V3.4_100m\dflowfm\network.nc",
    '50m': r"P:\HL-P24050\05_Analysis\02_Model\Combined_V3.3_50m\dflowfm\network.nc",
    'None': None
}

config_dhydro = r"combined_WBD_config"
config_list = [
    # {'config': r"ark_nzk_config",               'snap_dist': 10},
    # {'config': r"rijnmaasmonding_open_config",  'snap_dist': 10},
    # {'config': r"rijntakken_config",            'snap_dist': 10},
    # {'config': r"noordzee_config",              'snap_dist': 200},
    # {'config': r"markermeer_config",            'snap_dist': 10},
    # {'config': r"hhsk_config",                  'snap_dist': 10},
    {'config': r"hdsr_config",                  'snap_dist': 20},
    # {'config': r"hhd_config",                   'snap_dist': 10},
    # {'config': r"hhr_config",                   'snap_dist': 10},
    # {'config': r"wagv_config",                  'snap_dist': 10},
    {'config': r"ontbrekende_stuwen_config",    'snap_dist': 10},
    # {'config': r"randvoorwaarden_config",       'snap_dist': 0},
    # {'config': r"tunnel_config",                'snap_dist': 1},
    # {'config': r"underpass_config",             'snap_dist': 1},
]

defaults = r"defaults"

build_database = False
load_gpkgs = False
build_model = True
if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        

        dhd.hydamo_from_raw_data(
            defaults=defaults, config=config['config'], branch_snap_dist=config['snap_dist'], GIS_folder=folder_path_GIS, dhydro_config=config_dhydro
        )
        try:
            dhd.fixed_weirs_from_raw_data(config=config['config'], defaults=defaults)
        except AttributeError:
            pass

    dhd.clip_structures_by_branches(buffer=5)
    #dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults)
    #dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults)
    #dhd.fixed_weirs_from_raw_data(config="noordzeekeringen_config", defaults=defaults)
    #dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if load_gpkgs:
    dhd = DHydroData()
    for ix, gpkg in enumerate(gpkgs_list):
        print("\n" + gpkg['gpkg_file'])

        # 2. load data
        dhd.hydamo_from_gpkg(gpkg['gpkg_file'], branch_snap_dist=gpkg['snap_dist'])#, GIS_folder=folder_path_GIS)

    dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults, min_length = 50)
    dhd.fixed_weirs_from_raw_data(config="noordzeekeringen_config", defaults=defaults)

    # dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.hydamo_from_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    if "brug" in dhd.features:
        dhd.features.remove("brug")
    if "gemaal" in dhd.features:
        dhd.ddm.pomp["maximalecapaciteit"] = 0

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, 
                    load_mesh2d_path = existing_meshes['50m'], 
                    output_folder=output_folder
                )
