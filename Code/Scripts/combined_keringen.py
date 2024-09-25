# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"

import sys
import os

#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"
os.environ['GIS_folder_path'] = folder_path_GIS

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Combined_keringen.gpkg"
gpkgs_list = [
    folder_path_GIS + r"\GIS\HYDAMO\HHSK.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HDSR.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHD.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHR.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\WAGV.gpkg",
]

config_list = [
    r"hhsk_config",
    r"hdsr_config",
    r"hhd_config",
    r"hhr_config",
    r"wagv_config",
]
snap_dist_list = [0, 0, 10, 10, 50]

defaults = r"defaults"

build_database = True

if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        print("\n" + config)
        try:
            dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)
        except AttributeError:
            pass
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)
