# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"

import sys
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\noordzee_hoog.gpkg"
gpkg_file_2 = folder_path_GIS + r"\GIS\HYDAMO\noordzee_hoog_KW_open.gpkg"
output_folder = folder_path_output + r"\Models\noordzee\V0"

config = r"noordzee_hoog_config"
defaults = r"defaults"

build_database = True
build_model = False


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config)

    dhd.clip_structures_by_branches()
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

    del dhd.ddm.stuw, dhd.ddm.kunstwerkopening, dhd.ddm.regelmiddel
    dhd.features.remove("stuw")
    dhd.features.remove("kunstwerkopening")
    dhd.features.remove("regelmiddel")
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file_2)

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config, output_folder=output_folder)
