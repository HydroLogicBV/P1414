# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"

import sys
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)
from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\HHR_HM.gpkg"
output_folder = folder_path_output + r"\Models\HHR\HM_v4"

config = r"hhr_HM_config"
defaults = r"defaults"

build_database = False
build_model = True


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config)
    dhd.clip_structures_by_branches()
    dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)

    # 3. save data to gpkg
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)


if build_model:

    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.hydamo_from_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")
    dhd.ddm.pomp["maximalecapaciteit"] = 0

    # 3. save as dhydro model
    dhd.to_dhydro(config=config, output_folder=output_folder)
