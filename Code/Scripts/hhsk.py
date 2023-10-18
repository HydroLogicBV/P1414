import sys

sys.path.append("D:\Work\git\GIS_tools\Code")
# sys.path.append("D:\work\P1414_ROI\GitHub\GIS_tools\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\HHSK.gpkg"
output_folder = folder + r"\Models\HHSK\V0b"

config = r"hhsk_config"
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

    # 3. save as dhydro model
    dhd.to_dhydro(config=config, output_folder=output_folder)
