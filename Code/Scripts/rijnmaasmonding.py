import sys

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\RMM.gpkg"
output_folder = folder + r"\Models\RMM\V0"

config = r"rijnmaasmonding_config"
defaults = r"defaults"

build_database = True
build_model = True

if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. convert raw data to hydamo data
    dhd.from_raw_data(defaults=defaults, config=config)
    dhd.clip_structures_by_branches()

    # 3. save data to gpkg
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. load data
    dhd.from_dhydamo_gpkg(gpkg_file)

    # 3. save as dhydro model
    dhd.to_dhydro(config=config, output_folder=output_folder)
