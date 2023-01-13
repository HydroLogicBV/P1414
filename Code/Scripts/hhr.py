import sys

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\HHR_clipped.gpkg"
output_folder = folder + r"\Models\HHR\V1"

config = r"hhr_config"
defaults = r"defaults"

build_database = True
build_model = True

# 1. initialize an instance of DHydamoData
dhd = DHydamoData()
if build_database:
    # 2. convert raw data to hydamo data
    dhd.from_raw_data(defaults=defaults, config=config)
    dhd.clip_structures_by_branches()

    # 3. save data to gpkg
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    # # 2. load data
    dhd.from_dhydamo_gpkg(gpkg_file)

    # # 4. save as dhydro model
    dhd.to_dhydro(config=config, output_folder=output_folder)
