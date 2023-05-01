import sys

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\RMM_open.gpkg"
output_folder = folder + r"\Models\RMM\V0"

config = r"rijnmaasmonding_open_config"
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

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config, output_folder=output_folder)
