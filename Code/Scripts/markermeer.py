import sys

sys.path.append("D:\Work\git\GIS_tools\Code")
# sys.path.append("D:\work\P1414_ROI\GitHub\GIS_tools\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
# folder = r"D:\work\P1414_ROI"
gpkg_file = folder + r"\GIS\HYDAMO\markermeer.gpkg"
gpkg_file_2 = folder + r"\GIS\HYDAMO\markermeer_hoog.gpkg"
output_folder = folder + r"\Models\markermeer\V0"

config = r"markermeer_config"
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

    for name, wl in dhd.ddm.waterloop.iterrows():
        dhd.ddm.waterloop.loc[name, "peil"] = 3

    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file_2)

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config, output_folder=output_folder)