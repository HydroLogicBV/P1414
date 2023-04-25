import sys

sys.path.append("D:\Work\git\GIS_tools\Code")
# sys.path.append("D:\work\P1414_ROI\GitHub\GIS_tools\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\Ontbrekende_stuwen.gpkg"


config = r"ontbrekende_stuwen_config"
defaults = r"defaults"

build_database = True


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config)

    # 3. save data to gpkg
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)
