import sys

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\HHD_selection.gpkg"
output_folder = folder + r"\Models\HHD\V1"

config = r"hhd_config"
defaults = r"defaults"

# 1. initialize an instance of DHydamoData
dhd = DHydamoData()

# 2. convert raw data to hydamo data
dhd.from_raw_data(defaults=defaults, config=config)

# # 2. load data
# dhd.from_dhydamo_gpkg(gpkg_file)

# 3. save data to gpkg
dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

# # 4. save as dhydro model
# dhd.to_dhydro(config=config, output_folder=output_folder)
