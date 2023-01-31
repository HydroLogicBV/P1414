import sys

import geopandas as gpd
import pandas as pd

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\Combined_test.gpkg"

output_folder = folder + r"\Models\Combined\V5_2D_extent"

config_dhydro = r"hdsr_config"
config_list = [r"hdsr_config", r"hhd_config", r"hhr_config", r"hhsk_config", r"wagv_config"]

defaults = r"defaults"

build_database = True
build_model = False

dhd = DHydamoData()
if build_database:
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.from_raw_data(defaults=defaults, config=config)

    dhd.clip_structures_by_branches()
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd.from_dhydamo_gpkg(gpkg_file)
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
