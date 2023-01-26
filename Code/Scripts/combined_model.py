import sys

import geopandas as gpd
import pandas as pd

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\rhine_combined_test.gpkg"

output_folder = folder + r"\Models\Combined\rijn_V0"

# config_dhydro = r"hdsr_config"
# config_list = [r"hdsr_config", r"hhd_config", r"hhr_config", r"hhsk_config", r"wagv_config"]
config_dhydro = r"rijntakken_config"
config_list = [r"rijntakken_config", r"rijnmaasmonding_config"]
defaults = r"defaults"

build_database = True
build_model = True

if build_database:
    for ix, config in enumerate(config_list):
        print("\n" + config)

        # load raw data correspondingn to config file
        _dhd = DHydamoData()
        _dhd.from_raw_data(defaults=defaults, config=config)

        # if first, simply set data
        if ix == 0:
            dhd = _dhd

        # Else, append data to each geodataframe in the DHydamo Data Model
        else:
            for key, value in _dhd.ddm.__dict__.items():
                if value is not None:
                    if getattr(dhd.ddm, key) is None:
                        setattr(dhd.ddm, key, getattr(_dhd.ddm, key))
                    else:
                        new_gdf = gpd.GeoDataFrame(
                            pd.concat(
                                [getattr(dhd.ddm, key), getattr(_dhd.ddm, key)], ignore_index=True
                            )
                        )
                        setattr(dhd.ddm, key, new_gdf)

    dhd.clip_structures_by_branches()
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd = DHydamoData()
    dhd.from_dhydamo_gpkg(gpkg_file)
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
