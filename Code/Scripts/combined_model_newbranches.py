import sys

import geopandas as gpd
import pandas as pd

sys.path.append("D:\work\P1414_ROI\GitHub\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

#folder = r"D:\Work\Project\P1414"
folder = r"D:\work\P1414_ROI"
gpkg_file = folder + r"\GIS\HYDAMO\all_combined_test.gpkg"

output_folder = folder + r"\Models\Combined\combined_new_branches_v1"

config_dhydro = r"hdsr_config"
config_list = [r"hdsr_config", r"hhd_config", r"hhr_config", r"hhsk_config", r"wagv_config"]
#config_dhydro = r"rijntakken_config"
#config_list = [r"rijntakken_config", r"rijnmaasmonding_config"]
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
                if key == 'waterloop':
                    continue
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
    # Add waterloop as a gpkg
    new_branches = DHydamoData()
    new_branches.from_dhydamo_gpkg(r'D:\work\P1414_ROI\GIS\test_combined_branches.gpkg')

    # use setattr to set the attributes of dhd
    setattr(dhd.ddm, 'waterloop', getattr(new_branches.ddm, 'waterloop'))
    dhd.clip_structures_by_branches()
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

    

if build_model:
    dhd = DHydamoData()
    dhd.from_dhydamo_gpkg(gpkg_file)
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
