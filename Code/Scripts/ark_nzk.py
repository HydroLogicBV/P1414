import sys

sys.path.append("D:\Work\git\GIS_tools\Code")
# sys.path.append("D:\work\P1414_ROI\GitHub\GIS_tools\Code")
from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
# folder = r"D:\work\P1414_ROI"
gpkg_file = folder + r"\GIS\HYDAMO\ARKNZK.gpkg"
output_folder = folder + r"\Models\ARKNZK\V0"

config = r"ark_nzk_config"
defaults = r"defaults"

build_database = True
build_model = False


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config)
    dhd.clip_structures_by_branches(buffer=10)
    dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)

    zeesluis_list = [
        "ark__sl_Kleine sluis",
        "ark__sl_Middensluis",
        "ark__sl_Noordersluis",
        "ark__sl_Spuisluis",
    ]
    
    for sluis in zeesluis_list:
        stuwid = dhd.ddm.stuw.loc[dhd.ddm.stuw.code == sluis, "globalid"].values[0]
        dhd.ddm.kunstwerkopening.loc[
            dhd.ddm.kunstwerkopening.stuwid == stuwid, "laagstedoorstroomhoogte"
        ] = 6
        dhd.ddm.kunstwerkopening.loc[
            dhd.ddm.kunstwerkopening.stuwid == stuwid, "hoogstedoorstroomhoogte"
        ] = 6

    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config, output_folder=output_folder)
