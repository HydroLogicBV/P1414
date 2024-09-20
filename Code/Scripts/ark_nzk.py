# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"
import sys
import os
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"
os.environ['GIS_folder_path'] = folder_path_GIS

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\ARKNZK.gpkg"
output_folder = folder_path_output + r"\Models\ARKNZK\V0"

config = r"ark_nzk_config"
defaults = r"defaults"

build_database = True
build_model = False


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config, GIS_folder=folder_path_GIS)
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
