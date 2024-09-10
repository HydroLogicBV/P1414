import sys

# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)
from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Rijntakken.gpkg"
output_folder = folder_path_output + r"\Models\Rijntakken\V0"

config = r"rijntakken_config"
defaults = r"defaults"

build_database = True
build_model = False


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config)
    dhd.clip_structures_by_branches()

    arksluis_list = [
        "rijn_we_AR_62.0_C_HK_Keerschuif-Ravenswaaij",
        "rijn_we_AR_61.9_L_SS_Prinses-Marijkesluis-west",
    ]

    for sluis in arksluis_list:
        stuwid = dhd.ddm.stuw.loc[dhd.ddm.stuw.code == sluis, "globalid"].values[0]
        dhd.ddm.kunstwerkopening.loc[
            dhd.ddm.kunstwerkopening.stuwid == stuwid, "laagstedoorstroomhoogte"
        ] = 9
        dhd.ddm.kunstwerkopening.loc[
            dhd.ddm.kunstwerkopening.stuwid == stuwid, "hoogstedoorstroomhoogte"
        ] = 9

    arksluis_list = [
        "rijn_we_AR_71.2_L_SS_Prins-Bernhardsluis-oost",
        "rijn_we_AR_71.2_R_SS_Prins-Bernhardsluis-west",
    ]

    for sluis in arksluis_list:
        stuwid = dhd.ddm.stuw.loc[dhd.ddm.stuw.code == sluis, "globalid"].values[0]
        dhd.ddm.kunstwerkopening.loc[
            dhd.ddm.kunstwerkopening.stuwid == stuwid, "laagstedoorstroomhoogte"
        ] = 11
        dhd.ddm.kunstwerkopening.loc[
            dhd.ddm.kunstwerkopening.stuwid == stuwid, "hoogstedoorstroomhoogte"
        ] = 11

    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config, output_folder=output_folder)
