import sys

# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)
from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\tunnels.gpkg"
output_folder = folder_path_output + r"\Models\tunnels\V0"

config_list = [r"tunnel_config", "underpass_config"]
defaults = r"defaults"

build_database = True
build_model = True


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.hydamo_from_raw_data(defaults=defaults, config=config)

    dhd.ddm.waterloop = dhd.ddm.waterloop.loc[dhd.ddm.waterloop["code"] != "tunn_wl_124633861", :]
    dhd.ddm.duiker = dhd.ddm.duiker.loc[dhd.ddm.duiker["code"] != "tunn_cu_124633861", :]

    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config_list[0], output_folder=output_folder)
