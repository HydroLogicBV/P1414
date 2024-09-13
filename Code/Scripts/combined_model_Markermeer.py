# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"

import sys

#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Combined_test_v16_Markermeer.gpkg"
gpkgs_list = [
    folder_path_GIS + r"\GIS\HYDAMO\HHSK.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HDSR.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHD.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHR.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\WAGV.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\ARKNZK.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\Rijntakken.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\RMM_open.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\noordzee.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\markermeer_hoog.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\tunnels.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\Ontbrekende_stuwen.gpkg",
]
output_folder = folder_path_output + r"\Models\Combined\V16_Markermeer"

config_dhydro = r"combined_Markermeer_config"
config_list = [
    r"hhsk_config",
    r"hdsr_config",
    r"hhd_config",
    r"hhr_config",
    r"wagv_config",
    r"ark_nzk_config",
    r"rijntakken_config",
    r"rijnmaasmonding_config",
    r"noordzee_config",
    r"markermeer_config",
]
snap_dist_list = [0, 0, 10, 10, 50, 10, 10, 100, 100, 100, 0, 0]

defaults = r"defaults"

build_database = False
load_gpkgs = True
build_model = True


if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.hydamo_from_raw_data(
            defaults=defaults, config=config, branch_snap_dist=snap_dist_list[ix], GIS_folder=folder_path_GIS
        )
        try:
            dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)
        except AttributeError:
            pass

    dhd.clip_structures_by_branches()
    dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="noordzeekeringen_config", defaults=defaults)
    # dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if load_gpkgs:
    dhd = DHydroData()
    for ix, gpkg in enumerate(gpkgs_list):
        print("\n" + gpkg)

        # 2. load data
        dhd.hydamo_from_gpkg(gpkg, branch_snap_dist=snap_dist_list[ix])

    dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="noordzeekeringen_config", defaults=defaults)
    # dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.hydamo_from_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")
    dhd.ddm.pomp["maximalecapaciteit"] = 0

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
