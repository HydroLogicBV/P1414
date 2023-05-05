import sys

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\HHSKRMM.gpkg"
gpkgs_list = [
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHSK.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\Rijntakken.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\RMM_HIJ_closed.gpkg",
]
output_folder = folder + r"\Models\HHSK\HHSKRMM"

config_dhydro = r"hhsk_config"
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
snap_dist_list = [0, 10, 100]

defaults = r"defaults"

load_gpkgs = True
build_model = True


if load_gpkgs:
    dhd = DHydroData()
    for ix, gpkg in enumerate(gpkgs_list):
        print("\n" + gpkg)

        # 2. load data
        dhd.hydamo_from_gpkg(gpkg, branch_snap_dist=snap_dist_list[ix])

    # dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.hydamo_from_gpkg(gpkg_file)

    # remove brug as it needs a cs
    # del dhd.ddm.brug
    # dhd.features.remove("brug")
    dhd.ddm.pomp["maximalecapaciteit"] = 0

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
