import sys

from hydrolib.core.io.ext.models import ExtModel, Lateral

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\Combined_test_v7.6.gpkg"
gpkgs_list = [
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHSK_clipped_wnp.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\HDSR_clipped_test.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHD_clipped_wnp.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHR_clipped_wnp.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\WAGV_clipped.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\ARKNZK.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\Rijntakken.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\RMM.gpkg",
]
output_folder = folder + r"\Models\Combined\V7.6"

config_dhydro = r"combined_config"
config_list = [
    r"hhsk_config",
    r"hdsr_config",
    r"hhd_config",
    r"hhr_config",
    r"wagv_config",
    r"ark_nzk_config",
    r"rijntakken_config",
    r"rijnmaasmonding_config",
]
snap_dist_list = [0, 0, 10, 10, 50, 10, 10, 100]

defaults = r"defaults"

build_database = False
load_gpkgs = False
build_model = True


if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.hydamo_from_raw_data(
            defaults=defaults, config=config, branch_snap_dist=snap_dist_list[ix]
        )
        try:
            dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)
        except AttributeError:
            pass

    dhd.clip_structures_by_branches()
    dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults, min_length=500)
    dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults, min_length=500)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if load_gpkgs:
    dhd = DHydroData()
    for ix, gpkg in enumerate(gpkgs_list):
        print("\n" + gpkg)

        # 2. load data
        dhd.hydamo_from_gpkg(gpkg, branch_snap_dist=snap_dist_list[ix])

    dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults, min_length=500)
    dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults, min_length=500)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:

    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.hydamo_from_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
