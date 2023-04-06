import sys

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydro_data import DHydroData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\Combined_keringen.gpkg"
gpkgs_list = [
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHSK.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\HDSR.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHD.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\HHR.gpkg",
    r"D:\Work\Project\P1414\GIS\HYDAMO\WAGV.gpkg",
]

config_list = [
    r"hhsk_config",
    r"hdsr_config",
    r"hhd_config",
    r"hhr_config",
    r"wagv_config",
]
snap_dist_list = [0, 0, 10, 10, 50]

defaults = r"defaults"

build_database = True

if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        print("\n" + config)
        try:
            dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)
        except AttributeError:
            pass
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)
