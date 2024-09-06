# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"

import sys

#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"D:\Work\Project\P1414"
folder_path_output = r"D:\Work\Project\P1414"

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\rhine_combined_test.gpkg"

output_folder = folder_path_output + r"\Models\Rijn\V0"

config_dhydro = r"rijn_combined_config"
config_list = [r"rijntakken_config", r"rijnmaasmonding_config"]
snap_dist_list = [10, 50]

defaults = r"defaults"

build_database = False
build_model = True

if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.hydamo_from_raw_data(
            defaults=defaults, config=config, branch_snap_dist=snap_dist_list[ix]
        )

    dhd.clip_structures_by_branches()
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd = DHydroData()
    dhd.hydamo_from_gpkg(gpkg_file)
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder, write=False)

    dhd.write_dimr(output_folder=output_folder)
