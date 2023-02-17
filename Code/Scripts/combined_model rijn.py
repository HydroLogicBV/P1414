import sys

from hydrolib.core.io.ext.models import ExtModel, Lateral

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\rhine_combined_test.gpkg"

output_folder = folder + r"\Models\Combined\rijn_V0"

config_dhydro = r"rijntakken_config"
config_list = [r"ark_nzk_config", r"rijntakken_config", r"rijnmaasmonding_config"]
snap_dist_list = [0, 10, 50]

defaults = r"defaults"

build_database = True
build_model = True

if build_database:
    dhd = DHydamoData()
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.from_raw_data(defaults=defaults, config=config, branch_snap_dist=snap_dist_list[ix])

    dhd.clip_structures_by_branches()
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    dhd = DHydamoData()
    dhd.from_dhydamo_gpkg(gpkg_file)
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder, write=False)

    lateral2 = Lateral(
        id="LateralSource_1D_1",
        name="LateralSource_1D_1",
        branchId="rijn_DuitseRijn",
        chainage=30,
        discharge=15000,
    )
    extforcefilenew = ExtModel(lateral=[lateral2])

    dhd.fm.external_forcing.extforcefilenew = extforcefilenew

    dhd.write_dimr(output_folder=output_folder)
