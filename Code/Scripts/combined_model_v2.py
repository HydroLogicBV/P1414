import sys

from hydrolib.core.io.ext.models import ExtModel, Lateral

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\Combined_test_v5.gpkg"

output_folder = folder + r"\Models\Combined\V5"

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
build_model = True


if build_database:
    dhd = DHydamoData()
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.from_raw_data(defaults=defaults, config=config, branch_snap_dist=snap_dist_list[ix])

    dhd.clip_structures_by_branches()
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:

    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. load data
    dhd.from_dhydamo_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder, write=False)

    lateral1 = Lateral(
        id="LateralSource_1D_1",
        name="LateralSource_1D_1",
        branchId="hdsr_H012375",
        chainage=30,
        discharge=2500,
    )
    lateral2 = Lateral(
        id="LateralSource_1D_2",
        name="LateralSource_1D_2",
        branchId="rijn_DuitseRijn",
        chainage=30,
        discharge=15000,
    )
    extforcefilenew = ExtModel(lateral=[lateral1, lateral2])
    dhd.fm.external_forcing.extforcefilenew = extforcefilenew

    dhd.write_dimr(output_folder=output_folder)
