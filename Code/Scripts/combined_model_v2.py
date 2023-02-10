import sys

from hydrolib.core.io.ext.models import ExtModel, Lateral

sys.path.append("D:\Work\git\GIS_tools\Code")

from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\Combined_test.gpkg"

output_folder = folder + r"\Models\Combined\V6"

config_dhydro = r"combined_config"
config_list = [
    r"hdsr_config",
    r"hhd_config",
    r"hhr_config",
    r"hhsk_config",
    r"wagv_config",
    r"ark_nzk_config",
]

defaults = r"defaults"

build_database = True
build_model = True


if build_database:
    dhd = DHydamoData()
    for ix, config in enumerate(config_list):
        print("\n" + config)

        dhd.from_raw_data(defaults=defaults, config=config)

    dhd.clip_structures_by_branches()
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    lateral = Lateral(
        id="LateralSource_2D_1",
        name="LateralSource_2D_1",
        branchId="hdsr_H012375_0",
        chainage=30,
        discharge=10000,
    )
    extforcefilenew = ExtModel(lateral=[lateral])

    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. load data
    dhd.from_dhydamo_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder, write=False)
    dhd.fm.geometry.usecaching = 1
    dhd.fm.numerics.cflmax = 10
    dhd.fm.output.hisinterval = [0]

    dhd.fm.external_forcing.extforcefilenew = extforcefilenew

    dhd.write_dimr(output_folder=output_folder)
