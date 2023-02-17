import sys

from hydrolib.core.io.ext.models import ExtModel, Lateral

sys.path.append("D:\Work\git\GIS_tools\Code")
from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HYDAMO\HDSR_clipped_test.gpkg"
output_folder = folder + r"\Models\HDSR\V00"

config = r"hdsr_config"
defaults = r"defaults"

build_database = False
build_model = True


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. convert raw data to hydamo data
    dhd.from_raw_data(defaults=defaults, config=config)
    dhd.clip_structures_by_branches()

    # 3. save data to gpkg
    dhd.to_dhydamo_gpkg(output_gpkg=gpkg_file)

if build_model:
    lateral = Lateral(
        id="LateralSource_2D_1",
        name="LateralSource_2D_1",
        branchId="hdsr_H012375",
        chainage=30,
        discharge=10000,
    )
    extforcefilenew = ExtModel(lateral=[lateral])

    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. load data
    dhd.from_dhydamo_gpkg(gpkg_file)

    # 3. save as dhydro model
    dhd.to_dhydro(config=config, output_folder=output_folder, write=False)
    dhd.fm.geometry.usecaching = 1
    dhd.fm.numerics.cflmax = 10
    dhd.fm.output.hisinterval = [0]

    dhd.fm.external_forcing.extforcefilenew = extforcefilenew

    dhd.write_dimr(output_folder=output_folder)
