import importlib
import sys
from multiprocessing import Pool

from hydrolib.core.io.ext.models import ExtModel, Lateral

# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)


from data_structures.dhydro_data import DHydroData

def write_model(config_path: str, gpkg_file: str, output_path: str):
    # load configuration file
    model_config = getattr(importlib.import_module("dataset_configs." + config_path), "Models")

    # 1. initialize an instance of DHydamoData
    dhd = DHydamoData()

    # 2. load data
    dhd.from_dhydamo_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_path, output_folder=output_path, write=False)

    if model_config.FM.one_d_bool:
        lateral1 = Lateral(
            id="LateralSource_1D_1",
            name="LateralSource_1D_1",
            branchId="hdsr_H012375",
            chainage=30,
            discharge=2500,
        )
        # lateral2 = Lateral(
        #     id="LateralSource_1D_2",
        #     name="LateralSource_1D_2",
        #     branchId="rijn_DuitseRijn",
        #     chainage=30,
        #     discharge=15000,
        # )
        extforcefilenew = ExtModel(lateral=[lateral1])
        dhd.fm.external_forcing.extforcefilenew = extforcefilenew

    dhd.write_dimr(output_folder=output_path)


def main():

    folder = r"D:\Work\Project\P1414"
    gpkg_file = folder + r"\GIS\HYDAMO\Combined_test_v4.3.gpkg"

    output_folder = folder + r"\Models\Rekentijden\V3"

    config_dhydro = r"rekentijden_v3"

    pool = Pool(processes=2)
    res = []
    for ix in range(12, 14):
        _output_folder = output_folder + r"\{}".format(ix)
        _config_dhydro = config_dhydro + r".{}".format(ix)
        res.append(
            pool.apply_async(
                write_model,
                kwds={
                    "config_path": _config_dhydro,
                    "gpkg_file": gpkg_file,
                    "output_path": _output_folder,
                },
            )
        )

    for r in res:
        r.get()

    pool.join()
    pool.close()


if __name__ == "__main__":
    main()
