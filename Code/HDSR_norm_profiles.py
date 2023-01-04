from typing import List

import geopandas as gpd

from hydamo_helpers import convert_pp_to_hydamo, save_gpkg


def hdsr_norm_profiles(input_path: str, index_mapping: dict, output_path: str):
    branches_gdf = gpd.read_file(input_path)
    branches_gdf["ruwheidsty"] = 6
    branches_gdf["ruwheidhoo"] = 23.0
    branches_gdf["ruwheidlaa"] = 23.0

    hydroobject, hydroobject_normgp, normgeparamprofielwaarde = convert_pp_to_hydamo(
        branches_gdf=branches_gdf, index_mapping=index_mapping
    )

    layers = ["hydroobject", "hydroobject_normgp", "normgeparamprofielwaarde"]
    save_gpkg(
        input_gdfs=[hydroobject, hydroobject_normgp, normgeparamprofielwaarde],
        layers=layers,
        output_path=output_path,
    )


if __name__ == "__main__":
    # set paths
    branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HDSR_v4_test.shp"
    norm_profiles_output_path = r"D:\Work\Project\P1414\GIS\HDSR\norm_profielen_test.gpkg"

    # set column mapping
    index_mapping = dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("hoogte insteek linkerzijde", "IWS_W_INST"),
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
            ("taludhelling linkerzijde", "IWS_W_TALU"),
            ("taludhelling rechterzijde", "IWS_W_TA_1"),
            ("typeruwheid", "ruwheidsty"),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("water_width_index", "IWS_W_WATB"),
        ]
    )
    hdsr_norm_profiles(input_path=branches_path, index_mapping=index_mapping, output_path=norm_profiles_output_path)
