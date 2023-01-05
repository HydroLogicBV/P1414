from copy import copy

import geopandas as gpd

from utils_hydamo_profiles import convert_pp_to_hydamo


def hhsk_norm_profiles(input_path: str, output_path: str) -> None:
    # prepare HHSK data for function
    branches_gdf = gpd.read_file(input_path)
    branches_gdf["ruwheidsty"] = 6
    branches_gdf["ruwheidhoo"] = 23.0
    branches_gdf["ruwheidlaa"] = 23.0

    # set column mapping
    index_mapping = dict(
        [
            ("bodembreedte", "BODEMBREED"),
            ("bodemhoogte benedenstrooms", ""),
            ("bodemhoogte bovenstrooms", ""),
            ("hoogte insteek linkerzijde", ""),
            ("hoogte insteek rechterzijde", ""),
            ("taludhelling linkerzijde", "TALUDLINKS"),
            ("taludhelling rechterzijde", "TALUDRECHT"),
            ("typeruwheid", "ruwheidsty"),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("water_width_index", "WATERBREED"),
        ]
    )

    # convert branches with parameters to hydamo format
    out_branches_gdf, hydroobject_normgp, normgeparamprofielwaarde = convert_pp_to_hydamo(
        branches_gdf=branches_gdf,
        index_mapping=index_mapping,
    )

    # save hydamo data in geopackage
    layers = dict(
        [
            ("hydroobject", out_branches_gdf[["CODE", "globalid", "geometry"]]),
            ("hydroobject_normgp", hydroobject_normgp),
            ("normgeparamprofielwaarde", normgeparamprofielwaarde),
        ]
    )

    for name, layer in layers.items():
        layer.to_file(filename=output_path, driver="GPKG", layer=name)


if __name__ == "__main__":
    # set paths
    branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHSK_v4_test.shp"
    norm_profiles_output_path = r"D:\Work\Project\P1414\GIS\HDSR\norm_profielen_test.gpkg"
    hhsk_norm_profiles(input_path=branches_path, output_path=norm_profiles_output_path)
