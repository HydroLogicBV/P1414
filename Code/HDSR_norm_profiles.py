import geopandas as gpd

from utils_hydamo_profiles import convert_pp_to_hydamo


def hdsr_norm_profiles(input_path: str, output_path: str) -> None:
    branches_gdf = gpd.read_file(input_path)
    branches_gdf["ruwheidsty"] = 6
    branches_gdf["ruwheidhoo"] = 23.0
    branches_gdf["ruwheidlaa"] = 23.0

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

    # convert branches with parameters to hydamo format
    out_branches_gdf, hydroobject_normgp, normgeparamprofielwaarde = convert_pp_to_hydamo(
        branches_gdf=branches_gdf,
        index_mapping=index_mapping,
    )

    # save hydamo data in geopackage
    layers = dict(
        [
            ("hydroobject", out_branches_gdf[["code", "globalid", "geometry"]]),
            ("hydroobject_normgp", hydroobject_normgp),
            ("normgeparamprofielwaarde", normgeparamprofielwaarde),
        ]
    )
    for name, layer in layers.items():
        layer.to_file(filename=output_path, driver="GPKG", layer=name)


if __name__ == "__main__":
    # set paths
    branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HDSR.shp"
    norm_profiles_output_path = r"D:\Work\Project\P1414\GIS\HDSR\norm_profielen_test.gpkg"
    hdsr_norm_profiles(input_path=branches_path, output_path=norm_profiles_output_path)
