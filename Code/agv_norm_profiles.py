from copy import copy

import geopandas as gpd

from utils_hydamo_profiles import convert_pp_to_hydamo

if __name__ == "__main__":
    # prepare AGV data for function
    branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\AGV.shp"
    hydrovak_path = r"D:\Work\Project\P1414\GIS\WAGV\hydrovak\hydrovak.shp"
    norm_profiles_output_path = r"D:\Work\Project\P1414\GIS\WAGV\norm_profielen_test.gpkg"

    branches_gdf = gpd.read_file(branches_path)
    hydrovak_gdf = gpd.read_file(hydrovak_path)

    branches_gdf_new = copy(
        branches_gdf[["code", "ruwheidsty", "ruwheidhoo", "ruwheidlaa", "geometry"]]
    )
    branches_gdf_new["globalid"] = ""
    branches_gdf_new = branches_gdf_new.merge(
        hydrovak_gdf.drop(columns="geometry"), how="inner", left_on="code", right_on="OVKIDENT"
    )

    # set column mapping
    index_mapping = dict(
        [
            ("bodembreedte", "AVVBODDR"),
            ("bodemhoogte benedenstrooms", "AVVBODH"),
            ("bodemhoogte bovenstrooms", "AVVBODH"),
            ("hoogte insteek linkerzijde", "IWS_W_WATP"),
            ("hoogte insteek rechterzijde", "IWS_W_WATP"),
            ("taludhelling linkerzijde", "AVVTALUL"),
            ("taludhelling rechterzijde", "AVVTALUR"),
            ("typeruwheid", "ruwheidsty"),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("water_width_index", "IWS_W_WATB"),
        ]
    )

    # convert branches with parameters to hydamo format
    branches_gdf, hydroobject_normgp, normgeparamprofielwaarde = convert_pp_to_hydamo(
        branches_gdf=branches_gdf_new,
        index_mapping=index_mapping,
    )

    # save hydamo data in geopackage
    layers = dict(
        [
            ("hydroobject", branches_gdf[["code", "globalid", "geometry"]]),
            ("hydroobject_normgp", hydroobject_normgp),
            ("normgeparamprofielwaarde", normgeparamprofielwaarde),
        ]
    )

    for name, layer in layers.items():
        layer.to_file(filename=norm_profiles_output_path, driver="GPKG", layer=name)
