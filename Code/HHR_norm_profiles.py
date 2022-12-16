from copy import copy

import geopandas as gpd

from utils_hydamo_profiles import convert_pp_to_hydamo

if __name__ == "__main__":
    # prepare AGV data for function
    branches_path = (
        r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Watergang\Watergang_as_clipped_rm.shp"
    )
    norm_profiles_output_path = r"D:\Work\Project\P1414\GIS\HHRijnland\norm_profielen_test.gpkg"

    branches_gdf = gpd.read_file(branches_path)
    branches_gdf["ruwheidsty"] = 4
    # branches_gdf["ruwheidhoo"] = 23.0
    # branches_gdf["ruwheidlaa"] = 23.0

    # set column mapping
    index_mapping = dict(
        [
            ("bodembreedte", "BODEMBREED"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),  # onbekend
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),  # onbekend
            ("hoogte insteek linkerzijde", "IWS_W_INST"),  # onbekend
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),  # onbekend
            ("taludhelling linkerzijde", "TALUDHELLI"),
            ("taludhelling rechterzijde", "TALUDHEL_1"),
            ("typeruwheid", "ruwheidsty"),
            ("ruwheidhoog", "RUWHEIDSWA"),
            ("ruwheidlaag", "RUWHEIDS_1"),
            ("water_width_index", "BREEDTE"),
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
        layer.to_file(filename=norm_profiles_output_path, driver="GPKG", layer=name)
