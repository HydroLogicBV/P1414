import uuid
from copy import copy

import geopandas as gpd
import numpy as np
from tqdm import tqdm

roughness_mapping = {
    "Chezy": "Chezy",
    "Manning": "Manning",
    "StricklerKn": "StricklerNikuradse",
    "StricklerKs": "Strickler",
    "White Colebrook": "WhiteColebrook",
    "Bos en Bijkerk": "deBosBijkerk",
    "Onbekend": "Strickler",
    "Overig": "Strickler",
}
roughness_mapping = list(roughness_mapping)
dh = 0.5  # height difference between peil and insteekhoogte

branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\AGV.shp"
hydrovak_path = r"D:\Work\Project\P1414\GIS\WAGV\hydrovak\hydrovak.shp"
norm_profiles_output_path = r"D:\Work\Project\P1414\GIS\WAGV\norm_profielen.gpkg"

branches_gdf = gpd.read_file(branches_path)
hydrovak_gdf = gpd.read_file(hydrovak_path)

branches_gdf_new = copy(branches_gdf)
branches_gdf_new["globalid"] = ""
ho_ngp_list = []
ngp_list = []

for ix, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
    hydrovak = hydrovak_gdf.loc[hydrovak_gdf["OVKIDENT"] == branch["code"]]

    branch_gid = str(uuid.uuid4())
    branches_gdf_new.loc[ix, "globalid"] = branch_gid

    ngp_gid = str(uuid.uuid4())
    ngp = dict(
        [
            ("hydroobjectid", branch_gid),
            ("normgeparamprofielid", ngp_gid),
            ("globalid", ngp_gid),
            ("geometry", None),
        ]
    )
    ho_ngp_list.append(ngp)

    if (
        hydrovak[["AVVBODDR", "AVVBODH"]].empty
        | hydrovak[["AVVBODDR", "AVVBODH"]].isna().values.any()
    ):
        continue

    if (
        hydrovak[["IWS_W_WATP"]].empty
        | hydrovak[["AVVTALUL"]].empty
        | hydrovak[["AVVTALUR"]].empty
    ):
        if hydrovak["AVVBODDR"].to_numpy()[0] == 0:
            bbreedte = hydrovak["IWS_W_WATB"].to_numpy()[0]
        else:
            bbreedte = hydrovak["AVVBODDR"].to_numpy()[0]

        ## TODO: een dict per parameter...
        ngp_values = dict(
            [
                ("normgeparamprofielid", ngp_gid),
                ("typeruwheid", roughness_mapping[int(branch["ruwheidsty"]) - 1]),
                ("ruwheidhoog", branch["ruwheidhoo"]),
                ("ruwheidlaag", branch["ruwheidlaa"]),
                ("bodembreedte", bbreedte),
                ("bodemhoogte benedenstrooms", hydrovak["AVVBODH"].to_numpy()[0]),
                ("bodemhoogte bovenstrooms", hydrovak["AVVBODH"].to_numpy()[0]),
                ("hoogte insteek linkerzijde", np.nan),
                ("hoogte insteek rechterzijde", np.nan),
                ("taludhelling linkerzijde", np.nan),
                ("taludhelling rechterzijde", np.nan),
                ("geometry", None),
            ]
        )
    else:
        if hydrovak["AVVBODDR"].to_numpy()[0] == 0:
            bbreedte = hydrovak["IWS_W_WATB"].to_numpy()[0]
        else:
            bbreedte = hydrovak["AVVBODDR"].to_numpy()[0]

        ngp_values = dict(
            [
                ("normgeparamprofielid", ngp_gid),
                ("typeruwheid", roughness_mapping[int(branch["ruwheidsty"]) - 1]),
                ("ruwheidhoog", branch["ruwheidhoo"]),
                ("ruwheidlaag", branch["ruwheidlaa"]),
                ("bodembreedte", bbreedte),
                ("bodemhoogte benedenstrooms", hydrovak["AVVBODH"].to_numpy()[0]),
                ("bodemhoogte bovenstrooms", hydrovak["AVVBODH"].to_numpy()[0]),
                ("hoogte insteek linkerzijde", hydrovak["IWS_W_WATP"].to_numpy()[0] + dh),
                ("hoogte insteek rechterzijde", hydrovak["IWS_W_WATP"].to_numpy()[0] + dh),
                ("taludhelling linkerzijde", hydrovak["AVVTALUL"].to_numpy()[0]),
                ("taludhelling rechterzijde", hydrovak["AVVTALUR"].to_numpy()[0]),
                ("geometry", None),
            ]
        )

    ngp_list.append(ngp_values)


hydroobject_normgp = gpd.GeoDataFrame(ho_ngp_list, geometry="geometry", crs=28992)
normgeparamprofielwaarde = gpd.GeoDataFrame(ngp_list, geometry="geometry", crs=28992)

layers = dict(
    [
        ("hydroobject", branches_gdf_new),
        ("hydroobject_normgp", hydroobject_normgp),
        ("normgeparamprofielwaarde", normgeparamprofielwaarde),
    ]
)

for name, layer in layers.items():
    layer.to_file(filename=norm_profiles_output_path, driver="GPKG", layer=name)
