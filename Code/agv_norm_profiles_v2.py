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
index_mapping = dict(
    [
        ("bodembreedte", "AVVBODDR"),
        ("bodemhoogte benedenstrooms", "AVVBODH"),
        ("bodemhoogte bovenstrooms", "AVVBODH"),
        ("hoogte insteek linkerzijde", "IWS_W_WATP"),
        ("hoogte insteek rechterzijde", "IWS_W_WATP"),
        ("taludhelling linkerzijde", "AVVTALUL"),
        ("taludhelling rechterzijde", "AVVTALUR"),
    ]
)

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

    if (
        hydrovak[["AVVBODDR", "AVVBODH"]].empty
        | hydrovak[["AVVBODDR", "AVVBODH"]].isna().values.any()
    ):
        continue

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
        hydrovak[["IWS_W_WATP"]].empty
        | hydrovak[["AVVTALUL"]].empty
        | hydrovak[["AVVTALUR"]].empty
    ):
        prof_type = "rectangle"
    else:
        prof_type = "trapezium"

    for ix, (key, value) in enumerate(index_mapping.items()):
        # if bodembreedte == 0, choose waterwidth
        if (ix == 0) & (hydrovak[value].to_numpy()[0] == 0):
            _value = "IWS_W_WATB"
        else:
            _value = value

        # skip trapezium parameters if profile is rectangle
        if (ix == 3) & (prof_type == "rectangle"):
            break

        ngp_values = dict(
            [
                ("normgeparamprofielid", ngp_gid),
                ("typeruwheid", roughness_mapping[int(branch["ruwheidsty"]) - 1]),
                ("ruwheidhoog", branch["ruwheidhoo"]),
                ("ruwheidlaag", branch["ruwheidlaa"]),
                ("soortparameter", key),
                ("waarde", hydrovak[_value].to_numpy()[0]),
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
