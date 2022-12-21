import uuid
from typing import Tuple

import geopandas as gpd
import numpy as np
from tqdm import tqdm

# constants
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
ROUGHNESS_MAPPING = list(roughness_mapping)


def convert_pp_to_hydamo(
    branches_gdf: gpd.GeoDataFrame, index_mapping: dict
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:

    # script
    ho_ngp_list = []
    ngp_list = []

    # for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
    for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
        # Create unique ident for branch and add to branch
        branch_gid = str(uuid.uuid4())
        branches_gdf.loc[ix_1, "globalid"] = branch_gid
        branches_gdf.loc[ix_1, "code"] = branch_gid

        # Check if width and depth parameters are available, if not, skip
        # Also skip if no width is available
        if (
            (
                branch[
                    [index_mapping["bodembreedte"], index_mapping["bodemhoogte benedenstrooms"]]
                ].empty
            )
            or (
                branch[
                    [index_mapping["bodembreedte"], index_mapping["bodemhoogte benedenstrooms"]]
                ]
                .isna()
                .values.any()
            )
            or (branch[[index_mapping["bodembreedte"], index_mapping["water_width_index"]]].empty)
            or (
                branch[[index_mapping["bodembreedte"], index_mapping["water_width_index"]]]
                .isna()
                .values.any()
            )
        ):
            # print("no profile for branch {}".format(branch_gid))
            continue

        # add entry to hydroobject_normgp tabel
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

        # Check if all parameters are present for trapezium profile. If not, use rectangular profile
        width_ix = index_mapping["bodembreedte"]
        if branch[index_mapping["bodembreedte"]] == 0:
            prof_type = "rectangle"
            width_ix = index_mapping["water_width_index"]
        elif (
            np.isnan(branch[index_mapping["hoogte insteek linkerzijde"]])
            | np.isnan(branch[index_mapping["hoogte insteek rechterzijde"]])
            | np.isnan(branch[index_mapping["taludhelling linkerzijde"]])
            | np.isnan(branch[index_mapping["taludhelling rechterzijde"]])
        ):
            prof_type = "rectangle"
        elif (
            (branch[index_mapping["hoogte insteek linkerzijde"]] == 0)
            | (branch[index_mapping["hoogte insteek rechterzijde"]] == 0)
            | (branch[index_mapping["taludhelling linkerzijde"]] == 0)
            | (branch[index_mapping["taludhelling rechterzijde"]] == 0)
        ):
            prof_type = "rectangle"
        else:
            prof_type = "trapezium"

        # set required parameters for either rectangle or trapezium profile
        if prof_type == "rectangle":
            params = dict(
                [
                    ("bodembreedte", branch[width_ix]),
                    (
                        "bodemhoogte benedenstrooms",
                        branch[index_mapping["bodemhoogte benedenstrooms"]],
                    ),
                    (
                        "bodemhoogte bovenstrooms",
                        branch[index_mapping["bodemhoogte bovenstrooms"]],
                    ),
                ]
            )
        elif prof_type == "trapezium":
            params = dict(
                [
                    ("bodembreedte", branch[width_ix]),
                    (
                        "bodemhoogte benedenstrooms",
                        branch[index_mapping["bodemhoogte benedenstrooms"]],
                    ),
                    (
                        "bodemhoogte bovenstrooms",
                        branch[index_mapping["bodemhoogte bovenstrooms"]],
                    ),
                    (
                        "hoogte insteek linkerzijde",
                        branch[index_mapping["hoogte insteek linkerzijde"]],
                    ),
                    (
                        "hoogte insteek rechterzijde",
                        branch[index_mapping["hoogte insteek rechterzijde"]],
                    ),
                    (
                        "taludhelling linkerzijde",
                        branch[index_mapping["taludhelling linkerzijde"]],
                    ),
                    (
                        "taludhelling rechterzijde",
                        branch[index_mapping["taludhelling rechterzijde"]],
                    ),
                ]
            )

            # if hoogte insteek is lower than the bottom, swap bottom and hoogte insteek
            if np.amin(
                [params["hoogte insteek linkerzijde"], params["hoogte insteek rechterzijde"]]
            ) < np.amax(
                [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
            ):
                (
                    params["bodemhoogte benedenstrooms"],
                    params["bodemhoogte bovenstrooms"],
                    params["hoogte insteek linkerzijde"],
                    params["hoogte insteek rechterzijde"],
                ) = (
                    params["hoogte insteek linkerzijde"],
                    params["hoogte insteek rechterzijde"],
                    params["bodemhoogte benedenstrooms"],
                    params["bodemhoogte bovenstrooms"],
                )

        # loop over parameters to add to ngp_list
        for ix_2, (key, value) in enumerate(params.items()):
            ngp_values = dict(
                [
                    ("normgeparamprofielid", ngp_gid),
                    (
                        "typeruwheid",
                        ROUGHNESS_MAPPING[int(branch[index_mapping["typeruwheid"]]) - 1],
                    ),
                    ("ruwheidhoog", branch[index_mapping["ruwheidhoog"]]),
                    ("ruwheidlaag", branch[index_mapping["ruwheidlaag"]]),
                    ("soortparameter", key),
                    ("waarde", value),
                    ("geometry", None),
                ]
            )
            ngp_list.append(ngp_values)

    hydroobject_normgp = gpd.GeoDataFrame(ho_ngp_list, geometry="geometry", crs=28992)
    normgeparamprofielwaarde = gpd.GeoDataFrame(ngp_list, geometry="geometry", crs=28992)

    return (
        branches_gdf,
        hydroobject_normgp,
        normgeparamprofielwaarde,
    )
