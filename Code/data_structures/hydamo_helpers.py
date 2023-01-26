import importlib
import uuid
from copy import copy
from typing import List, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, MultiPoint, Point

from data_structures.dhydamo_data_model import DHydamoDataModel
from data_structures.hydamo_globals import (
    MANAGEMENT_DEVICE_TYPES,
    ROUGHNESS_MAPPING_LIST,
    WEIR_MAPPING,
)


def check_column_is_numerical(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    if (not isinstance(gdf, int)) | (not isinstance(gdf, float)):
        gdf = gdf.astype(float)
    return gdf


def check_roughness(structure: gpd.GeoSeries, rougness_map: List = ROUGHNESS_MAPPING_LIST):
    """ """
    type_ruwheid = structure["typeruwheid"]
    if isinstance(type_ruwheid, int) or isinstance(type_ruwheid, float):
        type_ruwheid = rougness_map[int(type_ruwheid) - 1]
    return type_ruwheid


def create_bridge_data(bridge_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """ """
    _bridge_gdf = copy(bridge_gdf)
    for ix, bridge in bridge_gdf.iterrows():
        # turn numerical roughnes types to strings
        _bridge_gdf.loc[ix, "typeruwheid"] = check_roughness(bridge)
    return _bridge_gdf


def create_culvert_data(culvert_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """ """
    culvert_gdf["doorstroomopening"] = None

    culvert_gdf["breedteopening"] = check_column_is_numerical(gdf=culvert_gdf["breedteopening"])
    culvert_gdf["hoogteopening"] = check_column_is_numerical(gdf=culvert_gdf["hoogteopening"])

    culvert_gdf["breedteopening"] = culvert_gdf["breedteopening"].replace(0, np.nan)
    culvert_gdf["hoogteopening"] = culvert_gdf["hoogteopening"].replace(0, np.nan)
    _culvert_gdf = copy(culvert_gdf)

    for ix, culvert in culvert_gdf.iterrows():
        if np.isnan(culvert["breedteopening"]) or np.isnan(culvert["hoogteopening"]):
            shape = "circle"

        # elif (culvert["vormkoker"].dtype == "int64") or (culvert["vormkoker"].dtype == "float64"):
        elif isinstance(culvert["vormkoker"], int) | isinstance(culvert["vormkoker"], float):
            if int(culvert["vormkoker"]) == 3:
                shape = "rectangle"
            else:
                shape = "circle"

        elif isinstance(culvert["vormkoker"], str):  # dtype not a number, so assuming string
            if str(culvert["vormkoker"]).lower() == "rechthoekig":
                shape = "rectangle"
                _culvert_gdf.loc[ix, "vormkoker"] = 3
            else:
                shape = "circle"
                _culvert_gdf.loc[ix, "vormkoker"] = 1
        else:
            raise ValueError("wrong datatype: {}".format(type(culvert["vormkoker"])))

        if shape == "rectangle":
            crosssection = {
                "shape": shape,
                "height": culvert["hoogteopening"],
                "width": culvert["breedteopening"],
                "closed": culvert["gesloten"],
            }
        elif shape == "circle":
            if not np.isnan(culvert["breedteopening"]):
                diameter = culvert["breedteopening"]
            elif not np.isnan(culvert["breedteopening"]):
                diameter = culvert["breedteopening"]
            else:
                raise ValueError("onbekende koker diameter")

            crosssection = {
                "shape": shape,
                "diameter": diameter,
            }
        _culvert_gdf.loc[ix, "doorstroomopening"] = str(crosssection)

        # turn numerical roughnes types to strings
        _culvert_gdf.loc[ix, "typeruwheid"] = check_roughness(culvert)

    _culvert_gdf = _culvert_gdf.drop(columns="gesloten")
    return _culvert_gdf


def create_measured_profile_data(
    profile_points_gdf: gpd.GeoDataFrame,
    dist_tol: float = 0.25,
    roughness_mapping: List = ROUGHNESS_MAPPING_LIST,
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:

    # Load shapefile with profile points & group them by metingprof attribute
    # profile_points = gpd.read_file(profile_points_path)
    grouped_points = profile_points_gdf.groupby(by="profiel nummer")

    # initialize empty lists
    profile_groups = []
    profile_lines = []
    profile_points = []
    roughnes_profiles = []

    # Loop over the grouped profile points to create a line from the points
    count = 0
    for name, group in grouped_points:
        # sort points based on codevolg number
        sorted_group = group.sort_values("codevolgnummer")

        # if both bottom and sludge measurements are present, select the first
        if (
            sorted_group["type meting"].isin([1]).any()
            and sorted_group["type meting"].isin([2]).any()
        ):
            sorted_group.drop(sorted_group[sorted_group["type meting"] == 2].index, inplace=True)

        # skip 'line' if it only has one point
        if sorted_group.shape[0] < 5:
            continue

        # Initialize name and list of points
        profile_group_id = str(uuid.uuid4())
        profile_line_id = str(uuid.uuid4())
        #

        list_of_points = []

        # Add info to profile_group table
        _profile_group = dict(
            [
                ("globalid", profile_group_id),
                ("brugid", None),
                ("stuwid", None),
                ("geometry", None),
            ]
        )
        # _profile_group = dict([("globalid", profile_line_id), ("geometry", None)])
        profile_groups.append(_profile_group)

        code_volg_nr = 0
        # add points to line
        for ix, row in sorted_group.iterrows():
            if type(row.geometry) == Point:
                l_points = [row.geometry]
            elif type(row.geometry) == MultiPoint:
                l_points = row.geometry.geoms

            for point in l_points:
                # Check if profile points are too close together, and skip if that is the case
                # if ix < sorted_group.shape[0]:
                if code_volg_nr > 0:
                    p_0 = list_of_points[code_volg_nr - 1]
                    p_1 = point
                    p_dist = p_0.distance(p_1)

                    if p_dist < dist_tol:
                        continue

                # Because sorted_group has been sorted on codevolgnr, we can assume sequentiallity
                code_volg_nr += 1

                # append point to list for line generation
                list_of_points.append(point)

                # create entry for point shape and append
                # point_id = "AGV_" + str(row["code"])
                point_code = profile_line_id + r"_" + str(code_volg_nr)
                point_id = str(uuid.uuid4())
                _point = dict(
                    [
                        ("code", point_code),
                        ("globalid", point_id),
                        ("profiellijnid", profile_line_id),
                        # ("codevolgnummer", row["codevolgnu"]),
                        ("codevolgnummer", code_volg_nr),
                        ("geometry", point),
                    ]
                )
                profile_points.append(_point)

                # create roughness table entry
                _roughness_profile = dict(
                    [
                        ("code", point_code),
                        ("profielpuntid", point_id),
                        ("typeruwheid", roughness_mapping[int(row["typeruwheid"]) - 1]),
                        ("ruwheidhoog", float(row["ruwheidhoog"])),
                        ("ruwheidlaag", float(row["ruwheidlaag"])),
                        ("geometry", None),
                    ]
                )
                roughnes_profiles.append(_roughness_profile)

        # Convert points to line
        profile_line = LineString(list_of_points)

        # add line to list
        _profile_line = dict(
            [
                ("globalid", profile_line_id),
                ("profielgroepid", profile_group_id),
                ("geometry", profile_line),
            ]
        )
        profile_lines.append(_profile_line)

    # Create GeoDataFrame from list of dicts
    profile_groups = gpd.GeoDataFrame(profile_groups, geometry="geometry", crs=28992)
    profile_lines = gpd.GeoDataFrame(profile_lines, geometry="geometry", crs=28992)
    profile_points = gpd.GeoDataFrame(profile_points, geometry="geometry", crs=28992)
    roughnes_profiles = gpd.GeoDataFrame(roughnes_profiles, geometry="geometry", crs=28992)

    return profile_groups, profile_lines, profile_points, roughnes_profiles


def create_norm_parm_profiles_v2(
    branches_gdf: gpd.GeoDataFrame, min_water_width: float = 0.1
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Function that converts parameterized profiles to HYDAMO compliant format.

    Args:
        branches_gdf (gpd.GeoDataFrame): input geodataframe containing branches
        index_mapping (dict): dictionary containing a mapping from required keys to values present in branches_gdf

    Returns:
        hydroobject (gpd.GeoDataFrame): ouput geodataframe containing branches
        hydroobject_normgp (gpd.GeoDataFrame): output geodataframe containing names of branches and corresponding profiles
        normgeparamprofielwaarde (gpd.GeoDataFrame): output geodataframe containing parameterized profiles
    """

    # script
    ho_ngp_list = []
    ngp_list = []

    branches_gdf.loc[
        branches_gdf["bodembreedte"] < min_water_width,
        "bodembreedte",
    ] = np.nan
    branches_gdf.loc[
        branches_gdf["water_width_index"] < min_water_width,
        "water_width_index",
    ] = np.nan

    branches_out_gdf = copy(branches_gdf)

    # Check for duplicate codes, if there are, they are replaced in the following loop
    duplicates = branches_gdf.duplicated(subset=["code"])

    # for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
    for jx, (ix_1, branch) in enumerate(branches_gdf.iterrows()):
        # Create unique ident for branch and add to branch
        branch_gid = str(uuid.uuid4())
        branches_out_gdf.loc[ix_1, "globalid"] = branch_gid

        # replace duplicate codes by addinng a uniqure number
        if duplicates[jx]:
            branches_out_gdf.loc[ix_1, "code"] = branches_out_gdf.loc[ix_1, "code"] + str(jx)
        # branches_out_gdf.loc[ix_1, "code"] = branch_gid

        # turn numerical roughnes types to strings
        type_ruwheid = check_roughness(branch)
        branches_out_gdf.loc[ix_1, "typeruwheid"] = type_ruwheid

        # check if linestring
        assert isinstance(branch.geometry, LineString)

        # Delete z-dimensison of linestring if it exists
        if np.array(branch.geometry.coords).shape[1] > 2:
            branches_out_gdf.loc[ix_1, "geometry"] = LineString(
                [xy[0:2] for xy in list(branch.geometry.coords)]
            )

        # Check if width and depth parameters are available, if not, skip
        # Also skip if no width is available
        if (
            (branch[["bodembreedte", "bodemhoogte benedenstrooms"]].empty)
            or (branch[["bodembreedte", "bodemhoogte benedenstrooms"]].isna().values.any())
            or (branch[["bodembreedte", "water_width_index"]].empty)
            or (branch[["bodembreedte", "water_width_index"]].empty)
            or (branch[["bodembreedte", "water_width_index"]].isna().values.any())
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
        width_ix = "bodembreedte"
        if np.isnan(branch["bodembreedte"]):
            prof_type = "rectangle"
            width_ix = "water_width_index"
        elif (
            np.isnan(branch["hoogte insteek linkerzijde"])
            | np.isnan(branch["hoogte insteek rechterzijde"])
            | np.isnan(branch["taludhelling linkerzijde"])
            | np.isnan(branch["taludhelling rechterzijde"])
        ):
            prof_type = "rectangle"
        elif (
            (branch["hoogte insteek linkerzijde"] == 0)
            | (branch["hoogte insteek rechterzijde"] == 0)
            | (branch["taludhelling linkerzijde"] == 0)
            | (branch["taludhelling rechterzijde"] == 0)
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
                        branch["bodemhoogte benedenstrooms"],
                    ),
                    (
                        "bodemhoogte bovenstrooms",
                        branch["bodemhoogte bovenstrooms"],
                    ),
                ]
            )
        elif prof_type == "trapezium":
            params = dict(
                [
                    ("bodembreedte", branch[width_ix]),
                    (
                        "bodemhoogte benedenstrooms",
                        branch["bodemhoogte benedenstrooms"],
                    ),
                    (
                        "bodemhoogte bovenstrooms",
                        branch["bodemhoogte bovenstrooms"],
                    ),
                    (
                        "hoogte insteek linkerzijde",
                        branch["hoogte insteek linkerzijde"],
                    ),
                    (
                        "hoogte insteek rechterzijde",
                        branch["hoogte insteek rechterzijde"],
                    ),
                    (
                        "taludhelling linkerzijde",
                        branch["taludhelling linkerzijde"],
                    ),
                    (
                        "taludhelling rechterzijde",
                        branch["taludhelling rechterzijde"],
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
                        type_ruwheid,
                    ),
                    ("ruwheidhoog", branch["ruwheidhoog"]),
                    ("ruwheidlaag", branch["ruwheidlaag"]),
                    ("soortparameter", key),
                    ("waarde", value),
                    ("geometry", None),
                ]
            )
            ngp_list.append(ngp_values)

    hydroobject_normgp = gpd.GeoDataFrame(ho_ngp_list, geometry="geometry", crs=28992)
    normgeparamprofielwaarde = gpd.GeoDataFrame(ngp_list, geometry="geometry", crs=28992)

    return (
        branches_out_gdf[["code", "globalid", "geometry", "typeruwheid"]],
        hydroobject_normgp,
        normgeparamprofielwaarde,
    )


def create_norm_parm_profiles(
    branches_gdf: gpd.GeoDataFrame, index_mapping: dict, min_water_width: float = 0.1
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Function that converts parameterized profiles to HYDAMO compliant format.

    Args:
        branches_gdf (gpd.GeoDataFrame): input geodataframe containing branches
        index_mapping (dict): dictionary containing a mapping from required keys to values present in branches_gdf

    Returns:
        hydroobject (gpd.GeoDataFrame): ouput geodataframe containing branches
        hydroobject_normgp (gpd.GeoDataFrame): output geodataframe containing names of branches and corresponding profiles
        normgeparamprofielwaarde (gpd.GeoDataFrame): output geodataframe containing parameterized profiles
    """

    # script
    ho_ngp_list = []
    ngp_list = []

    branches_gdf.loc[
        branches_gdf[index_mapping["bodembreedte"]] < min_water_width,
        index_mapping["bodembreedte"],
    ] = np.nan
    branches_gdf.loc[
        branches_gdf[index_mapping["water_width_index"]] < min_water_width,
        index_mapping["water_width_index"],
    ] = np.nan

    branches_out_gdf = copy(branches_gdf)

    # for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
    for ix_1, branch in branches_gdf.iterrows():
        # Create unique ident for branch and add to branch
        branch_gid = str(uuid.uuid4())
        branches_out_gdf.loc[ix_1, "globalid"] = branch_gid
        branches_out_gdf.loc[ix_1, "code"] = branch_gid

        # turn numerical roughnes types to strings
        type_ruwheid = check_roughness(branch)
        branches_out_gdf.loc[ix_1, "typeruwheid"] = type_ruwheid

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
        if np.isnan(branch[index_mapping["bodembreedte"]]):
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
                        type_ruwheid,
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
        branches_out_gdf[["code", "globalid", "geometry", "typeruwheid"]],
        hydroobject_normgp,
        normgeparamprofielwaarde,
    )


def create_pump_data(pump_gdf: gpd.GeoDataFrame) -> List[gpd.GeoDataFrame]:

    pump_station_list = []
    pump_list = []
    management_list = []

    # if not isinstance(pump_gdf["maximalecapaciteit"], float):
    #     pump_gdf["maximalecapaciteit"] = pump_gdf["maximalecapaciteit"].astype(float)
    pump_gdf["maximalecapaciteit"] = check_column_is_numerical(gdf=pump_gdf["maximalecapaciteit"])

    for ix, pump in pump_gdf.iterrows():
        pump_station = dict(
            [
                ("code", pump["code"]),
                ("geometry", pump["geometry"]),
                ("globalid", pump["globalid"]),
            ]
        )
        pump_station_list.append(pump_station)

        pump_gid = str(uuid.uuid4())
        _pump = dict(
            [
                ("code", pump["code"]),
                ("geometry", None),
                ("gemaalid", pump["globalid"]),
                ("globalid", pump_gid),
                ("maximalecapaciteit", pump["maximalecapaciteit"]),
            ]
        )
        pump_list.append(_pump)

        management = dict(
            [
                ("bovengrens", pump["streefwaarde"] + 0.5 * pump["peil_marge"]),
                ("code", pump["code"]),
                ("doelvariabele", pump["doelvariabele"]),
                ("geometry", None),
                ("globalid", str(uuid.uuid4())),
                ("ondergrens", pump["streefwaarde"] - 0.5 * pump["peil_marge"]),
                ("pompid", pump_gid),
                ("streefwaarde", pump["streefwaarde"]),
            ]
        )
        management_list.append(management)

    pump_station_gdf = gpd.GeoDataFrame(pump_station_list, geometry="geometry", crs=28992)
    _pump_gdf = gpd.GeoDataFrame(pump_list, geometry="geometry", crs=28992)
    management_gdf = gpd.GeoDataFrame(management_list, geometry="geometry", crs=28992)

    return pump_station_gdf, _pump_gdf, management_gdf


def create_river_profiles(riv_prof_df: pd.DataFrame) -> gpd.GeoDataFrame:
    ## TODO build list of dicts
    riv_prof_df["id"] = riv_prof_df["id"].astype(str)
    unique_ids = riv_prof_df["id"].unique()
    geom_df = pd.DataFrame(
        columns=[
            "name",
            "ix",
            "levels",
            "flowWidths",
            "totalWidths",
            "geometry",
        ]
    )
    meta_df = pd.DataFrame(
        columns=[
            "name",
            "thalweg",
            "leveecrestLevel",
            "leveebaselevel",
            "leveeflowarea",
            "leveetotalarea",
            "mainwidth",
            "fp1width",
            "fp2width",
            "branchid",
            "chainage",
            "geometry",
        ]
    )

    meta_df["chainage"] = meta_df["chainage"].astype(float)
    geom_df["levels"] = geom_df["levels"].astype(float)
    geom_df["flowWidths"] = geom_df["flowWidths"].astype(float)
    geom_df["totalWidths"] = geom_df["totalWidths"].astype(float)

    g_ix = 0
    for u_id in unique_ids:
        slice = riv_prof_df[riv_prof_df["id"] == u_id]
        meta_data = slice[slice["Data_type"] == "meta"]
        geom_data = slice[slice["Data_type"] == "geom"]
        # print(meta_data.dropna(axis=1))
        # print(geom_data.dropna(axis=1))
        for ix in range(geom_data.shape[0]):
            geom_df.at[g_ix, "name"] = u_id
            geom_df.at[g_ix, "ix"] = ix
            geom_df.at[g_ix, "levels"] = geom_data["level"].values[ix]
            geom_df.at[g_ix, "flowWidths"] = geom_data["Flow width"].values[ix]
            geom_df.at[g_ix, "totalWidths"] = geom_data["Total width"].values[ix]
            geom_df.at[g_ix, "geometry"] = None
            g_ix += 1

        meta_df.at[u_id, "numLevels"] = geom_data.shape[0]
        meta_df.at[u_id, "name"] = u_id
        meta_df.at[u_id, "thalweg"] = 0
        meta_df.at[u_id, "leveecrestLevel"] = meta_data["Crest level summerdike"].values[0]

        meta_df.at[u_id, "leveebaselevel"] = meta_data[
            "Floodplain baselevel behind summerdike"
        ].values[0]
        meta_df.at[u_id, "leveeflowarea"] = meta_data["Flow area behind summerdike"].values[0]
        meta_df.at[u_id, "leveetotalarea"] = meta_data["Total area behind summerdike"].values[0]
        meta_df.at[u_id, "mainwidth"] = meta_data["width main channel"].values[0]
        meta_df.at[u_id, "fp1width"] = meta_data["width floodplain 1"].values[0]
        meta_df.at[u_id, "fp2width"] = meta_data["width floodplain 2"].values[0]
        meta_df.at[u_id, "branchid"] = meta_data["branch"].values[0]
        meta_df.at[u_id, "chainage"] = meta_data["chainage"].values[0]
        meta_df.at[u_id, "geometry"] = None

        total_width = (
            meta_df.at[u_id, "mainwidth"]
            + meta_df.at[u_id, "fp1width"]
            + meta_df.at[u_id, "fp2width"]
        )
        max_width = geom_df.loc[geom_df["name"] == u_id, "flowWidths"].max()
        if total_width > max_width:
            diff = total_width - max_width

            max_ix = geom_df.loc[geom_df["name"] == u_id, "flowWidths"].idxmax()
            geom_df.at[max_ix, "flowWidths"] = geom_df.at[max_ix, "flowWidths"] + diff
            if geom_df.at[max_ix, "flowWidths"] > geom_df.at[max_ix, "totalWidths"]:
                geom_df.at[max_ix, "totalWidths"] = geom_df.at[max_ix, "flowWidths"]

    geom_gdf = gpd.GeoDataFrame(geom_df, geometry="geometry", crs="epsg:28992")
    meta_gdf = gpd.GeoDataFrame(meta_df, geometry="geometry", crs="epsg:28992")
    return geom_gdf, meta_gdf


def create_weir_data(weir_gdf: gpd.GeoDataFrame) -> List[gpd.GeoDataFrame]:

    weir_list = []
    opening_list = []
    management_device_list = []

    # if not "hoogstedoorstroombreedte" in weir_gdf.columns:
    #     weir_gdf["hoogstedoorstroombreedte"] = weir_gdf["laagstedoorstroombreedte"]

    # if not "hoogstedoorstroomhoogte" in weir_gdf.columns:
    #     weir_gdf["hoogstedoorstroomhoogte"] = weir_gdf["laagstedoorstroomhoogte"] + 10
    clist = [
        "hoogstedoorstroombreedte",
        "hoogstedoorstroomhoogte",
        "laagstedoorstroombreedte",
        "laagstedoorstroomhoogte",
    ]
    for column in clist:
        weir_gdf[column] = check_column_is_numerical(gdf=weir_gdf[column])

    weir_gdf["hoogstedoorstroombreedte"].fillna(weir_gdf["laagstedoorstroombreedte"], inplace=True)
    weir_gdf["hoogstedoorstroomhoogte"].fillna(
        weir_gdf["laagstedoorstroomhoogte"] + 10, inplace=True
    )
    for ix, weir in weir_gdf.iterrows():
        if isinstance(weir["soortstuw"], str):
            try:
                soort_stuw = WEIR_MAPPING[weir["soortstuw"].lower()]
            except KeyError as e:
                print(e)
                print("choosing overlaat")
                soort_stuw = 11

        else:
            soort_stuw = weir["soortstuw"]
        _weir = dict(
            [
                ("afvoercoefficient", weir["afvoercoefficient_stuw"]),
                ("code", weir["code"]),
                ("geometry", weir["geometry"]),
                ("globalid", weir["globalid"]),
                ("soortstuw", soort_stuw),
            ]
        )
        weir_list.append(_weir)

        if isinstance(weir["vormopening"], int) | isinstance(weir["vormopening"], float):
            shape = weir["vormopening"]

        elif isinstance(weir["vormopening"], str):  # dtype not a number, so assuming string
            if str(weir["vormopening"]).lower() == "rechthoekig":
                shape = 3
            else:
                shape = 1

        opening_gid = str(uuid.uuid4())
        opening = dict(
            [
                ("afvoercoefficient", weir["afvoercoefficient_opening"]),
                ("geometry", None),
                ("globalid", opening_gid),
                ("hoogstedoorstroombreedte", weir["hoogstedoorstroombreedte"]),
                ("hoogstedoorstroomhoogte", weir["hoogstedoorstroomhoogte"]),
                ("laagstedoorstroombreedte", weir["laagstedoorstroombreedte"]),
                ("laagstedoorstroomhoogte", weir["laagstedoorstroomhoogte"]),
                ("stuwid", weir["globalid"]),
                ("vormopening", shape),
            ]
        )
        opening_list.append(opening)

        if isinstance(weir["soortregelbaarheid"], int) | isinstance(
            weir["soortregelbaarheid"], float
        ):
            management = weir["soortregelbaarheid"]

        elif isinstance(weir["soortregelbaarheid"], str):  # dtype not a number, so assuming string
            try:
                management = MANAGEMENT_DEVICE_TYPES[weir["soortregelbaarheid"]]
            except KeyError as e:
                print(e)
                print("choosing niet regelbaar")
                management = 1

        management_device = dict(
            [
                ("code", weir["code"]),
                ("geometry", weir["geometry"]),
                ("globalid", weir["globalid"]),
                ("kunstwerkopeningid", opening_gid),
                ("overlaatonderlaat", weir["overlaatonderlaat"]),
                ("soortregelbaarheid", management),
                ("stuwid", weir["globalid"]),
            ]
        )
        management_device_list.append(management_device)

    _weir_gdf = gpd.GeoDataFrame(weir_list, geometry="geometry", crs=28992)
    opening_gdf = gpd.GeoDataFrame(opening_list, geometry="geometry", crs=28992)
    management_device_gdf = gpd.GeoDataFrame(
        management_device_list, geometry="geometry", crs=28992
    )

    return _weir_gdf, opening_gdf, management_device_gdf


def fill_empty_columns(
    defaults, gdf: gpd.GeoDataFrame, index_mapping: dict
) -> Tuple[gpd.GeoDataFrame, dict]:

    _index_mapping = copy(index_mapping)
    _gdf = copy(gdf)

    # loop over column names in index mapping, with DHYDAMO name in key and column name in shapefile in Value
    for key, value in index_mapping.items():
        _key = key.replace(" ", "_")
        # if data is not in shapefile, fill with default values
        if value is None:
            _gdf[key] = getattr(defaults, _key)
            _index_mapping[key] = key

        else:
            _gdf[key] = gdf[value]

            # if it is in shapefile, but contains missing values, fill them as well
            if _gdf[key].dtype == "geometry":
                pass

            elif (_gdf[key].dtype == "int64") or (_gdf[key].dtype == "float64"):
                if hasattr(defaults, _key):
                    _gdf[key] = _gdf[key].replace(
                        to_replace=[0, -999, np.nan],
                        value=getattr(defaults, _key),
                    )
                # print number of 0, -998 and nan if no default value is present
                else:
                    n_zero = _gdf[key][_gdf[key] == 0].count()
                    n_nnn = _gdf[key][_gdf[key] == -999].count()
                    n_nan = _gdf[key].isna().sum()

                    print("{}: entries with missing data: {}".format(key, n_zero + n_nnn + n_nan))

            elif (_gdf[key].dtype == "string") or (_gdf[key].dtype == "object"):
                if hasattr(defaults, _key):
                    _gdf[key] = _gdf[key].replace(
                        to_replace=[None],
                        value=getattr(defaults, _key),
                    )
                else:
                    n_nan = _gdf[key].isna().sum()
                    print("{}: entries with missing data: {}".format(key, n_nan))
            else:
                raise ValueError("type not yet implemented: {}".format(_gdf[key].dtype))

    return _gdf, _index_mapping


def fix_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """drops values without geometry and sets corrects crs"""
    gdf.dropna(axis=0, inplace=True, subset=["geometry"])
    try:
        gdf = gdf.to_crs("epsg:28992")
    except ValueError:
        gdf = gdf.set_crs("epsg:28992").to_crs("epsg:28992")
    return gdf


def load_geo_file(file_path: str, layer: str = None):
    """ """
    if file_path.endswith(r".shp"):
        gdf = gpd.read_file(file_path)
    elif file_path.endswith(r".gpkg"):
        if layer is not None:
            gdf = gpd.read_file(file_path, layer=layer)
        else:
            raise ValueError("provide a layer when loading a gpkg")
    else:
        raise ValueError("filetype not implemented")

    mls_struct_bool = (gdf.geometry.type == "MultiLineString") | (
        gdf.geometry.type == "MultiPoint"
    )

    if np.sum(mls_struct_bool) > 0:
        gdf = gdf.explode(ignore_index=True, index_parts=False)

    return gdf


def map_columns(defaults, gdf: gpd.GeoDataFrame, index_mapping: dict) -> gpd.GeoDataFrame:
    """ """
    # ensure no incorrect geometries exist and set a globalid
    gdf = fix_geometry(gdf)
    gdf["globalid"] = [str(uuid.uuid4()) for _ in range(gdf.shape[0])]

    # fill potential empty columns
    gdf, index_mapping = fill_empty_columns(
        defaults=defaults,
        gdf=gdf,
        index_mapping=index_mapping,
    )

    # return only columns that are in that dict
    return gdf[list(index_mapping.keys())], index_mapping


def convert_to_dhydamo_data(defaults: str, config: str) -> DHydamoDataModel:
    """ """
    ddm = DHydamoDataModel()
    defaults = importlib.import_module("dataset_configs." + defaults)
    data_config = getattr(importlib.import_module("dataset_configs." + config), "RawData")

    if hasattr(data_config, "branches_path"):
        if data_config.branches_path is not None:
            ## Branches
            branches_gdf = load_geo_file(data_config.branches_path, layer="waterloop")

            branches_gdf, index_mapping = map_columns(
                defaults=defaults.Branches,
                gdf=branches_gdf,
                index_mapping=data_config.branch_index_mapping,
            )

            (
                ddm.waterloop,
                ddm.hydroobject_normgp,
                ddm.normgeparamprofielwaarde,
            ) = create_norm_parm_profiles_v2(branches_gdf=branches_gdf)

    if hasattr(data_config, "bridges_path"):
        if data_config.bridges_path is not None:
            ## Bridges
            bridges_gdf, _ = map_columns(
                defaults=defaults.Bridges,
                gdf=load_geo_file(data_config.bridges_path, layer="brug"),
                index_mapping=data_config.bridge_index_mapping,
            )
            ddm.brug = create_bridge_data(bridge_gdf=bridges_gdf)

    if hasattr(data_config, "culvert_path"):
        if data_config.culvert_path is not None:
            ## Culverts
            culvert_gdf, _ = map_columns(
                defaults=defaults.Culverts,
                gdf=load_geo_file(data_config.culvert_path, layer="duiker"),
                index_mapping=data_config.culvert_index_mapping,
            )
            ddm.duiker = create_culvert_data(culvert_gdf=culvert_gdf)

    if hasattr(data_config, "measured_profile_path"):
        if data_config.measured_profile_path is not None:
            measure_profile_gdf = load_geo_file(
                data_config.measured_profile_path,
                layer="metingprofielpunt",
            )
            measure_profile_gdf, _ = map_columns(
                defaults=defaults.MeasuredProfiles,
                gdf=measure_profile_gdf,
                index_mapping=data_config.measured_profile_index_mapping,
            )
            (
                ddm.profielgroep,
                ddm.profiellijn,
                ddm.profielpunt,
                ddm.ruwheidsprofiel,
            ) = create_measured_profile_data(profile_points_gdf=measure_profile_gdf)

    if hasattr(data_config, "pump_path"):
        if data_config.pump_path is not None:
            ## Pumps
            pump_gdf = load_geo_file(data_config.pump_path, layer="gemaal")
            pump_gdf, _ = map_columns(
                defaults=defaults.Pumps,
                gdf=pump_gdf,
                index_mapping=data_config.pump_index_mapping,
            )

            ddm.gemaal, ddm.pomp, ddm.sturing = create_pump_data(pump_gdf=pump_gdf)

    if hasattr(data_config, "river_profile_path"):
        if data_config.river_profile_path is not None:
            riv_prof_df = pd.read_csv(data_config.river_profile_path)
            ddm.rivier_profielen, ddm.rivier_profielen_data = create_river_profiles(
                riv_prof_df=riv_prof_df
            )

    if hasattr(data_config, "weir_path"):
        if data_config.weir_path is not None:
            ## Weirs
            weir_gdf = load_geo_file(data_config.weir_path, layer="stuw")
            weir_gdf, _ = map_columns(
                defaults=defaults.Weirs, gdf=weir_gdf, index_mapping=data_config.weir_index_mapping
            )

            ddm.stuw, ddm.kunstwerkopening, ddm.regelmiddel = create_weir_data(weir_gdf=weir_gdf)

    return ddm


def save_gpkg(
    input_gdfs: List[gpd.GeoDataFrame], layers: List[str], output_path: str = None
) -> None:
    if len(input_gdfs) != len(layers):
        raise ValueError("expected lists of equal lengths")

    # save hydamo data in geopackage
    if output_path is not None:
        for name, gdf in zip(layers, input_gdfs):
            gdf.to_file(filename=output_path, driver="GPKG", layer=name)
