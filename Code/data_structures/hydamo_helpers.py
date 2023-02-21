import importlib
import uuid
import warnings
from copy import copy
from typing import List, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.geometry import LineString, MultiPoint, Point

from data_structures.dhydamo_data_model import DHydamoDataModel
from data_structures.hydamo_globals import (
    MANAGEMENT_DEVICE_TYPES,
    ROUGHNESS_MAPPING_LIST,
    WEIR_MAPPING,
)

warnings.filterwarnings(action="ignore", message="Mean of empty slice")


def check_and_fix_duplicate_code(gdf: gpd.GeoDataFrame, column="code") -> gpd.GeoDataFrame:
    if column in gdf.columns:
        gdf[column] = gdf[column].astype("str").str.strip()
        # _gdf = copy(gdf)
        duplicates = gdf.duplicated(subset=column, keep=False)
        duplicate_codes = gdf.loc[duplicates, column]

        for code in duplicate_codes.unique():
            n_duplicates = np.sum(gdf[column] == code)
            pad_list = [r"{}_{}".format(code, n) for n in np.arange(n_duplicates, dtype=np.int8)]

            gdf.loc[gdf[column] == code, column] = pad_list

    return gdf


def check_column_is_numerical(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    if (not isinstance(gdf, int)) | (not isinstance(gdf, float)):
        gdf = gdf.astype(float)
    return gdf


def check_is_not_na_number(input, zero_allowed=False) -> bool:
    if (input is None) or (np.isnan(input)):
        return False

    if (not zero_allowed) and (input == 0):
        return False

    return True


def check_roughness(structure: gpd.GeoSeries, rougness_map: List = ROUGHNESS_MAPPING_LIST) -> str:
    """ """
    type_ruwheid = structure["typeruwheid"]
    if isinstance(type_ruwheid, int) or isinstance(type_ruwheid, float):
        type_ruwheid = rougness_map[int(type_ruwheid) - 1]
    return type_ruwheid


def convert_to_dhydamo_data(defaults: str, config: str) -> DHydamoDataModel:
    """ """

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

        culvert_gdf["breedteopening"] = check_column_is_numerical(
            gdf=culvert_gdf["breedteopening"]
        )
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
                sorted_group.drop(
                    sorted_group[sorted_group["type meting"] == 2].index, inplace=True
                )

            # skip 'line' if it only has one point
            if sorted_group.shape[0] < 4:
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
                            # ("typeruwheid", roughness_mapping[int(row["typeruwheid"]) - 1]),
                            ("typeruwheid", check_roughness(row)),
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

    def create_mp_from_np(
        branches_gdf: gpd.GeoDataFrame,
        dist_tol: float = 0.25,
        min_water_width: float = 0.1,
        roughness_mapping: List = ROUGHNESS_MAPPING_LIST,
    ) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
        def rectangular_point_profile(
            branch: LineString,
            profiel_nummer: str,
            params: dict,
            def_height: float = 2,
            rect_offset: float = 0.1,
            interp_range: float = 0.1,
        ) -> List[Point]:

            bwidth = params["bodembreedte"]

            # l_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=(bwidth + rect_offset) / 2, side="left"
            # )
            # l_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="left")
            # r_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=(bwidth + rect_offset) / 2, side="right"
            # )
            # r_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="right")

            rotated_line = affinity.rotate(
                branch.geometry, 90, origin=branch.geometry.interpolate(0.5, normalized=True)
            )

            p1 = rotated_line.interpolate(0.5 - interp_range / 2, normalized=True)
            p2 = rotated_line.interpolate(0.5 + interp_range / 2, normalized=True)

            if p2.x > p1.x:
                dx_o = p2.x - p1.x
                dy_o = p2.y - p1.y
            else:
                dx_o = p1.x - p2.x
                dy_o = p1.y - p2.y

            if dx_o != 0:
                s = dy_o / dx_o
                # given a line slope dy/dx (i.e. a/b), pythagoras a**2 + b**2 = c**2, and a desired length c
                # we can determine that b_n = c_n/sqrt(1+(a/b)**2) and a_n = b_n * (a/b)
                dx_1 = bwidth / np.sqrt(1 + s**2)
                dy_1 = dx_1 * s

                dx_2 = (bwidth + rect_offset) / np.sqrt(1 + s**2)
                dy_2 = dx_2 * s
            else:
                dx_1 = dx_2 = 0
                dy_1 = dy_2 = bwidth

            centroid = branch.geometry.interpolate(0.5, normalized=True)
            # centroid = rotated_line.centroid
            c_x, c_y = centroid.x, centroid.y

            list_of_points = [
                Point(c_x - 0.5 * dx_2, c_y - 0.5 * dy_2),
                Point(c_x - 0.5 * dx_1, c_y - 0.5 * dy_1),
                Point(c_x + 0.5 * dx_1, c_y + 0.5 * dy_1),
                Point(c_x + 0.5 * dx_2, c_y + 0.5 * dy_2),
            ]

            bheight = np.nanmean(
                [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
            )
            if ("hoogte insteek linkerzijde" in params) and (
                "hoogte insteek rechterzijde" in params
            ):
                prof_height = [
                    params["hoogte insteek linkerzijde"],
                    bheight,
                    bheight,
                    params["hoogte insteek rechterzijde"],
                ]
            else:
                prof_height = [def_height, bheight, bheight, def_height]

            p_list = []
            for ix, point in enumerate(list_of_points):
                p_dict = dict(
                    [
                        ("codevolgnummer", ix + 1),
                        ("geometry", Point(point.x, point.y, prof_height[ix])),
                        ("profiel nummer", profiel_nummer),
                        ("type meting", 2),
                        ("typeruwheid", branch["typeruwheid"]),
                        ("ruwheidhoog", branch["ruwheidhoog"]),
                        ("ruwheidlaag", branch["ruwheidlaag"]),
                    ]
                )
                p_list.append(p_dict)

            return p_list

        def trapezium_point_profile(
            branch: LineString,
            profiel_nummer: str,
            params: dict,
            interp_range: float = 0.1,
        ) -> List[Point]:

            bheight = np.nanmean(
                [params["bodemhoogte benedenstrooms"], params["bodemhoogte bovenstrooms"]]
            )

            prof_height = [
                params["hoogte insteek linkerzijde"],
                bheight,
                bheight,
                params["hoogte insteek rechterzijde"],
            ]

            l_depth = params["hoogte insteek linkerzijde"] - bheight
            r_depth = params["hoogte insteek rechterzijde"] - bheight

            # talud is dx/dy
            left_offset = params["taludhelling linkerzijde"] * l_depth
            right_offset = params["taludhelling rechterzijde"] * r_depth

            bwidth = params["bodembreedte"]

            # l_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=bwidth / 2 + left_offset, side="left"
            # )
            # l_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="left")
            # r_offset_line_1 = branch.geometry.parallel_offset(
            #     distance=bwidth / 2 + right_offset, side="right"
            # )
            # r_offset_line_2 = branch.geometry.parallel_offset(distance=bwidth / 2, side="right")

            rotated_line = affinity.rotate(
                branch.geometry, 90, origin=branch.geometry.interpolate(0.5, normalized=True)
            )

            p1 = rotated_line.interpolate(0.5 - interp_range / 2, normalized=True)
            p2 = rotated_line.interpolate(0.5 + interp_range / 2, normalized=True)

            if p2.x > p1.x:
                dx_o = p2.x - p1.x
                dy_o = p2.y - p1.y
            else:
                dx_o = p1.x - p2.x
                dy_o = p1.y - p2.y

            if dx_o != 0:
                s = dy_o / dx_o
                # given a line slope dy/dx (i.e. a/b), pythagoras a**2 + b**2 = c**2, and a desired length c
                # we can determine that b_n = c_n/sqrt(1+(a/b)**2) and a_n = b_n * (a/b)
                dx_1 = bwidth / np.sqrt(1 + s**2)
                dy_1 = dx_1 * s

                dx_2 = (bwidth + left_offset + right_offset) / np.sqrt(1 + s**2)
                dy_2 = dx_2 * s
            else:
                dx_1 = dx_2 = 0
                dy_1 = dy_2 = bwidth

            centroid = branch.geometry.interpolate(0.5, normalized=True)
            # centroid = rotated_line.centroid
            c_x, c_y = centroid.x, centroid.y

            list_of_points = [
                Point(c_x - 0.5 * dx_2, c_y - 0.5 * dy_2),
                Point(c_x - 0.5 * dx_1, c_y - 0.5 * dy_1),
                Point(c_x + 0.5 * dx_1, c_y + 0.5 * dy_1),
                Point(c_x + 0.5 * dx_2, c_y + 0.5 * dy_2),
            ]

            p_list = []
            for ix, point in enumerate(list_of_points):
                try:
                    p_dict = dict(
                        [
                            ("codevolgnummer", ix + 1),
                            ("geometry", Point(point.x, point.y, prof_height[ix])),
                            ("profiel nummer", profiel_nummer),
                            ("type meting", 2),
                            ("typeruwheid", branch["typeruwheid"]),
                            ("ruwheidhoog", branch["ruwheidhoog"]),
                            ("ruwheidlaag", branch["ruwheidlaag"]),
                        ]
                    )
                except:
                    print(params)
                    raise
                p_list.append(p_dict)
            return p_list

        # Load shapefile with profile points & group them by metingprof attribute
        # profile_points = gpd.read_file(profile_points_path)
        branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ] = branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ].astype(
            float
        )

        branches_gdf.loc[
            branches_gdf["bodembreedte"] < min_water_width,
            "bodembreedte",
        ] = np.nan
        branches_gdf.loc[
            branches_gdf["water_width_index"] < min_water_width,
            "water_width_index",
        ] = np.nan

        branches_gdf["globalid"] = [str(uuid.uuid4()) for _ in range(branches_gdf.shape[0])]
        branches_gdf.set_index("globalid", drop=False, inplace=True)
        branches_out_gdf = copy(branches_gdf)

        # initialize empty lists
        pp_list = []
        # Loop over the grouped profile points to create a line from the points
        for jx, ix_1 in enumerate(branches_gdf.index):
            branch = branches_gdf.loc[ix_1, :]
            # Dont use jx with numerical index column

            # # Create unique ident for branch and add to branch
            branch_gid = ix_1

            # # replace duplicate codes by addinng a uniqure number
            # if duplicates.iloc[jx]:
            #     branches_out_gdf.at[ix_1, "code"] = branches_out_gdf.at[ix_1, "code"] + str(jx)

            # turn numerical roughnes types to strings
            type_ruwheid = check_roughness(branch)
            branches_out_gdf.at[ix_1, "typeruwheid"] = type_ruwheid

            # check if linestring
            try:
                assert isinstance(branch.geometry, LineString)
            except AssertionError:
                print(branch.geometry.type)
                raise

            # Delete z-dimensison of linestring if it exists
            if np.array(branch.geometry.coords).shape[1] > 2:
                branches_out_gdf.at[ix_1, "geometry"] = LineString(
                    [xy[0:2] for xy in list(branch.geometry.coords)]
                )

            # Check if width and depth parameters are available, if not, skip
            # Also skip if no width is available
            if (
                (
                    (branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]].empty)
                    or (
                        branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]]
                        .isna()
                        .values.all()
                    )
                )
                or (
                    (
                        (branch[["bodembreedte", "bodemhoogte benedenstrooms"]].empty)
                        or (
                            branch[["bodembreedte", "bodemhoogte benedenstrooms"]]
                            .isna()
                            .values.all()
                        )
                    )
                    and (
                        (branch[["water_width_index", "bodemhoogte benedenstrooms"]].empty)
                        or (
                            branch[["water_width_index", "bodemhoogte benedenstrooms"]]
                            .isna()
                            .values.all()
                        )
                    )
                )
                or (
                    (branch[["bodembreedte", "water_width_index"]].empty)
                    or (branch[["bodembreedte", "water_width_index"]].isna().values.all())
                )
            ):
                # print("no profile for branch {}".format(branch_gid))
                continue

            # Check if all parameters are present for trapezium profile. If not, use rectangular profile
            width_ix = "bodembreedte"
            if np.isnan(branch["bodembreedte"]):
                prof_type = "rectangle"
                width_ix = "water_width_index"
            elif (
                (not check_is_not_na_number(branch["hoogte insteek linkerzijde"]))
                | (not check_is_not_na_number(branch["hoogte insteek rechterzijde"]))
                | (not check_is_not_na_number(branch["taludhelling linkerzijde"]))
                | (not check_is_not_na_number(branch["taludhelling rechterzijde"]))
            ):
                prof_type = "rectangle"
            else:
                prof_type = "trapezium"

            width = branch[width_ix]
            if width < 0:
                width = 0
            elif width > 200:
                width = 0

            # set required parameters for either rectangle or trapezium profile
            if prof_type == "rectangle":
                params = dict(
                    [
                        ("bodembreedte", width),
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
                if check_is_not_na_number(branch["hoogte insteek linkerzijde"]):
                    params["hoogte insteek linkerzijde"] = branch["hoogte insteek linkerzijde"]

                if check_is_not_na_number(branch["hoogte insteek rechterzijde"]):
                    params["hoogte insteek rechterzijde"] = branch["hoogte insteek rechterzijde"]

                p_list = rectangular_point_profile(
                    branch=branch, profiel_nummer=ix_1, params=params
                )
                pp_list += p_list
            elif prof_type == "trapezium":
                params = dict(
                    [
                        ("bodembreedte", width),
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

                p_list = trapezium_point_profile(branch=branch, profiel_nummer=ix_1, params=params)
                pp_list += p_list

        # Create GeoDataFrame from list of dicts
        pp_gdf = gpd.GeoDataFrame(data=pp_list, geometry="geometry", crs=branches_gdf.crs)
        return create_measured_profile_data(profile_points_gdf=pp_gdf, dist_tol=0)

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
        branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ] = branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ].astype(
            float
        )

        branches_gdf.loc[
            branches_gdf["bodembreedte"] < min_water_width,
            "bodembreedte",
        ] = np.nan
        branches_gdf.loc[
            branches_gdf["water_width_index"] < min_water_width,
            "water_width_index",
        ] = np.nan

        branches_gdf["globalid"] = [str(uuid.uuid4()) for _ in range(branches_gdf.shape[0])]
        branches_gdf.set_index("globalid", drop=False, inplace=True)
        branches_out_gdf = copy(branches_gdf)

        # # Check for duplicate codes, if there are, they are replaced in the following loop
        # duplicates = branches_gdf.duplicated(subset=["code"])
        # for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
        # for jx, (ix_1, branch) in enumerate(branches_gdf.iterrows()):
        for jx, ix_1 in enumerate(branches_gdf.index):
            branch = branches_gdf.loc[ix_1, :]
            # Dont use jx with numerical index column

            # # Create unique ident for branch and add to branch
            branch_gid = ix_1

            # # replace duplicate codes by addinng a uniqure number
            # if duplicates.iloc[jx]:
            #     branches_out_gdf.at[ix_1, "code"] = branches_out_gdf.at[ix_1, "code"] + str(jx)

            # turn numerical roughnes types to strings
            type_ruwheid = check_roughness(branch)
            branches_out_gdf.at[ix_1, "typeruwheid"] = type_ruwheid

            # check if linestring
            assert isinstance(branch.geometry, LineString)

            # Delete z-dimensison of linestring if it exists
            if np.array(branch.geometry.coords).shape[1] > 2:
                branches_out_gdf.at[ix_1, "geometry"] = LineString(
                    [xy[0:2] for xy in list(branch.geometry.coords)]
                )

            # Check if width and depth parameters are available, if not, skip
            # Also skip if no width is available
            if (
                (
                    (branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]].empty)
                    or (
                        branch[["bodemhoogte benedenstrooms", "bodemhoogte bovenstrooms"]]
                        .isna()
                        .values.all()
                    )
                )
                or (
                    (
                        (branch[["bodembreedte", "bodemhoogte benedenstrooms"]].empty)
                        or (
                            branch[["bodembreedte", "bodemhoogte benedenstrooms"]]
                            .isna()
                            .values.all()
                        )
                    )
                    and (
                        (branch[["water_width_index", "bodemhoogte benedenstrooms"]].empty)
                        or (
                            branch[["water_width_index", "bodemhoogte benedenstrooms"]]
                            .isna()
                            .values.all()
                        )
                    )
                )
                or (
                    (branch[["bodembreedte", "water_width_index"]].empty)
                    or (branch[["bodembreedte", "water_width_index"]].isna().values.all())
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
            width_ix = "bodembreedte"
            if np.isnan(branch["bodembreedte"]):
                prof_type = "rectangle"
                width_ix = "water_width_index"
            # elif (
            #     (branch["hoogte insteek linkerzijde"] is None)
            #     | (branch["hoogte insteek rechterzijde"] is None)
            #     | (branch["taludhelling linkerzijde"] is None)
            #     | (branch["taludhelling rechterzijde"] is None)
            # ):
            #     prof_type = "rectangle"
            # elif (
            #     np.isnan(branch["hoogte insteek linkerzijde"])
            #     | np.isnan(branch["hoogte insteek rechterzijde"])
            #     | np.isnan(branch["taludhelling linkerzijde"])
            #     | np.isnan(branch["taludhelling rechterzijde"])
            # ):
            #     prof_type = "rectangle"
            # elif (
            #     (branch["hoogte insteek linkerzijde"] == 0)
            #     | (branch["hoogte insteek rechterzijde"] == 0)
            #     | (branch["taludhelling linkerzijde"] == 0)
            #     | (branch["taludhelling rechterzijde"] == 0)
            # ):
            elif (
                (not check_is_not_na_number(branch["hoogte insteek linkerzijde"]))
                | (not check_is_not_na_number(branch["hoogte insteek rechterzijde"]))
                | (not check_is_not_na_number(branch["taludhelling linkerzijde"]))
                | (not check_is_not_na_number(branch["taludhelling rechterzijde"]))
            ):
                prof_type = "rectangle"
            else:
                prof_type = "trapezium"

            width = branch[width_ix]
            if width < 0:
                width = 0
            elif width > 200:
                width = 0

            # set required parameters for either rectangle or trapezium profile
            if prof_type == "rectangle":
                params = dict(
                    [
                        ("bodembreedte", width),
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
                        ("bodembreedte", width),
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

        if len(ho_ngp_list) > 0:
            hydroobject_normgp = gpd.GeoDataFrame(
                ho_ngp_list, geometry="geometry", crs=branches_out_gdf.crs
            )
        else:
            hydroobject_normgp = None

        if len(ngp_list) > 0:
            normgeparamprofielwaarde = gpd.GeoDataFrame(
                ngp_list, geometry="geometry", crs=branches_out_gdf.crs
            )
        else:
            normgeparamprofielwaarde = None

        return (
            branches_out_gdf[["code", "globalid", "geometry", "typeruwheid"]].reset_index(
                drop=True
            ),
            hydroobject_normgp,
            normgeparamprofielwaarde,
        )

    def create_pump_data(pump_gdf: gpd.GeoDataFrame) -> List[gpd.GeoDataFrame]:
        pump_station_list = []
        pump_list = []
        management_list = []

        pump_gdf["maximalecapaciteit"] = check_column_is_numerical(
            gdf=pump_gdf["maximalecapaciteit"]
        )

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

    def create_river_profiles(
        branches_gdf: gpd.GeoDataFrame, riv_prof_df: pd.DataFrame, code_padding: str
    ) -> gpd.GeoDataFrame:
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
            name = code_padding + u_id
            slice = riv_prof_df[riv_prof_df["id"] == u_id]
            meta_data = slice[slice["Data_type"] == "meta"]
            geom_data = slice[slice["Data_type"] == "geom"]
            # print(meta_data.dropna(axis=1))
            # print(geom_data.dropna(axis=1))
            if np.sum(branches_gdf["code"].str.contains(meta_data["branch"].values[0])) == 0:
                continue

            for ix in range(geom_data.shape[0]):
                geom_df.at[g_ix, "name"] = name
                geom_df.at[g_ix, "ix"] = ix
                geom_df.at[g_ix, "levels"] = geom_data["level"].values[ix]
                geom_df.at[g_ix, "flowWidths"] = geom_data["Flow width"].values[ix]
                geom_df.at[g_ix, "totalWidths"] = geom_data["Total width"].values[ix]
                geom_df.at[g_ix, "geometry"] = None
                g_ix += 1

            meta_df.at[name, "numLevels"] = geom_data.shape[0]
            meta_df.at[name, "name"] = name
            meta_df.at[name, "thalweg"] = 0
            meta_df.at[name, "leveecrestLevel"] = meta_data["Crest level summerdike"].values[0]

            meta_df.at[name, "leveebaselevel"] = meta_data[
                "Floodplain baselevel behind summerdike"
            ].values[0]
            meta_df.at[name, "leveeflowarea"] = meta_data["Flow area behind summerdike"].values[0]
            meta_df.at[name, "leveetotalarea"] = meta_data["Total area behind summerdike"].values[
                0
            ]
            meta_df.at[name, "mainwidth"] = meta_data["width main channel"].values[0]
            meta_df.at[name, "fp1width"] = meta_data["width floodplain 1"].values[0]
            meta_df.at[name, "fp2width"] = meta_data["width floodplain 2"].values[0]
            # meta_df.at[u_id, "branchid"] = meta_data["branch"].values[0]
            # account for padded codes
            meta_df.at[name, "branchid"] = branches_gdf.loc[
                branches_gdf["code"].str.contains(meta_data["branch"].values[0]), "code"
            ].values[0]

            meta_df.at[name, "chainage"] = meta_data["chainage"].values[0]
            meta_df.at[name, "geometry"] = None

            total_width = (
                meta_df.at[name, "mainwidth"]
                + meta_df.at[name, "fp1width"]
                + meta_df.at[name, "fp2width"]
            )
            max_width = geom_df.loc[geom_df["name"] == name, "flowWidths"].max()
            if total_width > max_width:
                diff = total_width - max_width

                max_ix = geom_df.loc[geom_df["name"] == name, "flowWidths"].idxmax()
                geom_df.at[max_ix, "flowWidths"] = geom_df.at[max_ix, "flowWidths"] + diff
                if geom_df.at[max_ix, "flowWidths"] > geom_df.at[max_ix, "totalWidths"]:
                    geom_df.at[max_ix, "totalWidths"] = geom_df.at[max_ix, "flowWidths"] + 1

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

        weir_gdf["hoogstedoorstroombreedte"].fillna(
            weir_gdf["laagstedoorstroombreedte"], inplace=True
        )
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

            elif isinstance(
                weir["soortregelbaarheid"], str
            ):  # dtype not a number, so assuming string
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

    def fill_branch_norm_parm_profiles_data(
        defaults, in_branches_gdf: gpd.GeoDataFrame, data_config, insteek_marge=0.25
    ) -> gpd.GeoDataFrame:
        """ """

        def find_branch_width(branch_geom, buffer_list, name, watervlak_geometry):
            for jx, buffer in enumerate(buffer_list):
                buffered_branches = branch_geom.buffer(buffer, cap_style=2)
                polygon = buffered_branches.intersection(watervlak_geometry)

                if jx > 0:
                    overlap_area = polygon.area / buffered_branches.area
                    if overlap_area < 0.9:
                        width = buffer_list[jx - 1] * 2
                        break
                    else:
                        width = buffer * 2

            return dict([("id", name), ("overlap", overlap_area), ("width", width)])

        in_branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ] = in_branches_gdf[
            [
                "bodembreedte",
                "bodemhoogte benedenstrooms",
                "bodemhoogte bovenstrooms",
                "hoogte insteek linkerzijde",
                "hoogte insteek rechterzijde",
                "taludhelling linkerzijde",
                "taludhelling rechterzijde",
            ]
        ].astype(
            float
        )

        if (
            hasattr(data_config, "peil_gebieden_path")
            and ("diepte" in data_config.branch_index_mapping)
            and (data_config.branch_index_mapping["diepte"] is not None)
        ):
            peil_gebieden_gdf = gpd.read_file(data_config.peil_gebieden_path)
            peil_gebieden_gdf, _ = map_columns(
                defaults=defaults,
                gdf=peil_gebieden_gdf,
                index_mapping=data_config.peil_index_mapping,
            )

            out_branches_gdf = copy(in_branches_gdf)

            in_branches_with_peil = gpd.sjoin(
                in_branches_gdf, peil_gebieden_gdf, how="left", predicate="intersects"
            )

            for name, branch in in_branches_with_peil.iterrows():
                gid = branch["globalid"]
                bool_ix = out_branches_gdf["globalid"] == gid

                if in_branches_with_peil[in_branches_with_peil["globalid"] == gid].shape[0] > 1:
                    # check for duplicates
                    pass

                if (not check_is_not_na_number(branch["hoogte insteek linkerzijde"])) and (
                    check_is_not_na_number(branch["hoogte insteek rechterzijde"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "hoogte insteek linkerzijde"
                    ] = out_branches_gdf.at[bool_ix, "hoogte insteek rechterzijde"]
                elif (check_is_not_na_number(branch["hoogte insteek linkerzijde"])) and (
                    not check_is_not_na_number(branch["hoogte insteek rechterzijde"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "hoogte insteek rechterzijde"
                    ] = out_branches_gdf.at[bool_ix, "hoogte insteek linkerzijde"]
                elif (not check_is_not_na_number(branch["hoogte insteek linkerzijde"])) and (
                    not check_is_not_na_number(branch["hoogte insteek rechterzijde"])
                ):
                    if check_is_not_na_number(branch["vast peil"]):
                        peil = branch["vast peil"]
                    else:
                        peil = np.nanmean([branch["boven peil"], branch["onder peil"]])

                    insteek_hoogte = peil + insteek_marge

                    out_branches_gdf.loc[bool_ix, "hoogte insteek linkerzijde"] = insteek_hoogte
                    out_branches_gdf.loc[bool_ix, "hoogte insteek rechterzijde"] = insteek_hoogte

                if (not check_is_not_na_number(branch["bodemhoogte benedenstrooms"])) and (
                    check_is_not_na_number(branch["bodemhoogte bovenstrooms"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "bodemhoogte benedenstrooms"
                    ] = out_branches_gdf.at[bool_ix, "bodemhoogte bovenstrooms"]
                elif (check_is_not_na_number(branch["bodemhoogte benedenstrooms"])) and (
                    not check_is_not_na_number(branch["bodemhoogte bovenstrooms"])
                ):
                    out_branches_gdf.loc[
                        bool_ix, "bodemhoogte bovenstrooms"
                    ] = out_branches_gdf.at[bool_ix, "bodemhoogte benedenstrooms"]
                elif (not check_is_not_na_number(branch["bodemhoogte benedenstrooms"])) and (
                    not check_is_not_na_number(branch["bodemhoogte bovenstrooms"])
                ):
                    diepte = branch["diepte"]
                    insteek_hoogte = out_branches_gdf.loc[bool_ix, "hoogte insteek linkerzijde"]

                    bodem_hoogte = insteek_hoogte - diepte

                    if diepte == -99:
                        diepte = np.nan

                    out_branches_gdf.loc[bool_ix, "bodemhoogte benedenstrooms"] = bodem_hoogte
                    out_branches_gdf.loc[bool_ix, "bodemhoogte bovenstrooms"] = bodem_hoogte

            in_branches_gdf = copy(out_branches_gdf)

        if hasattr(data_config, "watervlak_path"):
            watervlak_gdf = gpd.read_file(data_config.watervlak_path, geometry="geometry")
            watervlak_geometry = watervlak_gdf.dissolve(by=None).geometry.values[0]

            buffer_list = [1.25, 2.5, 5, 10, 15, 20, 25, 50, 100]
            # cpus = 8
            # pool = Pool(processes=cpus)
            results = []
            save = False
            for ix, (name, branch) in enumerate(in_branches_gdf.iterrows()):
                if (branch["bodembreedte"] is None) and (branch["water_width_index"] is None):
                    if not save:
                        save = True
                    branch_geom = branch.geometry

                    kwds = dict(
                        [
                            ("branch_geom", branch_geom),
                            ("buffer_list", buffer_list),
                            ("name", name),
                            ("watervlak_geometry", watervlak_geometry),
                        ]
                    )
                    # results.append(pool.apply_async(find_branch_width, kwds=kwds))
                    results.append(find_branch_width(**kwds))

            for res in results:
                # r.wait()
                # res = r.get()
                out_branches_gdf.at[res["id"], "bodembreedte"] = res["width"]

            if save:
                out_branches_gdf.to_file(data_config.branches_path.replace(".shp", "_ww.shp"))

        if "out_branches_gdf" in locals():
            return out_branches_gdf
        else:
            return in_branches_gdf

    def load_geo_file(file_path: str, layer: str = None):
        """ """
        if type(file_path) == list:
            if len(file_path) <= 2:
                gdf_1 = gpd.read_file(file_path[0])
                gdf_2 = gpd.read_file(file_path[1])

                gdf_2_buffer = gdf_2.copy()
                gdf_2_buffer.geometry = gdf_2_buffer.geometry.buffer(1)
                gdf = gdf_1.sjoin(gdf_2_buffer, how="left", predicate="within")

            else:
                raise ValueError("Can only join two shapefiles together")
        elif file_path.endswith(r".shp"):
            gdf = gpd.read_file(file_path)
        elif file_path.endswith(r".gpkg"):
            if layer is not None:
                gdf = gpd.read_file(file_path, layer=layer)
            else:
                raise ValueError("provide a layer when loading a gpkg")
        else:
            raise ValueError("filetype not implemented")

        gdf.set_geometry(col="geometry", inplace=True)

        # check for multi-element geometries
        if np.sum(gdf.geometry.type.isin(["MultiPoint", "MultiLine", "MultiLineString"])) > 0:
            gdf = gdf.explode(ignore_index=True, index_parts=False)
        return gdf

    def map_columns(
        defaults, gdf: gpd.GeoDataFrame, index_mapping: dict, code_pad: str = None
    ) -> gpd.GeoDataFrame:
        """ """

        def fill_empty_columns(
            defaults, gdf: gpd.GeoDataFrame, index_mapping: dict
        ) -> Tuple[gpd.GeoDataFrame, dict]:
            """ """

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
                    if isinstance(value, list):
                        _gdf[key] = (
                            gdf[value]
                            .replace(to_replace=[None, "n.v.t."], value=np.nan)
                            .astype(float)
                            .apply(np.nanmean, axis=1)
                        )
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

                            print(
                                "{}: entries with missing data: {}".format(
                                    key, n_zero + n_nnn + n_nan
                                )
                            )

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

        # ensure no incorrect geometries exist and set a globalid
        _gdf = fix_geometry(gdf)
        _gdf["globalid"] = [str(uuid.uuid4()) for _ in range(_gdf.shape[0])]

        # fill potential empty columns
        _gdf, index_mapping = fill_empty_columns(
            defaults=defaults,
            gdf=_gdf,
            index_mapping=index_mapping,
        )

        if ("code" in _gdf.columns) and (code_pad is not None):
            _gdf["code"] = code_pad + _gdf["code"].astype(str)
            # check for dupliacte columns in code
            _gdf = check_and_fix_duplicate_code(_gdf)

        # gdf = gpd.GeoDataFrame(_gdf[list(index_mapping.keys())], geometry="geometry", crs=gdf.crs)
        gdf = _gdf[list(index_mapping.keys())]

        # return only columns that are in that dict
        return gdf, index_mapping

    def merge_to_ddm(
        ddm: DHydamoDataModel, feature: str, feature_gdf: gpd.GeoDataFrame
    ) -> DHydamoDataModel:
        if hasattr(ddm, feature):
            new_gdf = gpd.GeoDataFrame(
                data=pd.concat([getattr(ddm, feature), feature_gdf]),
                geometry="geometry",
                crs=profielgroep.crs,
            )
            setattr(ddm, feature, new_gdf)
        else:
            setattr(ddm, feature, feature_gdf)

        return ddm

    ddm = DHydamoDataModel()
    defaults = importlib.import_module("dataset_configs." + defaults)
    data_config = getattr(importlib.import_module("dataset_configs." + config), "RawData")
    code_padding = config[:4] + r"_"  # add prefix of length 4 to all objects with codes

    if hasattr(data_config, "branches_path"):
        if data_config.branches_path is not None:
            ## Branches
            branches_gdf = load_geo_file(data_config.branches_path, layer="waterloop")
            if hasattr(data_config, "branch_selection"):
                column, value = (
                    data_config.branch_selection["column"],
                    data_config.branch_selection["value"],
                )
                branches_gdf = branches_gdf.loc[branches_gdf[column] == value, :]

            branches_gdf, index_mapping = map_columns(
                code_pad=code_padding,
                defaults=defaults.Branches,
                gdf=branches_gdf,
                index_mapping=data_config.branch_index_mapping,
            )

            if hasattr(data_config, "norm_profile_path") and (
                data_config.norm_profile_path is not None
            ):
                np_gdf = load_geo_file(data_config.norm_profile_path)

                if hasattr(data_config, "np_selection"):
                    column, value = (
                        data_config.np_selection["column"],
                        data_config.np_selection["value"],
                    )
                    np_gdf = np_gdf.loc[np_gdf[column] == value, :]

                np_gdf, index_mapping = map_columns(
                    code_pad=code_padding,
                    defaults=defaults.Branches,
                    gdf=np_gdf,
                    index_mapping=data_config.np_index_mapping,
                )
                np_gdf = fill_branch_norm_parm_profiles_data(
                    defaults=defaults.Peil, in_branches_gdf=np_gdf, data_config=data_config
                )

                (
                    ddm.waterloop,
                    ddm.hydroobject_normgp,
                    ddm.normgeparamprofielwaarde,
                ) = create_norm_parm_profiles_v2(branches_gdf=branches_gdf)

                (
                    ddm.profielgroep,
                    ddm.profiellijn,
                    ddm.profielpunt,
                    ddm.ruwheidsprofiel,
                ) = create_mp_from_np(branches_gdf=np_gdf)
            else:
                (
                    ddm.waterloop,
                    _,
                    _,
                ) = create_norm_parm_profiles_v2(branches_gdf=branches_gdf)

    if hasattr(data_config, "bridges_path"):
        if data_config.bridges_path is not None:
            ## Bridges
            bridges_gdf, _ = map_columns(
                code_pad=code_padding,
                defaults=defaults.Bridges,
                gdf=load_geo_file(data_config.bridges_path, layer="brug"),
                index_mapping=data_config.bridge_index_mapping,
            )
            ddm.brug = create_bridge_data(bridge_gdf=bridges_gdf)

    if hasattr(data_config, "culvert_path"):
        if data_config.culvert_path is not None:
            ## Culverts
            culvert_gdf, _ = map_columns(
                code_pad=code_padding,
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
                code_pad=code_padding,
                defaults=defaults.MeasuredProfiles,
                gdf=measure_profile_gdf,
                index_mapping=data_config.measured_profile_index_mapping,
            )
            (
                profielgroep,
                profiellijn,
                profielpunt,
                ruwheidsprofiel,
            ) = create_measured_profile_data(profile_points_gdf=measure_profile_gdf)

            ddm = merge_to_ddm(ddm=ddm, feature="profielgroep", feature_gdf=profielgroep)
            ddm = merge_to_ddm(ddm=ddm, feature="profiellijn", feature_gdf=profiellijn)
            ddm = merge_to_ddm(ddm=ddm, feature="profielpunt", feature_gdf=profielpunt)
            ddm = merge_to_ddm(ddm=ddm, feature="ruwheidsprofiel", feature_gdf=ruwheidsprofiel)

            # if hasattr(ddm, "profielgroep"):
            #     ddm.profielgroep = gpd.GeoDataFrame(
            #         data=pd.concat([ddm.profielgroep, profielgroep]),
            #         geometry="geometry",
            #         crs=profielgroep.crs,
            #     )
            # else:
            #     ddm.profielgroep = profielgroep

            # if hasattr(ddm, "profiellijn"):
            #     ddm.profiellijn = gpd.GeoDataFrame(
            #         data=pd.concat([ddm.profiellijn, profiellijn]),
            #         geometry="geometry",
            #         crs=profiellijn.crs,
            #     )
            # else:
            #     ddm.profiellijn = profiellijn

            # if hasattr(ddm, "profielpunt"):
            #     ddm.profielpunt = gpd.GeoDataFrame(
            #         data=pd.concat([ddm.profielpunt, profielpunt]),
            #         geometry="geometry",
            #         crs=profielpunt.crs,
            #     )
            # else:
            #     ddm.profielpunt = profielpunt

            # if hasattr(ddm, "ruwheidsprofiel"):
            #     ddm.ruwheidsprofiel = gpd.GeoDataFrame(
            #         data=pd.concat([ddm.ruwheidsprofiel, ruwheidsprofiel]),
            #         geometry="geometry",
            #         crs=ruwheidsprofiel.crs,
            #     )
            # else:
            #     ddm.ruwheidsprofiel = ruwheidsprofiel

    if hasattr(data_config, "pump_path"):
        if data_config.pump_path is not None:
            ## Pumps
            pump_gdf = load_geo_file(data_config.pump_path, layer="gemaal")
            pump_gdf, _ = map_columns(
                code_pad=code_padding,
                defaults=defaults.Pumps,
                gdf=pump_gdf,
                index_mapping=data_config.pump_index_mapping,
            )

            ddm.gemaal, ddm.pomp, ddm.sturing = create_pump_data(pump_gdf=pump_gdf)

    if hasattr(data_config, "river_profile_path"):
        if data_config.river_profile_path is not None:
            riv_prof_df = pd.read_csv(data_config.river_profile_path)
            ddm.rivier_profielen, ddm.rivier_profielen_data = create_river_profiles(
                branches_gdf=branches_gdf,
                riv_prof_df=riv_prof_df,
                code_padding=code_padding,
            )

    if hasattr(data_config, "sluice_path"):
        if data_config.sluice_path is not None:
            sluice_gdf = load_geo_file(data_config.sluice_path, layer="sluis")

    if hasattr(data_config, "weir_path"):
        if data_config.weir_path is not None:
            ## Weirs
            weir_gdf = load_geo_file(data_config.weir_path, layer="stuw")
            weir_gdf, _ = map_columns(
                code_pad=code_padding,
                defaults=defaults.Weirs,
                gdf=weir_gdf,
                index_mapping=data_config.weir_index_mapping,
            )

            ddm.stuw, ddm.kunstwerkopening, ddm.regelmiddel = create_weir_data(weir_gdf=weir_gdf)

    return ddm
