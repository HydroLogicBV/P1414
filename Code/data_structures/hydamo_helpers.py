import importlib
import uuid
from copy import copy
from typing import List, Tuple

import geopandas as gpd
import numpy as np

from data_structures.dhydamo_data_model import DHydamoDataModel

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


def create_culvert_data(culvert_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """ """
    culvert_gdf["doorstroomopening"] = None
    culvert_gdf["hoogteopening"] = culvert_gdf["hoogteopening"].replace(0, np.nan)
    culvert_gdf["breedteopening"] = culvert_gdf["breedteopening"].replace(0, np.nan)
    _culvert_gdf = copy(culvert_gdf)
    for ix, culvert in culvert_gdf.iterrows():
        if np.isnan(culvert["breedteopening"]) or np.isnan(culvert["hoogteopening"]):
            shape = "circle"
        elif int(culvert["vormkoker"]) == 3:
            shape = "rectangle"
        else:
            shape = "circle"

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

    _culvert_gdf = _culvert_gdf.drop(columns="gesloten")
    return _culvert_gdf


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

    # for ix_1, branch in tqdm(branches_gdf.iterrows(), total=branches_gdf.shape[0]):
    for ix_1, branch in branches_gdf.iterrows():
        # Create unique ident for branch and add to branch
        branch_gid = str(uuid.uuid4())
        branches_gdf.loc[ix_1, "globalid"] = branch_gid
        branches_gdf.loc[ix_1, "code"] = branch_gid

        # change width of zero to nan
        # branches_gdf[index_mapping["bodembreedte"]] = branches_gdf[
        #     index_mapping["bodembreedte"]
        # ].replace(0, np.nan)
        # branches_gdf[index_mapping["water_width_index"]] = branches_gdf[
        #     index_mapping["water_width_index"]
        # ].replace(0, np.nan)

        ## also under 0.03

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

        # turn numerical roughnes types to strings
        type_ruwheid = branch[index_mapping["typeruwheid"]]
        if isinstance(type_ruwheid, int) or isinstance(type_ruwheid, float):
            type_ruwheid = ROUGHNESS_MAPPING[int(type_ruwheid) - 1]

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
        branches_gdf[["code", "globalid", "geometry", "typeruwheid"]],
        hydroobject_normgp,
        normgeparamprofielwaarde,
    )


def create_pump_data(pump_gdf: gpd.GeoDataFrame) -> List[gpd.GeoDataFrame]:

    pump_station_list = []
    pump_list = []
    management_list = []

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


def create_weir_data(weir_gdf: gpd.GeoDataFrame) -> List[gpd.GeoDataFrame]:

    weir_list = []
    opening_list = []
    management_device_list = []

    # if not "hoogstedoorstroombreedte" in weir_gdf.columns:
    #     weir_gdf["hoogstedoorstroombreedte"] = weir_gdf["laagstedoorstroombreedte"]

    # if not "hoogstedoorstroomhoogte" in weir_gdf.columns:
    #     weir_gdf["hoogstedoorstroomhoogte"] = weir_gdf["laagstedoorstroomhoogte"] + 10

    weir_gdf["hoogstedoorstroombreedte"].fillna(weir_gdf["laagstedoorstroombreedte"], inplace=True)
    weir_gdf["hoogstedoorstroomhoogte"].fillna(
        weir_gdf["laagstedoorstroomhoogte"] + 10, inplace=True
    )
    for ix, weir in weir_gdf.iterrows():
        _weir = dict(
            [
                ("afvoercoefficient", weir["afvoercoefficient_stuw"]),
                ("code", weir["code"]),
                ("geometry", weir["geometry"]),
                ("globalid", weir["globalid"]),
                ("soortstuw", weir["soortstuw"]),
            ]
        )
        weir_list.append(_weir)

        opening_gid = str(uuid.uuid4())
        opening = dict(
            [
                ("afvoercoefficient", weir["afvoercoefficient_opening"]),
                ("geometry", weir["geometry"]),
                ("globalid", opening_gid),
                ("hoogstedoorstroombreedte", weir["hoogstedoorstroombreedte"]),
                ("hoogstedoorstroomhoogte", weir["hoogstedoorstroomhoogte"]),
                ("laagstedoorstroombreedte", weir["laagstedoorstroombreedte"]),
                ("laagstedoorstroomhoogte", weir["laagstedoorstroomhoogte"]),
                ("stuwid", weir["globalid"]),
                ("vormopening", weir["vormopening"]),
            ]
        )
        opening_list.append(opening)

        management_device = dict(
            [
                ("code", weir["code"]),
                ("geometry", weir["geometry"]),
                ("globalid", weir["globalid"]),
                ("kunstwerkopeningid", opening_gid),
                ("overlaatonderlaat", weir["overlaatonderlaat"]),
                ("soortregelbaarheid", weir["soortregelbaarheid"]),
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
    default_values, gdf: gpd.GeoDataFrame, index_mapping: dict
) -> Tuple[gpd.GeoDataFrame, dict]:

    _index_mapping = copy(index_mapping)

    # loop over column names in index mapping, with DHYDAMO name in key and column name in shapefile in Value
    for key, value in index_mapping.items():
        # if data is not in shapefile, fill with default values
        if value is None:
            gdf[key] = getattr(default_values, key)
            _index_mapping[key] = key

        else:
            # if it is in shapefile, but contains missing values, fill them as well
            if (gdf[value].dtype == "int64") or (gdf[value].dtype == "float64"):
                if hasattr(default_values, key):
                    gdf[value] = gdf[value].replace(
                        to_replace=[0, -999, np.nan],
                        value=getattr(default_values, key),
                    )
                # print number of 0, -998 and nan if no default value is present
                else:
                    n_zero = gdf[value][gdf[value] == 0].count()
                    n_nnn = gdf[value][gdf[value] == -999].count()
                    n_nan = gdf[value].isna().sum()
                    # print(
                    #     "{}: number of zeros: {}, number of -999: {}, number of nan: {}".format(
                    #         key, n_zero, n_nnn, n_nan
                    #     )
                    # )
                    print("{}: entries with missing data: {}".format(key, n_zero + n_nnn + n_nan))
    return gdf, _index_mapping


def fix_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """drops values without geometry and sets corrects crs"""
    gdf.dropna(axis=0, inplace=True, subset=["geometry"])
    gdf = gdf.to_crs("epsg:28992")
    return gdf


def map_columns(defaults, gdf: gpd.GeoDataFrame, index_mapping: dict) -> gpd.GeoDataFrame:
    """ """
    # ensure no incorrect geometries exist and set a globalid
    gdf = fix_geometry(gdf)
    gdf["globalid"] = [str(uuid.uuid4()) for _ in range(gdf.shape[0])]

    # fill potential empty columns
    gdf, index_mapping = fill_empty_columns(
        default_values=defaults,
        gdf=gdf,
        index_mapping=index_mapping,
    )

    # rename columns to keys of dictionary inde_mapping
    inv_map = {v: k for k, v in index_mapping.items()}
    gdf.rename(columns=inv_map, inplace=True)

    # return only columns that are in that dict
    return gdf[list(index_mapping.keys())]


def convert_to_dhydamo_data(defaults: str, config: str) -> DHydamoDataModel:
    """ """
    ddm = DHydamoDataModel()
    defaults = importlib.import_module("dataset_configs." + defaults)
    data_config = getattr(importlib.import_module("dataset_configs." + config), "RawData")
    features = []

    if hasattr(data_config, "branches_path"):
        ## Branches
        branches_gdf = gpd.read_file(data_config.branches_path)

        branches_gdf, index_mapping = fill_empty_columns(
            default_values=defaults.Branches,
            gdf=branches_gdf,
            index_mapping=data_config.branch_index_mapping,
        )

        (
            ddm.waterloop,
            ddm.hydroobject_normgp,
            ddm.normgeparamprofielwaarde,
        ) = create_norm_parm_profiles(branches_gdf=branches_gdf, index_mapping=index_mapping)

    if hasattr(data_config, "bridges_path"):
        ## Bridges
        ddm.brug = map_columns(
            defaults=defaults.Bridges,
            gdf=gpd.read_file(data_config.bridges_path),
            index_mapping=data_config.bridge_index_mapping,
        )

    if hasattr(data_config, "culvert_path"):
        ## Culverts
        culvert_gdf = map_columns(
            defaults=defaults.Culverts,
            gdf=gpd.read_file(data_config.culvert_path),
            index_mapping=data_config.culvert_index_mapping,
        )
        ddm.duiker = create_culvert_data(culvert_gdf=culvert_gdf)

    if hasattr(data_config, "pump_path"):
        ## Pumps
        pump_gdf = gpd.read_file(data_config.pump_path)
        pump_gdf = map_columns(
            defaults=defaults.Pumps,
            gdf=pump_gdf,
            index_mapping=data_config.pump_index_mapping,
        )

        ddm.gemaal, ddm.pomp, ddm.sturing = create_pump_data(pump_gdf=pump_gdf)

    if hasattr(data_config, "weir_path"):
        ## Weirs
        weir_gdf = gpd.read_file(data_config.weir_path)
        weir_gdf = map_columns(
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
