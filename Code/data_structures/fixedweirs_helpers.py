import importlib

import geopandas as gpd

from data_structures.dhydro_data_model import ROIDataModel
from data_structures.hydamo_helpers import load_geo_file, map_columns, merge_to_dm


def merge_to_fddm(ddm: ROIDataModel, feature: str, feature_gdf: gpd.GeoDataFrame):
    return merge_to_dm(dm=ddm, feature=feature, feature_gdf=feature_gdf)


def create_fixed_weir_data(
    config: str, defaults: str, ddm: ROIDataModel, min_length: float = None
) -> ROIDataModel:

    defaults = importlib.import_module("dataset_configs." + defaults)
    fw_data_config = getattr(importlib.import_module("dataset_configs." + config), "FixedWeirs")

    code_padding = config[:4] + r"_"  # add prefix of length 4 to all objects with codes
    if hasattr(fw_data_config, "flood_defences_path") and (
        fw_data_config.flood_defences_path is not None
    ):
        fd_gdf = load_geo_file(fw_data_config.flood_defences_path, layer="flood_defence")
        fd_gdf, _ = map_columns(
            code_pad=code_padding,
            defaults=defaults.FixedWeirs,
            gdf=fd_gdf,
            index_mapping=fw_data_config.fixed_weir_index_mapping,
        )
        if min_length is not None:
            # fd_gdf = fd_gdf.loc[fd_gdf.geometry.length > min_length, :]
            fd_gdf = fd_gdf.drop(fd_gdf[fd_gdf.geometry.length < min_length].index)

        ddm.keringen = fd_gdf

    return ddm
