import importlib

import geopandas as gpd

from data_structures.hydamo_helpers import load_geo_file, map_columns, merge_to_dm, select_features
from data_structures.roi_data_model import ROIDataModel as DataModel


def merge_to_fddm(ddm: DataModel, feature: str, feature_gdf: gpd.GeoDataFrame):
    return merge_to_dm(dm=ddm, feature=feature, feature_gdf=feature_gdf)


def create_dambreak_data(config: str, defaults: str, dm: DataModel):
    defaults = importlib.import_module("dataset_configs." + defaults)
    db_data_config = getattr(importlib.import_module("dataset_configs." + config), "Dambreak")

    code_padding = config[:4] + r"_"  # add prefix of length 4 to all objects with codes
    if hasattr(db_data_config, "dambreak_path") and (db_data_config.dambreak_path is not None):
        db_gdf = load_geo_file(db_data_config.dambreak_path, layer="flood_defence")
        db_gdf, _ = map_columns(
            code_pad=code_padding,
            defaults=defaults.Dambreak,
            gdf=db_gdf,
            index_mapping=db_data_config.dambreak_index_mapping,
        )

        dm.doorbraak = db_gdf

    return dm


def create_fixed_weir_data(
    config: str,
    defaults: str,
    dm: DataModel,
    min_length: float = None,
) -> DataModel:

    defaults = importlib.import_module("dataset_configs." + defaults)
    fw_data_config = getattr(importlib.import_module("dataset_configs." + config), "FixedWeirs")

    code_padding = config[:4] + r"_"  # add prefix of length 4 to all objects with codes
    if hasattr(fw_data_config, "flood_defences_path") and (
        fw_data_config.flood_defences_path is not None
    ):
        fd_gdf = load_geo_file(fw_data_config.flood_defences_path, layer="flood_defence")
        
        if hasattr(fw_data_config, "fixed_weir_selection"):
            fd_gdf = select_features(fw_data_config.fixed_weir_selection, fd_gdf)
        fd_gdf, _ = map_columns(
            code_pad=code_padding,
            defaults=defaults.FixedWeirs,
            gdf=fd_gdf,
            index_mapping=fw_data_config.fixed_weir_index_mapping,
        )
        if min_length is not None:
            # fd_gdf = fd_gdf.loc[fd_gdf.geometry.length > min_length, :]
            fd_gdf = fd_gdf.drop(fd_gdf[fd_gdf.geometry.length < min_length].index)

        if hasattr(fw_data_config, "dm_attribute"):
            dm_attribute = getattr(fw_data_config, "dm_attribute")
        else:
            dm_attribute = "keringen"

        setattr(dm, dm_attribute, fd_gdf)

    return dm


def create_underpass_data():
    pass
