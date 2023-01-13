from copy import copy

import geopandas as gpd
import numpy as np


def _clip_structures_by_branches(self, buffer: float = 1, min_overlap: float = 0.5):
    if (not hasattr(self, "ddm")) | (not hasattr(self, "features")):
        raise AttributeError("Modeldatabase not loaded")

    buffered_branches = gpd.GeoDataFrame(
        self.ddm.waterloop.dissolve(by=None).buffer(buffer),
        columns=["geometry"],
        crs=self.ddm.waterloop.crs,
        geometry="geometry",
    )

    for feature in self.features:
        if (feature == "waterloop") or (feature == "profiellijn") or (feature == "profielpunt"):
            continue

        ddm_feature = getattr(self.ddm, feature)
        ddm_feature = gpd.GeoDataFrame(
            ddm_feature[ddm_feature["geometry"].notna()],
            crs=ddm_feature.crs,
            geometry="geometry",
        )
        if ddm_feature.shape[0] > 0:
            print(feature)
            print(ddm_feature.shape[0])
            clipped_gdf = gpd.overlay(
                ddm_feature, buffered_branches, how="intersection", keep_geom_type=True
            )
            print(clipped_gdf.shape[0])

            mls_struct_bool = clipped_gdf.geometry.type == "MultiLineString"
            if np.sum(mls_struct_bool) > 0:
                _clipped_gdf = copy(clipped_gdf)

                for ix, struct in clipped_gdf.iterrows():
                    new_length = np.sum(clipped_gdf.loc[struct.name].geometry.length)
                    old_length = np.sum(ddm_feature.loc[struct.name].geometry.length)

                    if mls_struct_bool[ix]:
                        if (new_length / old_length) >= min_overlap:
                            _clipped_gdf.loc[struct.name] = ddm_feature.loc[struct.name]
                        else:
                            _clipped_gdf.drop(index=struct.name, inplace=True)

                    elif not mls_struct_bool[ix]:
                        if (new_length / old_length) < min_overlap:
                            _clipped_gdf.drop(index=struct.name, inplace=True)
                clipped_gdf = _clipped_gdf

            setattr(self.ddm, feature, clipped_gdf)
    return self.ddm
