import uuid

import geopandas as gpd
from shapely.geometry import LineString, MultiPoint, Point


def geometry_check(geometry: gpd.GeoSeries) -> bool:
    # Allow points and linestrings
    if isinstance(geometry, Point) or isinstance(geometry, LineString):
        return True
    # # Verify that multipoints are infact normal points
    # elif isinstance(geometry, MultiPoint):
    #     if len(geometry.geoms) == 1:
    #         return True
    #     else:
    #         return False
    else:
        return False


def globalid_check(globalid: str) -> bool:
    try:
        uuid_obj = uuid.UUID(globalid)
        return True
    except ValueError:
        return False


def none_geometry_check(geometry: gpd.GeoSeries) -> bool:
    # Allow only None types
    if geometry is None:
        return True
    else:
        return False
