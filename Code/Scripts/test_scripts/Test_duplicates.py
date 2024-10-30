import importlib
import uuid
import warnings
from copy import copy
from typing import Any, List, Tuple
import os
from time import time

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.ops import linemerge, nearest_points

from scipy.spatial import KDTree
from tqdm import tqdm
from rtree import index

warnings.filterwarnings(action="ignore", message="Mean of empty slice")

in_branches = gpd.read_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HYDAMO\testfiles\branches_gdf_nearestV2.shp", crs='EPSG:28992')

def normalize_line_direction(geom):
    if geom.geom_type == 'LineString':
        coords = list(geom.coords)
        # Sort to ensure consistent direction: if reversed is smaller, use reversed
        if coords[0] > coords[-1]:  
            coords = coords[::-1]
        return LineString(coords)
    return geom  # Handle cases where geometry might not be LineString

def snap_nodes(in_branches: gpd.GeoDataFrame, geometry_accuracy: float):
    in_branches["t_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]
    in_branches.set_index("t_id", inplace=True)
    out_branches = copy(in_branches)
    for ix, branch in in_branches.iterrows():

        point_list = []
        for x, y, *_ in branch.geometry.coords[:]:
            point_list.append(
                (np.around(x, int(geometry_accuracy)), np.around(y, int(geometry_accuracy)))
            )

        geometry = LineString(point_list)
        if geometry.length > 0:
            out_branches.loc[ix, "geometry"] = LineString(point_list)
        else:
            out_branches = out_branches.drop(index=ix)

    # remove geometries that are double. Use WKT format to speed up
    starttime2 = time()

    # Apply the normalization and create a new column for WKT
    out_branches['normalized_geometry_wkt'] = out_branches['geometry'].apply(lambda geom: normalize_line_direction(geom).wkt)
    #out_branches['geometry_wkt'] = out_branches['geometry'].apply(lambda geom: geom.wkt)
    duplicate_rows = out_branches[out_branches.duplicated(subset='normalized_geometry_wkt', keep=False)]
    out_branches_no_dups = out_branches.drop_duplicates(subset='normalized_geometry_wkt', keep='first').drop(columns='normalized_geometry_wkt')
    endtime2 = time()
    if len(duplicate_rows) > 0:
        print(f"There are {len(duplicate_rows)} duplicate rows, in {round(endtime2-starttime2,1)} sec:")
        print(duplicate_rows)
        print('Only kept 1 entry for each duplicate geometry')
        duplicate_rows = []   # reset for the next objects

    print('Branches that remained: ')
    print(out_branches_no_dups[['geometry', 'TARGET_FID']])
    return out_branches_no_dups

relevant_branch = 'OWA-31044'
relevant_gdf = in_branches[in_branches['CODE']==relevant_branch]
print('Relevant branches: ')
print(relevant_gdf)
branches_gdf_snapped = snap_nodes(in_branches=relevant_gdf, geometry_accuracy=2)

#branches_gdf_snapped['Volgnummer'] = [0,1,2]
branches_gdf_snapped.to_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HYDAMO\testfiles\relevante_branches.shp")

start_points = []
end_points = []

# Loop through each line in the GeoDataFrame
for line in branches_gdf_snapped.geometry:
    if line.geom_type == 'LineString':  # Make sure it's a single LineString
        start_points.append(Point(line.coords[0]))    # Start point
        end_points.append(Point(line.coords[-1]))     # End point
    elif line.geom_type == 'MultiLineString':  # Handle MultiLineString if necessary
        for part in line:
            start_points.append(Point(part.coords[0]))
            end_points.append(Point(part.coords[-1]))

# Create GeoDataFrames for start and end points
start_gdf = gpd.GeoDataFrame(geometry=start_points, crs=branches_gdf_snapped.crs)
end_gdf = gpd.GeoDataFrame(geometry=end_points, crs=branches_gdf_snapped.crs)

# Save as shapefiles
start_gdf.to_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HYDAMO\testfiles\relevante_branches_startpoints.shp")
end_gdf.to_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HYDAMO\testfiles\relevante_branches_endpoints.shp")


gm_list = branches_gdf_snapped.geometry.to_list()
for ix in range(len(gm_list)):
    print(f'Remaining location {ix+1}: {gm_list[ix].coords[:]}')

#difference_1and2 = [(abs(gm_list[0].coords[:][0]-gm_list[1].coords[:][0]),abs(gm_list[0].coords[:][1]-gm_list[1].coords[:][1]))]
#print(difference_1and2)