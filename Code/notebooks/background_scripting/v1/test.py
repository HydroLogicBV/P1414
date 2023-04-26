import re
from scipy.spatial import KDTree
import numpy as np
import os
import shutil
# from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import json 
from shapely.geometry import Point, LineString
import geopandas as gpd
from ipyleaflet import Map, Marker, projections, basemaps
import ipyleaflet as ipl
import pyproj
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points
import xarray as xr
import matplotlib.pyplot as plt
import netCDF4 as nc
# keringen = gpd.read_file(r"C:\Werk\Projecten\P1414_ROI\Github_P1414\P1414\Code\notebooks\data\Combined_test_v14_WBD.gpkg",
#                     layer = 'keringen',
#                     crs=28992)

# strings_to_match = ['hhsk', 'hhd', 'hdsr', 'hhr', 'wagv', 'ark']
# keringen = keringen[keringen['code'].str.startswith(tuple(strings_to_match))]

# keringen.to_file(r'C:\Werk\Projecten\P1414_ROI\Github_P1414\P1414\Code\notebooks\data\keringen.shp')

# node = "121494.000000_443715.000000"

# waterLevelDownstreamLocationX      = 109664.77608820342
# waterLevelDownstreamLocationY      = 437140.22771585674

# ds = nc.Dataset(r"D:\work\Project\P1414\Models_SAS\Model_runs\V14_krimpenerwaard_v7_run_2023-04-26T12-03-01\dflowfm\DFM_OUTPUT_test\test_map.nc")

# # ids = ds['network_node_id'][:]

# # for i, id in enumerate(ids):
# #     id = ''.join([x.decode('UTF-8') for x in id])
# #     if id.strip(' ') == node:
# #         index = i

# # # fig, ax = plt.subplots(figsize = (10,10))

# print(list(ds.variables.keys()))
location = r"D:\work\Project\P1414\Models_SAS\Model_runs\demo\dflowfm\network.nc"
ds = xr.open_dataset(location)
ds = nc.Dataset(location)

print(list(ds.variables.keys()))
print(ds['mesh1d'])