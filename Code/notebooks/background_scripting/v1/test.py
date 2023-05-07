# import xarray as xr
# import geopandas as gpd
# import pandas as pd


# ds = xr.open_dataset(r"C:\Werk\Projecten\P1414_ROI\Cursus\StandAloneServiceZipfile\Model_database\V20_WBD_v1\dflowfm\network.nc")
# network_x = [ds.network_node_x.values[i] for i in range(len(ds.network_node_x.values))]
# network_y = [ds.network_node_y.values[i] for i in range(len(ds.network_node_y.values))]


# df = pd.DataFrame({'x': network_x , 'y': network_y})
# geometry = gpd.points_from_xy(df.x, df.y, crs="EPSG:28992")

# gdf = gpd.GeoDataFrame(data = {'index': range(len(network_x))}, geometry = geometry)
# print(gdf)

# gdf.to_file(r"C:\Temp\network_x_y_ROI.shp")

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_inundation_volume(nc_path):
    ncfile = nc.Dataset(nc_path)

    wd = ncfile.variables['Mesh2d_waterdepth'][:]

    wd_over_time = np.sum(wd, axis = 1)[:36]

    # fig, ax = plt.subplots(figsize = (12,12))
    # ax.plot(
    #     [x/3 for x in range(len(wd_over_time))],
    #     wd_over_time
    #     )
    # savefigpath = os.path.dirname(os.path.dirname(os.path.dirname(nc_path)))
    # plt.savefig(os.path.join(savefigpath, 'inundation_total.png'))
    return np.sum(wd_over_time)

runs = r"C:\Werk\Projecten\P1414_ROI\Cursus\StandAloneServiceZipfile\Model_runs"
for folder in os.listdir(runs):
    map_file = f"{folder}/dflowfm/output/DFM_map.nc"
    path = os.path.join(runs, map_file)
    if os.path.exists(path):
        vol = plot_inundation_volume(path)
        # print(vol)
        if vol > 24900:
            print(path)
            print(vol)