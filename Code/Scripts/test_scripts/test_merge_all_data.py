#%%
# # Testfile for overlapping lines
import geopandas as gpd
import shapely
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import shapelyTools as sT
from importlib import reload
reload(sT)

# Load the watergangen per waterschap
HDSR_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HDSR_v4_test.shp")
AGV_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\AGV_v4_test.shp")
HHDL_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHD_v4_test.shp")
HHSK_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHSK_v4_test.shp")
HHRL_water = gpd.read_file(r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHR_v4_test.shp")

# Create a list of geometries for both datasets
AGV_geom = [x[1].geometry for x in AGV_water.iterrows()]
HDSR_geom = [x[1].geometry for x in HDSR_water.iterrows()]
HHDL_geom = [x[1].geometry for x in HHDL_water.iterrows()]
HHSK_geom = [x[1].geometry for x in HHSK_water.iterrows()]
HHRL_geom = [x[1].geometry for x in HHRL_water.iterrows()]

agv_mls = shapely.geometry.MultiLineString(AGV_geom)
hdsr_mls = shapely.geometry.MultiLineString(HDSR_geom)
hhdl_mls = shapely.geometry.MultiLineString(HHDL_geom)
hhsk_mls = shapely.geometry.MultiLineString(HHSK_geom)
hhrl_mls = shapely.geometry.MultiLineString(HHRL_geom)

# %%
# Connect all the different datasets at the 10 connection points
# Nr 1: HHDL - HHSK in R'dam Noord, Gordelbrug
hhdl_to_hhsk = sT.snap_endpoints(hhdl_mls, hhsk_mls, max_dist = 5)
print('HHDL en HHSK verbonden in R\'dam Noord')
#%%
# Nr 2: HHDL - HHRL in Noordelijke Sluisbrug Leidschendam
hhdl_to_hhrl_hhsk = sT.snap_endpoints(hhdl_to_hhsk, hhrl_mls, max_dist = 5)
print('HHDL en HHRL verbonden bij Noordelijke Sluisbrug, Leidschendam')
#%%
# Nr 3, 4 & 5: AGV - HDSR op de Vecht, bij Breukelen (ARK) en bij Kockengen
HDSR_buffer = hdsr_mls.buffer(20.0, cap_style=2)
AGV_diff = agv_mls.difference(HDSR_buffer)
agv_to_hdsr = sT.snap_endpoints(AGV_diff, hdsr_mls, max_dist = 80)
print('AGV en HDSR verbonden op de Vecht, bij Breukelen en bij Kockengen')
#%%
# Nr 6: RL- HDSR bij Zwammerdam
RL_buffer = hhrl_mls.buffer(10.0, cap_style=2)
HDSR_diff = hdsr_mls.difference(RL_buffer)
hdsr_to_rl = sT.snap_endpoints(HDSR_diff,hhrl_mls,max_dist = 5)
print('HHRL en HDSR verbonden bij Zwammerdam')

#%%
# Nr 7, 8 & 9: RL - AGV bij de Langeraarsche Plas, Haarlemmervaart en Nieuwe Meer
RL_buffer_15 = hhrl_mls.buffer(15.0, cap_style=2)
AGV_diff_rl = agv_to_hdsr.difference(RL_buffer_15)
agv_to_hdsr_hhrl = sT.snap_endpoints(AGV_diff_rl, hhrl_mls, max_dist=25)
print('HHRL en AGV verbonden bij Langeraarsche Plassen, Haarlemmervaart en Nieuwe Meer')
# %%
all_together = list(hhdl_to_hhrl_hhsk) + list(agv_to_hdsr_hhrl) + \
               list(hdsr_to_rl) + list(hhsk_mls) + list(hhrl_mls)

# %%

all_data_buffered = pd.concat([AGV_water, HDSR_water, HHDL_water, HHRL_water, HHSK_water])
all_data_buffered.geometry = all_data_buffered.geometry.buffer(1, cap_style=2)
gdf = gpd.GeoDataFrame(all_together, columns=["geometry"], geometry="geometry", crs=28992)
#%%
# add data from buffered branches if lines in gdf fall within. But keep geometry of branches in gdf
intersected_gdf = gdf.sjoin(all_data_buffered, how="left", predicate="within")
list_of_null = []
for i in intersected_gdf.iterrows():
    if np.isnan(i[1].index_right): # It is a float type if the value is NaN
        subgdf = gdf.iloc[[i[0]]]
        
        # Intersected_geo might consists of multiple polygons, especially when connecting at a node
        # Therefore, another selection procedure is needed to get the right branch
        overlay = gpd.overlay(subgdf, all_data_buffered, how='intersection', keep_geom_type=False)
        max_len = 0
        max_idx = 0
        for j in overlay.iterrows():
            if j[1].geometry.length > max_len:
                max_idx = j[0] 
                max_len = j[1].geometry.length
        
        selected_branch_gdf = overlay.iloc[[max_idx]]

        # With then newly selected branches: perform the intersect based sjoin
        intersected_geo = subgdf.sjoin(selected_branch_gdf, how='left',predicate='intersects')
        intersected_gdf = pd.concat([intersected_gdf, intersected_geo])
    
# Remove the NULL entries that caused a problem in the first place
intersected_gdf = intersected_gdf.dropna(axis=0, subset=['index_right'])     
intersected_gdf.to_file(r'D:\work\P1414_ROI\GIS\test_all_intersected.shp')

print('All data added to the right snapped branches.')
# %%
