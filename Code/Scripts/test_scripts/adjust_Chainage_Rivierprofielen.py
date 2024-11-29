#%% This script adjust the chainages of the cross-sections of the Rhine branches Nederrijn and Lek after being split in multiple parts.
# 

import geopandas as gpd
import pandas as pd

# The old rivieren shape:
gdf_old = gpd.read_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken\RTK_Branches.shp")
gdf_new = gpd.read_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken\RTK_BranchesV2.shp")

gdf_old.set_crs('EPSG:28992')
gdf_new.set_crs('EPSG:28992')

csv_ZW_cross = pd.read_csv(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken\ZW_cross_v2.csv")
csv_ruwheid = pd.read_csv(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken\roughness.csv")

# Update the length in the GeomLength Column. 
# This will be used to compute the new chainage
gdf_new['GeomLengthV2'] = round(gdf_new['geometry'].length,3)
print(gdf_new)

csv_ZW_cross_new = csv_ZW_cross.copy()

csv_ZW_cross_new["id"] = csv_ZW_cross_new["id"].astype(str)
unique_ids = csv_ZW_cross_new["id"].unique()

#geom_list = []
#meta_list = []

#%%
# Compute the boundaries for each section of the cut rivier branches:
#idx1 = gdf_new.loc['Lek_1'].index
index = gdf_new[gdf_new["Name"] == "Lek_2"].index[0]

#%%
nederrijn_boundaries = [gdf_new[gdf_new['Name']=='Lek_1', 'GeomLengthV2']]

#%%
for u_id in unique_ids:
    slice = csv_ZW_cross_new[csv_ZW_cross_new["id"] == u_id, :]
    meta_data = slice.loc[slice["Data_type"] == "meta", :]
    geom_data = slice.loc[slice["Data_type"] == "geom", :]




#%%
import geopandas as gpd
agv = gpd.read_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\WAGV\AGV_Onderdoorgangen_Extra\AGV_Onderdoorgangen_Extra.shp")


# %%
