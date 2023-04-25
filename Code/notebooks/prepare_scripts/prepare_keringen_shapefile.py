import geopandas as gpd
import os
import pandas as pd



p_folder = r"D:\Work\Project\P1414\GIS"

arknzk = os.path.join(p_folder, r"Keringen_met_hoogte\ARK.shp")
hdsr =  os.path.join(p_folder, r"Keringen_met_hoogte\hdsr.shp")
hhd = os.path.join(p_folder, r"Keringen_met_hoogte\hhd_zeewering.shp")
hhr = os.path.join(p_folder, r"Keringen_met_hoogte\hhr_primaire_kering.shp")
hhsk =  os.path.join(p_folder, r"Keringen_met_hoogte\hhsk_primaire_kering.shp")

df_list = []
for a in [arknzk, hdsr, hhd, hhr, hhsk]:
    keringen = gpd.read_file(a)
    df_list.append(keringen)

merged = pd.concat(df_list)
merged = merged[['geometry']]
merged.reset_index(inplace = True)

merged.to_file(r"C:\Werk\Projecten\P1414_ROI\Github_P1414\P1414\Code\notebooks\data\keringen_for_interactive_map.shp")
print("done")

