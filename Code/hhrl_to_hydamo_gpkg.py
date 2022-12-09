#%%
# Dit script zet HH Rijnland data om naar een geschikt geopackage format voor HyDAMO.
# Stappen die genomen worden
# - Voor verschillende elementen gegevens inladen en ordenen naar HyDAMO formaat
# - Opslaan als nieuwe .shp file
# - Samenvoegen in een gpkg (opslagtechnisch en overzichtelijker)

# Importeer de benodige packages
import geopandas as gpd
import numpy as np
import pandas as pd
import data_functions as daf

# Definieer locatie waar bestanden staan 
folder_data = r"D:\work\P1414_ROI\GIS\HHRijnland\Legger"

# Definieer ruwheidtypes
roughness_list = [
    "Bos en Bijkerk",
    "Chezy",
    "Manning",
    "StricklerKn",
    "StricklerKs",
    "White Colebrook",
]
#### WATERGANGEN ####
filename = r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHR.shp"
hydroobject = gpd.read_file(filename)

hydroobjectsub = hydroobject[['CODE','RUWHEIDSWA', 'TYPERUWHEI','geometry']]
hydroobjectsub = hydroobjectsub.rename(columns={'CODE':'globalid',
                                                'RUWHEIDSWA':'ruwheid',
                                                'TYPERUWHEI':'typeruwheid'})
for n,typeruwheid in enumerate(hydroobjectsub['typeruwheid']):
    hydroobjectsub['typeruwheid'][n] = roughness_list.index(typeruwheid)

hydroobjectsub['code'] = daf.getuniquecode('HHRL_watergang_',len(hydroobjectsub['globalid']))

#### BRUGGEN ####
filename = r"\Brug\brug_3.shp"
bridges = gpd.read_file(folder_data+filename)
bridgessub = bridges[['CODE','DOORSTRO_1','geometry']]
bridgessub = bridgessub.rename(columns={'CODE':'code',
                                        'DOORSTRO_1':'lengte'})
bridgessub['globalid'] = daf.getuniquecode('HHRL_brug_',len(bridgessub['code']))                                       
bridgessub['typeruwheid'] = 4
bridgessub['ruwheid'] = 75.
bridgessub['intreeverlies'] = 0.5
bridgessub['uittreeverlies'] = 0.7
#%%
#### DUIKERS ####
filename_duiker = "\Duiker\duiker.shp"
filename_sifon = "\Sifon\sifon_2.shp"

culverts = gpd.read_file(folder_data + filename_duiker)
sifons = gpd.read_file(folder_data + filename_sifon)

culvertssub = culverts[['CODE','VORMKOKER','LENGTE',
                        'HOOGTEOPEN','BREEDTEOPE','HOOGTEBINN',
                        'HOOGTEBI_1','geometry']]
culvertssub = culvertssub.rename(columns={'HOOGTEBINN':'hoogtebinnenonderkantbene',
                                          'HOOGTEBI_1':'hoogtebinnenonderkantbov',
                                          'CODE':'code',
                                          'HOOGTEOPEN':'hoogteopening',
                                          'BREEDTEOPE':'breedteopening',
                                          'VORMKOKER':'vormkoker',
                                          'LENGTE':'lengte'})
culvertssub['globalid'] = np.nan
for i in range(len(culvertssub['globalid'])):
    culvertssub['globalid'].iloc[i] = 'HHRL_duiker_'+str(culverts['ORACLE_OBJ'].iloc[i])

sifonssub = sifons[['CODE','VORMKOKER','LENGTE',
                        'HOOGTEOPEN','BREEDTEOPE','HOOGTEBINN',
                        'HOOGTEBI_1','geometry']]
sifonssub = sifonssub.rename(columns={'HOOGTEBINN':'hoogtebinnenonderkantbene',
                                          'HOOGTEBI_1':'hoogtebinnenonderkantbov',
                                          'CODE':'code',
                                          'HOOGTEOPEN':'hoogteopening',
                                          'BREEDTEOPE':'breedteopening',
                                          'VORMKOKER':'vormkoker',
                                          'LENGTE':'lengte'})
sifonssub['globalid'] = np.nan
for i in range(len(sifonssub['globalid'])):
    sifonssub['globalid'].iloc[i] = 'HHRL_duiker_'+str(sifons['ORACLE_OBJ'].iloc[i])

culverts_total = culvertssub.append(sifonssub)
culverts_total['typeruwheid'] = 4
culverts_total['ruwheid'] = 75.
culverts_total['intreeverlies'] = 0.6
culverts_total['uittreeverlies'] = 0.8

# Check if vormkoker is an int or a str
if type(culverts_total.vormkoker.iloc[0]) == str:
    culverts_total.vormkoker = daf.vormkoker_str2int(culverts_total.vormkoker)

# %%
