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
#### STUWEN ####
filename = "\Stuw\stuw.shp"
weirs = gpd.read_file(folder_data+filename)

weirssub = weirs[['CODE','gml_id','geometry']]
weirssub = weirssub.rename(columns={'gml_id':'globalid',
                                    'CODE':'code'})
weirssub['afvoercoefficient'] = 1.
 
opening = weirs[['gml_id','DOORSTROOM','LAAGSTEDOO', 'geometry']]
opening = opening.rename(columns={'gml_id':'stuwid',
                                  'DOORSTROOM':'laagstedoorstroombreedte',
                                  'LAAGSTEDOO':'laagstedoorstroomhoogte'})
opening['globalid'] = np.nan
for i in range(len(opening['globalid'])):
    opening['globalid'].iloc[i] = 'HHRL_opening_'+str(weirs['ORACLE_OBJ'].iloc[i])

management_device = opening[['globalid','geometry']]
management_device = management_device.rename(columns={'globalid':'kunstwerkopeningid'})
management_device['overlaatonderlaat'] = 'overlaat'
# %%
#### GEMALEN ####
filename = "\Gemaal\gemaal_peil.shp"
pumpingstations = gpd.read_file(folder_data+filename)

pumpingstationssub = pumpingstations[['CODE','geometry']]
pumpingstationssub['globalid'] = daf.getuniquecode('HHRL_gemaal_',len(pumpingstations['CODE']))
pumpingstationssub = pumpingstationssub.rename(columns={'CODE':'code'})

pumps = pumpingstations[['ORACLE_OBJ','MAXIMALECA','geometry']]
pumps = pumps.rename(columns={'MAXIMALECA':'maximalecapaciteit',
                              'ORACLE_OBJ':'code'})
pumps['gemaalid'] = pumpingstationssub['globalid']
pumps['globalid'] = daf.getuniquecode('HDSR_pomp_',len(pumps['maximalecapaciteit']))

sturing = pumpingstations[['streefpeil','geometry']]
sturing = sturing.rename(columns={ 'streefpeil':'streefwaarde'})
sturing['pompid'] = pumps['globalid']
sturing['ondergrens'] = sturing['streefwaarde'] - 0.05
sturing['bovengrens'] = sturing['streefwaarde'] + 0.05
sturing['doelvariabele'] = 'waterstand'
sturing['code'] = daf.getuniquecode('HDSR_sturing_', len(sturing['pompid']))
sturing['globalid'] = daf.getuniquecode('HDSR_sturing_glob_',len(sturing['pompid']))

# Sla de verschillende kunstwerken en watergangen op in een geopackage per waterschap
file_gpkg = "D:\work\P1414_ROI\GIS\HHRijnland\HHRL_hydamo.gpkg"

hydroobjectsub.to_file(file_gpkg, layer='waterloop', driver="GPKG")
bridgessub.to_file(file_gpkg, layer='brug', driver="GPKG")
culverts_total.to_file(file_gpkg, layer='duiker', driver="GPKG")
weirssub.to_file(file_gpkg, layer='stuw', driver="GPKG")
opening.to_file(file_gpkg, layer='kunstwerkopening', driver="GPKG")
management_device.to_file(file_gpkg, layer='regelmiddel', driver="GPKG")
pumpingstationssub.to_file(file_gpkg, layer='gemaal', driver="GPKG")
pumps.to_file(file_gpkg, layer='pomp', driver="GPKG")
sturing.to_file(file_gpkg, layer='sturing', driver="GPKG")
# %%
