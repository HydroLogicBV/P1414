#%%
# Dit script zet HDSR data om naar een geschikt geopackage format voor HyDAMO.
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
folder_data = r"D:\work\P1414_ROI\GIS\HDSR\Legger"

#### WATERGANGEN ####
filename = "\Hydro_Objecten(2)\HydroObject.shp"
hydroobject = gpd.read_file(folder_data + filename)
hydroobjectsub = hydroobject[['CODE','GLOBALID','LENGTE','geometry']]
hydroobjectsub = hydroobjectsub.rename(columns={'CODE':'code',
                                                'GLOBALID':'globalid',
                                                'LENGTE':'lengte'})
#### BRUGGEN ####
filename = "\Bruggen\Bruggen.shp"
bridges = gpd.read_file(folder_data+filename)
bridgessub = bridges[['OBJECTID','CODE','WS_DOORVAA','geometry']]
bridgessub = bridgessub.rename(columns={'OBJECTID':'globalid',
                                        'CODE':'code',
                                        'WS_DOORVAA':'lengte'})
bridgessub['typeruwheid'] = 4
bridgessub['ruwheid'] = 75.
bridgessub['intreeverlies'] = 0.5
bridgessub['uittreeverlies'] = 0.7

#### DUIKERS ####
filename = "\Kokers_Lijnen\Kokers_Lijnen.shp"
culverts = gpd.read_file(folder_data + filename) 
culvertssub = culverts[['OBJECTID','CODE','VORMKOKER','LENGTE',
                        'HOOGTEOPEN','BREEDTEOPE','HOOGTEBOKB',
                        'HOOGTEBO_1','geometry']]
culvertssub = culvertssub.rename(columns={'HOOGTEBOKB':'hoogtebinnenonderkantbene',
                                          'HOOGTEBO_1':'hoogtebinnenonderkantbov',
                                          'OBJECTID':'globalid',
                                          'CODE':'code',
                                          'HOOGTEOPEN':'hoogteopening',
                                          'BREEDTEOPEN':'breedteopening',
                                          'VORMKOKER':'vormkoker',
                                          'LENGTE':'lengte'})
culvertssub['typeruwheid'] = 4
culvertssub['ruwheid'] = 75.
culvertssub['intreeverlies'] = 0.6
culvertssub['uittreeverlies'] = 0.8

#### STUWEN ####
filename = "\Stuwen\BR_Stuwen.shp"
weirs = gpd.read_file(folder_data+filename)

weirssub = weirs[['CODE','STUWID','geometry']]
weirssub = weirssub.rename(columns={'STUWID':'globalid'})
weirssub['afvoercoefficient'] = 1.
 
opening = weirs[['STUWID','DOORSTROOM','LAAGSTEDOO', 'geometry']]
opening = opening.rename(columns={'STUWID':'stuwid',
                                  'DOORSTROOM':'laagstedoorstroombreedte',
                                  'LAAGSTEDOO':'laagstedoorstroomhoogte'})
opening['globalid'] = daf.getuniquecode('HDSR_OPEN_', len(weirssub['globalid']))

management_device = opening[['globalid','geometry']]
management_device = management_device.rename(columns={'globalid':'kunstwerkopeningid'})
management_device['overlaatonderlaat'] = 'overlaat'

#### GEMALEN ####
filename = "\Gemalen\Gemalen_peil.shp"
pumpingstations = gpd.read_file(folder_data+filename)

pumpingstationssub = pumpingstations[['GEMAALID','geometry']]
pumpingstationssub['code'] = daf.getuniquecode('HDSR_gemaal_',len(pumpingstations['GEMAALID']))
pumpingstationssub = pumpingstationssub.rename(columns={'GEMAALID':'globalid'})

pumps = pumpingstations[['GEMAALID','CODE','MAXIMALECA','geometry']]
pumps = pumps.rename(columns={'GEMAALID':'gemaalid',
                              'MAXIMALECA':'maximalecapaciteit',
                              'CODE':'globalid'})
pumps['code'] = daf.getuniquecode('HDSR_pomp_',len(pumps['gemaalid']))

sturing = pumpingstations[['CODE','streefpeil','geometry']]
sturing = sturing.rename(columns={ 'CODE':'pompid',
                                   'streefpeil':'streefwaarde'})
sturing['ondergrens'] = sturing['streefwaarde'] - 0.05
sturing['bovengrens'] = sturing['streefwaarde'] + 0.05
sturing['doelvariabele'] = 'waterstand'

# Sla de verschillende kunstwerken en watergangen op in een geopackage per waterschap
file_gpkg = "D:\work\P1414_ROI\GIS\HDSR\HDSR_hydamo.gpkg"

hydroobjectsub.to_file(file_gpkg, layer='waterloop', driver="GPKG")
bridgessub.to_file(file_gpkg, layer='brug', driver="GPKG")
culvertssub.to_file(file_gpkg, layer='duiker', driver="GPKG")
weirssub.to_file(file_gpkg, layer='stuw', driver="GPKG")
opening.to_file(file_gpkg, layer='kunstwerkopening', driver="GPKG")
management_device.to_file(file_gpkg, layer='regelmiddel', driver="GPKG")
pumpingstationssub.to_file(file_gpkg, layer='gemaal', driver="GPKG")
pumps.to_file(file_gpkg, layer='pomp', driver="GPKG")
sturing.to_file(file_gpkg, layer='sturing', driver="GPKG")
