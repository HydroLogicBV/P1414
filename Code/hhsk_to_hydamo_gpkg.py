#%%
# Dit script zet HH Schieland en de Krimpenerwaard data om naar een geschikt geopackage format voor HyDAMO.
# Stappen die genomen worden
# - Voor verschillende elementen gegevens inladen en ordenen naar HyDAMO formaat
# - Opslaan als nieuwe .shp file
# - Samenvoegen in een gpkg (opslagtechnisch en overzichtelijker)

# Importeer de benodige packages
import geopandas as gpd
import numpy as np
import pandas as pd
import data_functions as daf
from importlib import reload
reload(daf)
# Definieer locatie waar bestanden staan 
folder_data = r"D:\work\P1414_ROI\GIS\HHSK\Legger"

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
filename = r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHSK.shp"
#filename = r"D:\work\P1414_ROI\GIS\HHRijnland\Legger\Watergang\Watergang_as.shp"
hydroobject = gpd.read_file(filename)

hydroobjectsub = hydroobject[['CODE', 'GlobalID','LENGTE','geometry']]
hydroobjectsub = hydroobjectsub.rename(columns={'GlobalID':'globalid',
                                                'LENGTE':'lengte',
                                                'CODE':'code',
                                                })
hydroobjectsub['ruwheid'] = 23.
hydroobjectsub['typeruwheid'] = 4

#### BRUGGEN ####

# Tot nu toe geen bruggen gevonden in de legger van HHSK: niet meegenomen

#### DUIKERS ####
filename_duiker = "\Duiker.shp"

culverts = gpd.read_file(folder_data + filename_duiker)

culvertssub = culverts[['CODE','VORM','LENGTE', 'GlobalID',
                        'HOOGTE','BREEDTE','HOOGTEBOK',
                        'geometry']]
culvertssub = culvertssub.rename(columns={'HOOGTEBOK':'hoogtebinnenonderkantbene',
                                          'CODE':'code',
                                          'HOOGTE':'hoogteopening',
                                          'BREEDTE':'breedteopening',
                                          'VORM':'vormkoker',
                                          'LENGTE':'lengte'})

culvertssub['hoogtebinnenonderkantbov'] = culvertssub['hoogtebinnenonderkantbene']
culvertssub['typeruwheid'] = 4
culvertssub['ruwheid'] = 75.
culvertssub['intreeverlies'] = 0.6
culvertssub['uittreeverlies'] = 0.8

# Check if vormkoker is an int or a str
if type(culvertssub.vormkoker.iloc[0]) == str:
    culvertssub.vormkoker = daf.vormkoker_str2int(culvertssub.vormkoker)

#### STUWEN ####
filename = "\Stuw.shp"
weirs = gpd.read_file(folder_data+filename)

weirssub = weirs[['CODE','GlobalID','geometry']]
weirssub = weirssub.rename(columns={'GlobalID':'globalid',
                                    'CODE':'code',
                                    })
weirssub['afvoercoefficient'] = 0.85
weirssub['soortstuw'] = 'overlaat'

opening = weirs[['GlobalID','KRUINBREED','KRUINHOOGT','geometry']]
opening = opening.rename(columns={'GlobalID':'stuwid',
                                  'KRUINBREED':'laagstedoorstroombreedte',
                                  'KRUINHOOGT':'laagstedoorstroomhoogte',})
opening['hoogstedoorstroombreedte'] = opening['laagstedoorstroombreedte']
opening['hoogstedoorstroomhoogte'] = opening['laagstedoorstroomhoogte']
opening['vormopening'] = 3

opening['globalid'] = np.nan
for i in range(len(opening['globalid'])):
    opening['globalid'].iloc[i] = 'HHSK_opening_'+str(weirs['OBJECTID'].iloc[i])

management_device = opening[['globalid','geometry']]
management_device = management_device.rename(columns={'globalid':'kunstwerkopeningid'})
management_device['overlaatonderlaat'] = 'overlaat'

#### GEMALEN ####
filename = "\Gemaal_peil.shp"
pumpingstations = gpd.read_file(folder_data+filename)

pumpingstationssub = pumpingstations[['CODE','GlobalID','geometry']]
pumpingstationssub = pumpingstationssub.rename(columns={'CODE':'code',
                                                        'GlobalID':'globalid'})

pumps = pumpingstations[['OBJECTID','CAPACITEIT','geometry']]
pumps = pumps.rename(columns={'CAPACITEIT':'maximalecapaciteit',
                              'OBJECTID':'code'})
pumps['gemaalid'] = pumpingstationssub['globalid']
pumps['globalid'] = daf.getuniquecode('HHSK_pomp_',len(pumps['maximalecapaciteit']))

sturing = pumpingstations[['streefpeil','geometry']]
sturing = sturing.rename(columns={ 'streefpeil':'streefwaarde'})
sturing['pompid'] = pumps['globalid']
sturing['ondergrens'] = sturing['streefwaarde'] - 0.05
sturing['bovengrens'] = sturing['streefwaarde'] + 0.05
sturing['doelvariabele'] = 'waterstand'
sturing['code'] = daf.getuniquecode('HHSK_sturing_', len(sturing['pompid']))
sturing['globalid'] = daf.getuniquecode('HHSK_sturing_glob_',len(sturing['pompid']))

# Sla de verschillende kunstwerken en watergangen op in een geopackage per waterschap
file_gpkg = "D:\work\P1414_ROI\GIS\HHSK\HHSK_hydamo.gpkg"

hydroobjectsub.to_file(file_gpkg, layer='waterloop', driver="GPKG")
#bridgessub.to_file(file_gpkg, layer='brug', driver="GPKG")
culvertssub.to_file(file_gpkg, layer='duiker', driver="GPKG")
weirssub.to_file(file_gpkg, layer='stuw', driver="GPKG")
opening.to_file(file_gpkg, layer='kunstwerkopening', driver="GPKG")
management_device.to_file(file_gpkg, layer='regelmiddel', driver="GPKG")
pumpingstationssub.to_file(file_gpkg, layer='gemaal', driver="GPKG")
pumps.to_file(file_gpkg, layer='pomp', driver="GPKG")
sturing.to_file(file_gpkg, layer='sturing', driver="GPKG")
# %%
