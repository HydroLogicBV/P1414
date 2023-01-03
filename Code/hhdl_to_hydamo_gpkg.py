#%%
# Dit script zet HH Delfland data om naar een geschikt geopackage format voor HyDAMO.
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
folder_data = r"D:\work\P1414_ROI\GIS\HHDelfland\Legger_Delfland_shp"

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
filename = r"D:\work\P1414_ROI\GIS\Uitgesneden watergangen\HHD.shp"
#filename = r"D:\work\P1414_ROI\GIS\HHRijnland\Legger\Watergang\Watergang_as.shp"
hydroobject = gpd.read_file(filename)

hydroobjectsub = hydroobject[['CODE','LENGTE','geometry']]
hydroobjectsub = hydroobjectsub.rename(columns={'LENGTE':'lengte',
                                                'CODE':'globalid',
                                                })
hydroobjectsub['code'] = daf.getuniquecode('HHDL_watergang_',len(hydroobjectsub['lengte']))                                                
hydroobjectsub['ruwheid'] = 23.
hydroobjectsub['typeruwheid'] = 4

#### BRUGGEN ####

# Tot nu toe geen bruggen gevonden in de legger van HHDL: niet meegenomen

#### DUIKERS ####
filename_duiker_list = ["\Ondersteunende kunstwerken\Open duiker.shp",
                        "\Ondersteunende kunstwerken\Sifon.shp",
                        "\Ondersteunende kunstwerken\Stuwende duiker.shp",
                        "\Ondersteunende kunstwerken\Vispassageduiker.shp",
                        "\Ondersteunende kunstwerken\Inlaatduiker.shp"]

for filename_duiker in filename_duiker_list:
    culverts = gpd.read_file(folder_data + filename_duiker)
    if filename_duiker != filename_duiker_list[0]:
        culvertssub = culverts[['CODE','VORMKOKER_','LENGTE',
                                'HOOGTEOPEN','BREEDTEOPE',
                                'geometry']]
        culvertssub = culvertssub.rename(columns={'CODE':'code',
                                                'HOOGTEOPEN':'hoogteopening',
                                                'BREEDTEOPE':'breedteopening',
                                                'VORMKOKER_':'vormkoker',
                                                'LENGTE':'lengte'})
        if 'HOOGTEBINN' in culverts.columns:
            culvertssub['hoogtebinnenonderkantbov'] = culverts['HOOGTEBINN']
            culvertssub['hoogtebinnenonderkantbene'] = culverts['HOOGTEBI00']
        else: 
            culvertssub['hoogtebinnenonderkantbov'] = None #TODO FIX 
            culvertssub['hoogtebinnenonderkantbene'] = None #TODO FIX    
      
        culvert_total = culvert_total.append(culvertssub)
    else: 
        culvertssub = culverts[['CODE','VORM','LENGTE',
                                'HOOGTEOPEN','BREEDTEOPE','HOOGBOKBOV',
                                'HOOGBOKBEN','geometry']]
        culvertssub = culvertssub.rename(columns={'CODE':'code',
                                                'HOOGTEOPEN':'hoogteopening',
                                                'BREEDTEOPE':'breedteopening',
                                                'VORM':'vormkoker',
                                                'LENGTE':'lengte',
                                                'HOOGBOKBOV':'hoogtebinnenonderkantbov',
                                                'HOOGBOKBEN': 'hoogtebinnenonderkantbene',
                                                })        
        culvert_total = culvertssub
        

culvert_total['typeruwheid'] = 4
culvert_total['ruwheid'] = 75.
culvert_total['intreeverlies'] = 0.6
culvert_total['uittreeverlies'] = 0.8
culvert_total['globalid'] = daf.getuniquecode('HHRL_duiker_', len(culvert_total['geometry']))

# Check if vormkoker is an int or a str
if type(culvert_total.vormkoker.iloc[0]) == str:
    culvert_total.vormkoker = daf.vormkoker_str2int(culvert_total.vormkoker)

#### STUWEN ####
filename = "\Ondersteunende kunstwerken\Stuw.shp"
weirs = gpd.read_file(folder_data+filename)

weirssub = weirs[['CODE','geometry']]
weirssub = weirssub.rename(columns={'CODE':'code',
                                    })
weirssub['afvoercoefficient'] = 1
weirssub['soortstuw'] = 'overlaat'
weirssub['globalid'] = daf.getuniquecode('HHDL_stuw_',len(weirssub['geometry']))

opening = weirs[['DOORSTROOM','LAAGSTEDOO','HOOGSTEDOO','geometry']]
opening = opening.rename(columns={'GlobalID':'stuwid',
                                  'DOORSTROOM':'laagstedoorstroombreedte',
                                  'LAAGSTEDOO':'laagstedoorstroomhoogte',
                                  'HOOGSTEDOO':'hoogstedoorstroomhoogte',})
opening['hoogstedoorstroombreedte'] = opening['laagstedoorstroombreedte']
opening['vormopening'] = 3
opening['stuwid'] = weirssub['globalid']
opening['afvoercoefficient'] = 0.85

opening['globalid'] = daf.getuniquecode('HHDL_opening_',len(opening['geometry']))

management_device = opening[['globalid','geometry']]
management_device = management_device.rename(columns={'globalid':'kunstwerkopeningid'})
management_device['code'] = daf.getuniquecode('HHDL_regelmiddel_',len(management_device['geometry']))
management_device["soortregelbaarheid"] = weirs['REGELBAARH']
management_device['overlaatonderlaat'] = 'Overlaat'

#### GEMALEN ####
filename = "\Ondersteunende kunstwerken\Gemaal_peil.shp"
pumpingstations = gpd.read_file(folder_data+filename)

pumpingstationssub = pumpingstations[['CODE','geometry']]
pumpingstationssub = pumpingstationssub.rename(columns={'CODE':'code',})
pumpingstationssub['globalid'] = daf.getuniquecode('HHDL_gemaal_',len(pumpingstationssub['geometry']))

pumps = pumpingstations[['MAXCAPACIT','geometry']]
pumps = pumps.rename(columns={'MAXCAPACIT':'maximalecapaciteit'})
pumps['gemaalid'] = pumpingstationssub['globalid']
pumps['code'] = daf.getuniquecode('HHDL_pomp_',len(pumps['geometry']))
pumps['globalid'] = daf.getuniquecode('HHDL_pomp_globid_',len(pumps['geometry']))

sturing = pumpingstations[['streefpeil','geometry']]
sturing = sturing.rename(columns={ 'streefpeil':'streefwaarde'})
sturing['pompid'] = pumps['globalid']
sturing['ondergrens'] = sturing['streefwaarde'] - 0.05
sturing['bovengrens'] = sturing['streefwaarde'] + 0.05
sturing['doelvariabele'] = 'waterstand'
sturing['code'] = daf.getuniquecode('HHDL_sturing_', len(sturing['pompid']))
sturing['globalid'] = daf.getuniquecode('HHDL_sturing_glob_',len(sturing['pompid']))

# Sla de verschillende kunstwerken en watergangen op in een geopackage per waterschap
file_gpkg = "D:\work\P1414_ROI\GIS\HHDelfland\HHDL_hydamo.gpkg"

hydroobjectsub.to_file(file_gpkg, layer='waterloop', driver="GPKG")
#bridgessub.to_file(file_gpkg, layer='brug', driver="GPKG")
culvert_total.to_file(file_gpkg, layer='duiker', driver="GPKG")
weirssub.to_file(file_gpkg, layer='stuw', driver="GPKG")
opening.to_file(file_gpkg, layer='kunstwerkopening', driver="GPKG")
management_device.to_file(file_gpkg, layer='regelmiddel', driver="GPKG")
pumpingstationssub.to_file(file_gpkg, layer='gemaal', driver="GPKG")
pumps.to_file(file_gpkg, layer='pomp', driver="GPKG")
sturing.to_file(file_gpkg, layer='sturing', driver="GPKG")
# %%
