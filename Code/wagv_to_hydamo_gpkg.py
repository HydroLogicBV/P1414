#%%
# Dit script zet WAGV data om naar een geschikt geopackage format voor HyDAMO.
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
p_folder = r"D:\work\P1414_ROI"
folder_data = p_folder + r"\GIS\WAGV"
file_gpkg = p_folder + r"\GIS\WAGV\WAGV_hydamo.gpkg"

#### WATERGANGEN ####
filename = p_folder + r"\GIS\Uitgesneden watergangen\AGV_v4_test.shp"

# filename_test = p_folder + r"\GIS\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
# hydroobject = gpd.read_file(filename_test)
# hydroobjectsub = hydroobject[["CODE", "GLOBALID", "LENGTE", "geometry"]]
# hydroobjectsub = hydroobjectsub.rename(
#     columns={"CODE": "code", "GLOBALID": "globalid", "LENGTE": "lengte"}
# )

filename_test = r"D:\work\P1414_ROI\GIS\AGV\norm_profielen_test.gpkg"
hydroobject = gpd.read_file(filename)#, layer="hydroobject")
hydroobjectsub = hydroobject[["code", "ruwheidsty", 'ruwheidlaa',"geometry"]]

hydroobjectsub = hydroobjectsub.rename(columns={'ruwheidsty': 'typeruwheid',
                                                'ruwheidlaa':'ruwheid'})

hydroobjectsub['globalid'] = daf.getuniquecode('WAGV_watergang_', len(hydroobjectsub['geometry']))

#### BRUGGEN ####
filename = r"\brug_v13\brug_v13_clipped.shp"
bridges = gpd.read_file(folder_data + filename)
bridgessub = bridges[["code", "lengte", 'ruwheidsty','ruwheid',
                      'intreeverl', 'uittreever', "geometry"]]
bridgessub = bridgessub.rename(
    columns={'ruwheidsty':'typeruwheid',
             'intreeverl':'intreeverlies',
             'uittreever':'uittreeverlies'}
)
bridgessub['globalid'] = daf.getuniquecode('WAGV_brug_glob_', len(bridgessub['geometry']))
bridgessub.dropna(axis=0, inplace=True, subset=["geometry"])

#### DUIKERS ####
filename = r"\duikersifonhevel_v13\duikersifonhevel_v13_clipped.shp"
culverts = gpd.read_file(folder_data + filename)
culvertssub = culverts[
    [
        "code",
        "vormkokeri",
        "lengte",
        "hoogteopen",
        "breedteope",
        "hoogtebinn",
        "hoogtebin0",
        "ruwheidsty",
        "ruwheid",
        "intreeverl",
        "uittreever",
        "geometry",
    ]
]
culvertssub = culvertssub.rename(
    columns={
        "hoogtebinn": "hoogtebinnenonderkantbene",
        "hoogtebin0": "hoogtebinnenonderkantbov",
        "hoogteopen": "hoogteopening",
        "breedteope": "breedteopening",
        "vormkokeri": "vormkoker",
        "ruwheidsty": "typeruwheid",
        "intreeverl": "intreeverlies",
        "uittreever": "uittreeverlies",
    }
)
culvertssub['globalid'] = daf.getuniquecode('WAGV_duiker_', len(culvertssub['geometry']))

#### STUWEN ####
filename = r"\stuw_v13\stuw_v13_clipped.shp"
weirs = gpd.read_file(folder_data + filename)

weirssub = weirs[["code", "soortstuwi", "afvoercoef","geometry"]]
weirssub = weirssub.rename(columns={"code": "globalid", 
                                    "soortstuwi": "soortstuw",
                                    "afvoercoef": "afvoercoefficient"})
weirssub['code'] = daf.getuniquecode('WAGV_stuw_',len(weirssub['geometry']))


filename = "\doorstroomopening_stuw_v13\doorstroomopening_stuw_v13_clipped.shp"
opening_total = gpd.read_file(folder_data + filename)
opening = opening_total[["stuwcode", "laagstedoo", "laagstedo0", "hoogstedoo",
                         "hoogstedo0", "afvoercoef","vormopenin","geometry"]]
opening = opening.rename(
    columns={
        "stuwcode": "stuwid",
        "laagstedoo": "laagstedoorstroombreedte",
        "laagstedo0": "laagstedoorstroomhoogte",
        "hoogstedo0": "hoogstedoorstroomhoogte",
        "afvoercoef": "afvoercoefficient",
    }
)
opening["globalid"] = daf.getuniquecode("HDSR_OPEN_", len(weirssub["globalid"]))


management_device = opening[["globalid", "geometry"]]
management_device = management_device.rename(columns={"globalid": "kunstwerkopeningid"})
management_device["soortregelbaarheid"] = weirs["soortregel"]
management_device["code"] = daf.getuniquecode("HDSR_regelmiddel_", len(management_device["geometry"]))
management_device["overlaatonderlaat"] = "Overlaat"


#### GEMALEN ####
filename = r"\pomp_gemaal_v13\pomp_gemaal_v13_clipped_streefpeil.shp"
pumpingstations = gpd.read_file(folder_data + filename)

pumpingstationssub = pumpingstations[["gemaalcode", "code", "geometry"]]
pumpingstationssub = pumpingstationssub.rename(columns={"gemaalcode": "globalid"})
pumpingstationssub = pumpingstationssub.to_crs("epsg:28992")

pumps = pumpingstations[["gemaalcode", "CODE", "MAXIMALECA", "geometry"]]
pumps = pumps.rename(
    columns={"gemaalcode": "gemaalid", "MAXIMALECA": "maximalecapaciteit", "CODE": "globalid"}
)
pumps["code"] = daf.getuniquecode("HDSR_pomp_", len(pumps["gemaalid"]))
pumps = pumps.to_crs("epsg:28992")

sturing = pumpingstations[["CODE", "streefpeil", "geometry"]]
sturing = sturing.rename(columns={"CODE": "pompid", "streefpeil": "streefwaarde"})
sturing["ondergrens"] = sturing["streefwaarde"] - 0.05
sturing["bovengrens"] = sturing["streefwaarde"] + 0.05
sturing["doelvariabele"] = "waterstand"
sturing["code"] = daf.getuniquecode("HDSR_sturing_", len(sturing["pompid"]))
sturing["globalid"] = daf.getuniquecode("HDSR_sturing_glob_", len(sturing["pompid"]))
sturing = sturing.to_crs("epsg:28992")

# Sla de verschillende kunstwerken en watergangen op in een geopackage per waterschap


hydroobjectsub.to_file(file_gpkg, layer="waterloop", driver="GPKG")
bridgessub.to_file(file_gpkg, layer="brug", driver="GPKG")
culvertssub.to_file(file_gpkg, layer="duiker", driver="GPKG")
weirssub.to_file(file_gpkg, layer="stuw", driver="GPKG")
opening.to_file(file_gpkg, layer="kunstwerkopening", driver="GPKG")
management_device.to_file(file_gpkg, layer="regelmiddel", driver="GPKG")
pumpingstationssub.to_file(file_gpkg, layer="gemaal", driver="GPKG")
pumps.to_file(file_gpkg, layer="pomp", driver="GPKG")
sturing.to_file(file_gpkg, layer="sturing", driver="GPKG")

# %%
