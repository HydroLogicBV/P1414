#%%
# Dit script zet HDSR data om naar een geschikt geopackage format voor HyDAMO.
# Stappen die genomen worden
# - Voor verschillende elementen gegevens inladen en ordenen naar HyDAMO formaat
# - Opslaan als nieuwe .shp file
# - Samenvoegen in een gpkg (opslagtechnisch en overzichtelijker)

# Importeer de benodige packages
import uuid
from copy import copy

import geopandas as gpd

from HDSR_norm_profiles import hdsr_norm_profiles

# Definieer locatie waar bestanden staan
# p_folder = r"D:\work\P1414_ROI"
p_folder = r"D:\Work\Project\P1414"
branches_path = p_folder + r"\GIS\HDSR\hydro_object_w_norm_profielen.gpkg"
bridges_path = p_folder + r"\GIS\HDSR\Legger\Bruggen\Bruggen.shp"
culvert_path = p_folder + r"\GIS\HDSR\Legger\Kokers_Lijnen\Kokers_Lijnen.shp"
pump_path = p_folder + r"\GIS\HDSR\Legger\Gemalen\Gemalen_peil.shp"
weir_path = p_folder + r"\GIS\HDSR\Legger\Stuwen\BR_Stuwen.shp"

output_gpkg = p_folder + r"\GIS\HDSR\HDSR_hydamo.gpkg"

#### WATERGANGEN ####
# First try to find geopacke with branches and norm-profiles, if it does not exist, create it

try:
    hydroobject = gpd.read_file(branches_path, layer="hydroobject")
except (FileNotFoundError, ValueError):
    # set paths
    old_branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HDSR_v4_test.shp"

    # set column mapping
    index_mapping = dict(
        [
            ("bodembreedte", "IWS_W_BODB"),
            ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
            ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
            ("hoogte insteek linkerzijde", "IWS_W_INST"),
            ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
            ("taludhelling linkerzijde", "IWS_W_TALU"),
            ("taludhelling rechterzijde", "IWS_W_TA_1"),
            ("typeruwheid", "ruwheidsty"),
            ("ruwheidhoog", "ruwheidhoo"),
            ("ruwheidlaag", "ruwheidlaa"),
            ("water_width_index", "IWS_W_WATB"),
        ]
    )
    hdsr_norm_profiles(input_path=old_branches_path, index_mapping=index_mapping, output_path=branches_path)
    hydroobject = gpd.read_file(branches_path, layer="hydroobject")


hydroobjectsub = copy(hydroobject[["code", "globalid", "geometry"]])

hydroobjectsub["typeruwheid"] = 6
hydroobjectsub["ruwheid"] = 23.0

hydroobjectsub = hydroobjectsub.to_crs("epsg:28992")

## normprofielen ##
hydroobject_normgp = gpd.read_file(branches_path, layer="hydroobject_normgp")
normgeparamprofielwaarde = gpd.read_file(branches_path, layer="normgeparamprofielwaarde")


#### BRUGGEN ####
bridges = gpd.read_file(bridges_path)
bridgessub = copy(bridges[["OBJECTID", "CODE", "WS_DOORVAA", "geometry"]])
bridgessub = bridgessub.rename(
    columns={"OBJECTID": "globalid", "CODE": "code", "WS_DOORVAA": "lengte"}
)
bridgessub["globalid"] = [str(uuid.uuid4()) for _ in range(bridgessub.shape[0])]
bridgessub["typeruwheid"] = 4
bridgessub["ruwheid"] = 75.0
bridgessub["intreeverlies"] = 0.5
bridgessub["uittreeverlies"] = 0.7

bridgessub.dropna(axis=0, inplace=True, subset=["geometry"])

bridgessub = bridgessub.to_crs("epsg:28992")

#### DUIKERS ####
culverts = gpd.read_file(culvert_path)
culvertssub = copy(
    culverts[
        [
            "OBJECTID",
            "CODE",
            "VORMKOKER",
            "LENGTE",
            "HOOGTEOPEN",
            "BREEDTEOPE",
            "HOOGTEBOKB",
            "HOOGTEBO_1",
            "geometry",
        ]
    ]
)
culvertssub = culvertssub.rename(
    columns={
        "HOOGTEBOKB": "hoogtebinnenonderkantbene",
        "HOOGTEBO_1": "hoogtebinnenonderkantbov",
        "OBJECTID": "globalid",
        "CODE": "code",
        "HOOGTEOPEN": "hoogteopening",
        "BREEDTEOPE": "breedteopening",
        "VORMKOKER": "vormkoker",
        "LENGTE": "lengte",
    }
)
culvertssub["globalid"] = [str(uuid.uuid4()) for _ in range(culvertssub.shape[0])]
culvertssub["typeruwheid"] = 4
culvertssub["ruwheid"] = 75.0
culvertssub["intreeverlies"] = 0.6
culvertssub["uittreeverlies"] = 0.8
culvertssub = culvertssub.to_crs("epsg:28992")

#### STUWEN ####
weirs = gpd.read_file(weir_path)
weirs["globalid"] = [str(uuid.uuid4()) for _ in range(weirs.shape[0])]

weirssub = copy(weirs[["CODE", "globalid", "SOORTSTUW", "geometry"]])
weirssub = weirssub.rename(columns={ "SOORTSTUW": "soortstuw", "CODE": "code"})

weirssub["afvoercoefficient"] = 1.0
weirssub = weirssub.to_crs("epsg:28992")

opening = copy(weirs[["globalid", "DOORSTROOM", "LAAGSTEDOO", "HOOGSTEDOO", "geometry"]])
opening = opening.rename(
    columns={
        "globalid": "stuwid",
        "DOORSTROOM": "laagstedoorstroombreedte",
        "LAAGSTEDOO": "laagstedoorstroomhoogte",
        "HOOGSTEDOO": "hoogstedoorstroomhoogte",
    }
)
opening["globalid"] = [str(uuid.uuid4()) for _ in range(opening.shape[0])]
opening["hoogstedoorstroombreedte"] = opening["laagstedoorstroombreedte"]
opening["vormopening"] = 3
opening["afvoercoefficient"] = 0.85
opening["geometry"] = None
# opening = opening.to_crs("epsg:28992")


management_device = copy(opening[["globalid", "geometry", "stuwid"]])
management_device = management_device.rename(columns={"globalid": "kunstwerkopeningid"})
management_device["globalid"] = [str(uuid.uuid4()) for _ in range(management_device.shape[0])]
management_device["code"] = management_device["globalid"]
management_device["soortregelbaarheid"] = weirs["SOORTREGEL"]
management_device["overlaatonderlaat"] = "Overlaat"
# management_device = management_device.to_crs("epsg:28992")

#### GEMALEN ####
pumpingstations = gpd.read_file(pump_path)

pumpingstationssub = copy(pumpingstations[["GEMAALID", "geometry"]])
pumpingstationssub["globalid"] = [str(uuid.uuid4()) for _ in range(pumpingstationssub.shape[0])]
pumpingstationssub = pumpingstationssub.rename(columns={"GEMAALID": "code"})
pumpingstationssub = pumpingstationssub.astype({"code": str})
pumpingstationssub = pumpingstationssub.to_crs("epsg:28992")

pumps = copy(pumpingstations[["MAXIMALECA", "geometry", "GEMAALID"]])
pumps = pumps.rename(columns={"MAXIMALECA": "maximalecapaciteit", "GEMAALID": "code"})
pumps = pumps.astype({"code": str})
pumps["gemaalid"] = pumpingstationssub["globalid"]
pumps["globalid"] = [str(uuid.uuid4()) for _ in range(pumps.shape[0])]
# pumps["code"] = pumps["globalid"]
pumps["geometry"] = None
# pumps = pumps.to_crs("epsg:28992")

sturing = copy(pumpingstations[["streefpeil", "geometry"]])
sturing = sturing.rename(columns={"streefpeil": "streefwaarde"})
sturing["pompid"] = pumps["globalid"]
sturing["ondergrens"] = sturing["streefwaarde"] - 0.05
sturing["bovengrens"] = sturing["streefwaarde"] + 0.05
sturing["doelvariabele"] = "waterstand"
sturing["code"] = sturing["pompid"]
sturing = sturing.astype({"code": str})
sturing["globalid"] = [str(uuid.uuid4()) for _ in range(sturing.shape[0])]
sturing["geometry"] = None
# sturing = sturing.to_crs("epsg:28992")

# Sla de verschillende kunstwerken en watergangen op in een geopackage per waterschap
hydroobjectsub.to_file(output_gpkg, layer="waterloop", driver="GPKG")
bridgessub.to_file(output_gpkg, layer="brug", driver="GPKG")
culvertssub.to_file(output_gpkg, layer="duiker", driver="GPKG")
weirssub.to_file(output_gpkg, layer="stuw", driver="GPKG")
opening.to_file(output_gpkg, layer="kunstwerkopening", driver="GPKG")
management_device.to_file(output_gpkg, layer="regelmiddel", driver="GPKG")
pumpingstationssub.to_file(output_gpkg, layer="gemaal", driver="GPKG")
pumps.to_file(output_gpkg, layer="pomp", driver="GPKG")
sturing.to_file(output_gpkg, layer="sturing", driver="GPKG")
hydroobject_normgp.to_file(output_gpkg, layer="hydroobject_normgp", driver="GPKG")
normgeparamprofielwaarde.to_file(output_gpkg, layer="normgeparamprofielwaarde", driver="GPKG")

# %%
 