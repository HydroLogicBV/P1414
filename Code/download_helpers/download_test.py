# %%
from download_helpers.utilities import json_download_rijnland, wfs_download_rijnland

# %% Stuw
stuw_url = r"https://rijnland.enl-mcs.nl/arcgis/services/Stuw/MapServer/WFSServer?service=wfs&request=GetFeature&typeNames=stuw"
stuw_output_path = r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Stuw\stuw_v2.shp"

wfs_download_rijnland(url=stuw_url, output_file_path=stuw_output_path)

# %% gemaal
gemaal_url = r"https://rijnland.enl-mcs.nl/arcgis/rest/services/Gemaal/MapServer/0/query?"
gemaal_output_path = r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Gemaal\gemaal.shp"

json_download_rijnland(url=gemaal_url, output_file_path=gemaal_output_path)

# %% brug
brug_url = r"https://rijnland.enl-mcs.nl/arcgis/services/Brug/MapServer/WFSServer?service=wfs&request=GetFeature&typeNames=brug"
brug_output_path = r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Brug\brug.shp"

wfs_download_rijnland(url=brug_url, output_file_path=brug_output_path)

# %% duiker
duiker_url = r"https://rijnland.enl-mcs.nl/arcgis/rest/services/Duiker/MapServer/0/query?"
duiker_output_path = r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Duiker\duiker.shp"

json_download_rijnland(url=duiker_url, output_file_path=duiker_output_path)

# %% sifon
sifon_url = r"https://rijnland.enl-mcs.nl/arcgis/services/Sifon/MapServer/WFSServer?service=wfs&request=GetFeature&typeNames=sifon"
sifon_output_path = r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Sifon\sifon.shp"

wfs_download_rijnland(url=sifon_url, output_file_path=sifon_output_path)

# %% keringen AGV
kering_agv_url = r"https://maps.waternet.nl/arcgis/rest/services/AGV_Legger/AGV_Legger_Keringen_GEONIS/MapServer?f=jsapi"
kering_agv_output_path = r"D:\Work\Project\P1414\GIS\WAGV\kering\kering.shp"

json_download_rijnland(url=kering_agv_url, output_file_path=kering_agv_output_path)

# %%