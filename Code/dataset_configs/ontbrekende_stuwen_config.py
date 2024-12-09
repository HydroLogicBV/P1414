import os

# Specify the default path to the GIS folder in the case it is not defined in the environment/main script
default_GIS_path = r"D:\Work\Project\P1414_default"
folder_path_GIS = os.environ.get('GIS_folder_path', default_GIS_path)

class Name:
    name = 'Ontbrekende stuwen'
    
class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    weir_path = folder_path_GIS + r"\GIS\Ontbrekende stuwen in keringen\stuwenV2.shp"

    weir_selection = dict([("column", "OPNEMEN"), ("value", "JA")])

    ## Weirs
    weir_index_mapping = dict(
        [
            ("afvoercoefficient_stuw", None),
            ("afvoercoefficient_opening", None),
            ("code", "id"),
            ("geometry", "geometry"),
            ("globalid", "globalid"),
            ("hoogstedoorstroombreedte", "breedte"),
            ("hoogstedoorstroomhoogte", "hoogte"),
            ("laagstedoorstroombreedte", "breedte"),
            ("laagstedoorstroomhoogte", "hoogte"),
            ("overlaatonderlaat", None),
            ("soortregelbaarheid", None),
            ("soortstuw", None),
            ("vormopening", None),
        ]
    )
