# Specify here the folder (location) in which the GIS folder is located
folder_path_GIS = r"D:\Work\Project\P1414"

class RawData:
    ## PATHS
    p_folder = folder_path_GIS + r"\GIS"
    weir_path = folder_path_GIS + r"\GIS\Keringen\stuwen.shp"

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
