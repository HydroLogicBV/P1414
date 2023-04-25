class RawData:
    ## PATHS
    p_folder = r"D:\Work\Project\P1414\GIS"
    weir_path = r"D:\Work\Project\P1414\GIS\Keringen\stuwen.shp"

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
