from copy import copy

import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import LineString


def add_height_to_linestrings(gdf: gpd.GeoDataFrame, ahn_path: str) -> gpd.GeoDataFrame:
    ahn = rasterio.open(ahn_path)
    out_gdf = copy(gdf)
    for name, row in gdf.iterrows():
        line = row.geometry
        _line = np.asarray(line.coords)[:, 0:2]
        # np_line = np.insert(_line, 2, values=np.nan, axis=1)

        heights = ahn.sample(_line)
        z = np.array([z for z in heights]).flatten()
        z_mean = np.nanmean(z)
        z[np.isnan(z)] = z_mean

        if np.isnan(z_mean):
            out_gdf.drop(name, inplace=True)
            print("dropped {}".format(name))
            continue

        np_line = np.insert(_line, 2, values=z, axis=1)
        # np_line[:, 2] = [h for h in heights]

        # h_mean = np.nanmean(np_line[:, 2])
        # np_line[np.isnan(np_line[:, 2]), 2] = h_mean

        out_gdf.loc[name, "geometry"] = LineString(np_line)
    ahn.close()
    return out_gdf


if __name__ == "__main__":
    ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_merged.tif"
    # input_paths = [
    #     r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Zeewering landwaartse begrenzing.shp",
    #     r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Regionale waterkering buitenkruinlijn.shp",
    #     r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Overige waterkering polderkade middenkruinlijn.shp",
    #     r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Overige waterkering landscheiding middenkruinlijn.shp",
    #     r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Waterkeringen\Delflandsedijk buitenkruinlijn.shp",
    # ]
    # output_paths = [
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_zeewering.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_regionale_kering.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_polderkade.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_landscheiding.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhd_delflandsedijk.shp"
    # ]
    # input_paths = [
    #     r"D:\Work\Project\P1414\GIS\HHRijnland\Primaire_keringen\Primaire_kering.shp",
    #     r"D:\Work\Project\P1414\GIS\HHRijnland\Regionale_keringen\Regionale_keringen.shp",
    # ]
    # output_paths = [
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhr_primaire_kering.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhr_regional_kering.shp",
    # ]
    # input_paths = [
    #     r"D:\Work\Project\P1414\GIS\HHSK\Keringen\Primaire_waterkering.shp",
    #     r"D:\Work\Project\P1414\GIS\HHSK\Keringen\Regionale_waterkering_1.shp",
    #     r"D:\Work\Project\P1414\GIS\HHSK\Keringen\Overige_waterkering_1.shp",
    # ]
    # output_paths = [
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhsk_primaire_kering.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhsk_regionale_kering.shp",
    #     r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hhsk_overige_kering.shp",
    # ]
    # input_paths=[r"D:\Work\Project\P1414\GIS\WAGV\Keringen\keringen.shp"]
    # output_paths = [ r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\wagv.shp"]
    input_paths = [r"D:\Work\Project\P1414\GIS\Wegen\NBWHoofdwegen.shp"]
    output_paths = [r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hoofdwegen.shp"]

    for input_path, output_path in zip(input_paths, output_paths):
        gdf = gpd.read_file(input_path).to_crs(crs="EPSG:28992").explode()
        print(np.sum(gdf.geometry.type == "MultiLineString"))

        out_gdf = add_height_to_linestrings(gdf=gdf, ahn_path=ahn_path)
        out_gdf.to_file(output_path)
