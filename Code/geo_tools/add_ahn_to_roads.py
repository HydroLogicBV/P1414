import geopandas as gpd
import networkx as nx
import numpy as np

from add_ahn_height_to_fw import add_height_to_linestrings
from networkx_tools import combine_straight_branches, gdf_to_nx

if __name__ == "__main__":
    ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_merged.tif"
    input_path = r"D:\Work\Project\P1414\GIS\Wegen\Top10NLwegen_en_spoorwegen_clipped.shp"
    output_path = r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\hoofd_en_spoorwegen_v2.shp"
    output_tunnel_path = r"D:\Work\Project\P1414\GIS\Keringen_met_hoogte\tunnel.shp"

    gdf = gpd.read_file(input_path).to_crs(crs="EPSG:28992").explode()

    gdf["fysiekvoor"] = gdf["fysiekvoor"].fillna("")

    brug_bool = gdf["fysiekvoor"].str.contains("brug")
    tunnel_bool = gdf["fysiekvoor"].str.contains("tunnel")

    comb_bool = brug_bool | tunnel_bool

    print(np.sum(comb_bool), np.sum(brug_bool), np.sum(tunnel_bool))

    in_gdf = gdf.loc[~comb_bool, :]
    G = nx.Graph(gdf_to_nx(gdf_network=in_gdf))
    H = combine_straight_branches(G=G)

    out_branches_list = []
    for x, y, data in H.edges(data=True):
        out_branches_list.append(data)

    in_gdf = gpd.GeoDataFrame(data=out_branches_list, geometry="geometry", crs=in_gdf.crs)

    out_gdf = add_height_to_linestrings(gdf=in_gdf, ahn_path=ahn_path)
    out_gdf.to_file(output_path)

    tunnel_gdf = gdf.loc[tunnel_bool, :]
    out_tunnel_gdf = add_height_to_linestrings(gdf=tunnel_gdf, ahn_path=ahn_path)
    out_tunnel_gdf.to_file(output_tunnel_path)
