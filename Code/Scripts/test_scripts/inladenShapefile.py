# inladen geofile

import geopandas as gpd

# Load shapefile
gdf = gpd.read_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HDSR\HydroObject\HydroObject_v8.shp")

# Print column names and CRS
print("Columns:", gdf.columns)
print("CRS:", gdf.crs)

gdf.to_crs('EPSG:28992')

# Make the adjustment to the opnemen column.

names = ['hdsr_wl_H040953',
        'hdsr_wl_H074394',
        'hdsr_wl_H015860',
        'hdsr_wl_H077402',
        'hdsr_wl_H077401',
        'hdsr_wl_H074412',
        'hdsr_wl_H072474',
        'hdsr_wl_H014605',
        'hdsr_wl_H000017',
        'hdsr_wl_H013080',
        'hdsr_wl_H076583',
        'hdsr_wl_H077310',
        'hdsr_wl_H012909',
        'hdsr_wl_H012198',
        'hdsr_wl_H012642',
        'hdsr_wl_H074394',
        'hdsr_wl_H015860',
        'hdsr_wl_H000325',
        'hdsr_wl_H000326',
        'hdsr_wl_H074391',
        'hdsr_wl_21f30d35',
        'hdsr_wl_390b1467',
        'hdsr_wl_H076209',
        'hdsr_wl_H000351',
        'hdsr_wl_H012658',
        'hdsr_wl_H012897',
        'hdsr_wl_H022088',
        ]  # Replace with your actual list

names_add = [
    'H055566',
    'H040103',
    'H043136',
]

names_add = [
    'H077130'
]

# Strip 'hdsr_wl_' from each name
stripped_names = [name.replace('hdsr_wl_', '') for name in names]

# Set the desired values to NEE
#gdf.loc[gdf['CODE'].isin(stripped_names), 'OPNEMEN'] = 'NEE'
#gdf.loc[gdf['CODE'].isin(names_add), 'OPNEMEN'] = 'JA'
gdf.loc[gdf['CODE'].isin(names_add), 'OPNEMEN'] = 'JA'

# Save to a new shapefile
gdf.to_file(r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HDSR\HydroObject\HydroObject_v9.shp")

print("Shapefile export complete.")

"""
import geopandas as gpd
from shapely.ops import unary_union

# Load layers
polygons = gpd.read_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\peilen\boezem_v2.shp")
lines = gpd.read_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\ARKNZK\netwerk_selectie_v4.shp")

lines.set_crs('epsg:28992')

polygons.to_crs('epsg:28992')
#lines.to_crs('epsg:28992')

# Ensure same CRS
#if polygons.crs != lines.crs:
#    lines = lines.to_crs(polygons.crs)

# Filter polygons that intersect with the lines
intersecting_polygons = polygons[polygons.geometry.apply(lambda poly: lines.geometry.unary_union.intersects(poly))]

# Merge the intersecting polygons
merged_polygon = intersecting_polygons.geometry.unary_union

# Create a GeoDataFrame for the result
merged_gdf = gpd.GeoDataFrame({"geometry": [merged_polygon]}, crs=polygons.crs)

# Save the result to a shapefile
merged_gdf.to_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\ARKNZK\netwerk_selectie_v4_polygons_merged.shp")

print("Processing complete. Merged shapefile saved.")
"""