# inladen geofile

import geopandas as gpd
from shapely.geometry import LineString

# Load shapefile
gdf = gpd.read_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HHSK\Hoofdwatergang\Hoofdwatergangen_model_V6.3.shp")

# Print column names and CRS
print("Columns:", gdf.columns)
print("CRS:", gdf.crs)
#print(gdf)
#gdf.to_crs('EPSG:28992')

gdf.to_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HHSK\Hoofdwatergang\Hoofdwatergangen_model_V6.4.shp")

"""
# Replace the coordinates of branch id 24, index 10
target_geometry = gdf.loc[10, 'geometry']
coords = list(target_geometry.coords)
#print(target_geometry)
#print(coords)
print(f'First coordinate: {target_geometry.coords[0]}')
print(f'Last coordinate: {target_geometry.coords[-1]}')

# coords to map towards:
left_target = gdf.loc[11, 'geometry']
left_coords = list(left_target.coords)
print(f'Left coords: {left_coords}')

# coords to map towards:
right_target = gdf.loc[2, 'geometry']
right_coords = list(right_target.coords)
print(f'Right coords: {right_coords}')

coords_adjusted = [right_coords[0], left_coords[0]]

# Create a new LineString with the adjusted coordinates
new_geometry = LineString(coords_adjusted)

# Update the geometry in the GeoDataFrame
gdf.at[10, 'geometry'] = new_geometry

# _____________________________________________________________________
# Replace the coordinates of branch id 21, index 8
target_geometry = gdf.loc[8, 'geometry']
coords = list(target_geometry.coords)
#print(target_geometry)
#print(coords)
print(f'First coordinate: {target_geometry.coords[0]}')
print(f'Last coordinate: {target_geometry.coords[-1]}')

# coords to map towards:
left_target = gdf.loc[3, 'geometry']
left_coords = list(left_target.coords)
print(f'Left coords: {left_coords}')

# coords to map towards:
right_target = gdf.loc[9, 'geometry']
right_coords = list(right_target.coords)
print(f'Right coords: {right_coords}')

coords_adjusted = [left_coords[-1], right_coords[0]]

# Create a new LineString with the adjusted coordinates
new_geometry = LineString(coords_adjusted)

# Update the geometry in the GeoDataFrame
gdf.at[8, 'geometry'] = new_geometry

gdf.to_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Markermeer\MarkermeerV4.shp")

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