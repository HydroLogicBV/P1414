# inladen geofile
"""
import geopandas as gpd

# Load shapefile
gdf = gpd.read_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HDSR\HydroObject\HydroObject_v3.shp")

# Print column names and CRS
print("Columns:", gdf.columns)
print("CRS:", gdf.crs)

gdf.to_crs('EPSG:28992')

# Save to a new shapefile
gdf.to_file(r"\\dc02\Project\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HDSR\HydroObject\HydroObject_v4.shp")

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