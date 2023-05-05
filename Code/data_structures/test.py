from ipyleaflet import Map, Marker, projections
import ipyleaflet as ipl
import pyproj as proj
import geopandas as gpd
import json
keringen = gpd.read_file(r"C:\Werk\Projecten\P1414_ROI\Github_P1414\P1414\Code\notebooks\data\Combined_test_v14_WBD.gpkg",
                    layer = 'keringen',
                    crs=28992)
a = 1