import xarray as xr
import pyproj
import geopandas as gpd
import os
import pandas as pd
import time
from shapely.geometry import Point
import re
from scipy.spatial import KDTree
import numpy as np
import os
import shutil
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import json 
from shapely.geometry import Point, LineString
import geopandas as gpd
from ipyleaflet import Map, Marker, projections, basemaps, GeoData
import ipyleaflet as ipl
import pyproj
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points



class NetcdfNetwork():
    def __init__(self, location, crs_rd, crs_map, rd_to_map, map_to_rd):
        self.location = location
        self.crs_rd = crs_rd
        self.crs_map = crs_map
        self.rd_to_map = rd_to_map
        self.map_to_rd = map_to_rd

        self.generate_geodfs()

    def generate_geodfs(self):
        ds = xr.open_dataset(self.location)
        network = [(ds.network_node_x.values[i], ds.network_node_y.values[i]) for i in range(len(ds.network_node_x.values))]
        mesh= [(ds.Mesh2d_face_x.values[i], ds.Mesh2d_face_y.values[i]) for i in range(len(ds.Mesh2d_face_x.values))]

        df = pd.DataFrame()
        df['coors'] = network
        df['geometry'] = df['coors'].apply(Point)
        self.network_rd = gpd.GeoDataFrame(df['geometry'], crs = self.crs_rd)
        self.network_map = self.network_rd.to_crs(self.crs_map)

        df = pd.DataFrame()
        df['coors'] = mesh
        df['geometry'] = df['coors'].apply(Point)
        self.mesh_rd = gpd.GeoDataFrame(df['geometry'], crs = self.crs_rd) 
        self.mesh_map = self.mesh_rd.to_crs(self.crs_map)
    




class DambreakWidget():
    def __init__(self, model_folder):
        self.path_keringen = r"C:\Werk\Projecten\P1414_ROI\Github_P1414\P1414\Code\notebooks\data\Combined_test_v14_WBD.gpkg"
        self.network_loc = os.path.join(model_folder, r'dflowfm\network.nc')
        
        self.crs_rd = 28992
        self.crs_map = 4326
        crs_rd = pyproj.CRS.from_epsg(self.crs_rd)
        crs_map = pyproj.CRS.from_epsg(self.crs_map)

        # Create a transformer object to convert between the two coordinate systems
        self.map_to_rd = pyproj.Transformer.from_crs(crs_map, crs_rd)
        self.rd_to_map = pyproj.Transformer.from_crs(crs_rd, crs_map)

        self.netcdf = NetcdfNetwork(self.network_loc, self.crs_rd, self.crs_map, self.rd_to_map, self.map_to_rd)

        self.current_step = 0

    def gdf_to_geojson(self, input_gdf, layer_name):
        json_data = json.loads(input_gdf.to_json())
        geojson = ipl.GeoJSON(data=json_data, name = layer_name)
        return geojson
    
    def filter_gdf_based_on_location(self, location, dataframe):
        dist_max_x = 0.05
        dist_max_y = 0.03
        dataframe = dataframe[(location[1] - dist_max_x < dataframe.geometry.x) & (dataframe.geometry.x < location[1] + dist_max_x)]
        dataframe = dataframe.loc[(location[0] - dist_max_y < dataframe.geometry.y) & (dataframe.geometry.y < location[0] + dist_max_y)]
        return dataframe

    def draw_map(self, center = [51.970682, 4.64013599]):
        if self.current_step <= 0:
            self.m = Map(center=center, zoom=12)
        
        if self.current_step <= 1:
            
            network_1d_points = GeoData(geo_dataframe = network,
                style={'color': 'black', 'radius':2, 'fillColor': '#3366cc', 'opacity':0.5, 'weight':1.9, 'dashArray':'2', 'fillOpacity':0.6},
                hover_style={'fillColor': 'red' , 'fillOpacity': 0.2},
                point_style={'radius': 1, 'color': 'red', 'fillOpacity': 0.8, 'fillColor': 'blue', 'weight': 1},
                name = 'NETWORK')
        
            self.m.add_layer(network_1d_points)
    
            mesh2d = GeoData(geo_dataframe = mesh,
                style={'color': 'black', 'radius':2, 'fillColor': '#8366cc', 'opacity':0.5, 'weight':1.9, 'dashArray':'2', 'fillOpacity':0.6},
                hover_style={'fillColor': 'red' , 'fillOpacity': 0.2},
                point_style={'radius': 1, 'color': 'red', 'fillOpacity': 0.8, 'fillColor': 'blue', 'weight': 1},
                name = 'MESH')

            self.m.add_layer(mesh2d)
           
        # self.m.add_layer(self.gdf_to_geojson(self.netcdf.mesh_map.iloc[0:100], 'test'))
        display(self.m)


    def remove_layer_from_map(self):
        for index, layer in enumerate(self.m.layers):
            if layer.name == 'MESH':
                 self.m.remove_layer(self.m.layers[index])
                
    def add_network(self):
        network = self.filter_gdf_based_on_location(self.m.center, self.netcdf.network_map)
        network_1d_points = GeoData(geo_dataframe = network,
            style={'color': 'black', 'radius':2, 'fillColor': '#3366cc', 'opacity':0.5, 'weight':1.9, 'dashArray':'2', 'fillOpacity':0.6},
            hover_style={'fillColor': 'red' , 'fillOpacity': 0.2},
            point_style={'radius': 1, 'color': 'red', 'fillOpacity': 0.8, 'fillColor': 'blue', 'weight': 1},
            name = 'NETWORK')
        self.m.add_layer(network_1d_points)
    
    def add_mesh(self):
        mesh = self.filter_gdf_based_on_location(self.m.center, self.netcdf.mesh_map)
        mesh2d = GeoData(geo_dataframe = mesh,
            style={'color': 'black', 'radius':2, 'fillColor': '#8366cc', 'opacity':0.5, 'weight':1.9, 'dashArray':'2', 'fillOpacity':0.6},
            hover_style={'fillColor': 'red' , 'fillOpacity': 0.2},
            point_style={'radius': 1, 'color': 'red', 'fillOpacity': 0.8, 'fillColor': 'blue', 'weight': 1},
            name = 'MESH')
        self.m.add_layer(mesh2d)
 
    
    def next_step(self):
        pass
        



    




