import xarray as xr
import pyproj
import geopandas as gpd
import os
import pandas as pd
import time
from shapely.geometry import Point, LineString, Polygon
import re
from scipy.spatial import KDTree
import numpy as np
import os
import shutil
import ipywidgets as ipy
from IPython.display import display, clear_output, Image
import json 
from shapely.geometry import Point, LineString
import geopandas as gpd
from ipyleaflet import Map, Marker, projections, basemaps, GeoData, ImageOverlay, Popup
from ipywidgets import HTML
import ipyleaflet as ipl
import pyproj
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
from branca.colormap import linear

class NetcdfNetwork():
    def __init__(self, location, crs_rd, crs_map, rd_to_map, map_to_rd):
        self.location = location
        self.crs_rd = crs_rd
        self.crs_map = crs_map
        self.rd_to_map = rd_to_map
        self.map_to_rd = map_to_rd
        self.generate_geodfs()
    
    def generate_mesh2d(self, point:Point, max_distance:int):

        """Creates both a geodataframe with rectangles that reperesent the 2d mesh, and two geodataframes that contain all points of the 2d mesh (filtered around breach location) 

        Args:
            point (Point): Central point around which to draw rectangles
            max_distance (int): maximum distance

        Returns:
            gpd.GeoDataFrame: GeoDataFrame with the rectangles and bed level as attribute
        """
        ds = xr.open_dataset(self.location)
        distance_x = ds.Mesh2d_face_x - point.x
        distance_y = ds.Mesh2d_face_y - point.y
        distance = np.sqrt(np.square(distance_x) + np.square(distance_y)) 
        max_rectangles = 1000
        max_distance_treshold = np.partition(distance, max_rectangles)[max_rectangles]
        treshold = min(max_distance, max_distance_treshold)
        mask =  distance < treshold

        if mask.sum() == 0:
            raise Exception("You selected a point that is not near the 2D grid of the model.")

        mesh2d = np.column_stack([ds.Mesh2d_face_x.values[mask], ds.Mesh2d_face_y.values[mask]]) 
        df = pd.DataFrame()
        df['coors'] = [tuple(row) for row in mesh2d]
        df['geometry'] = df['coors'].apply(Point)
        self.mesh2d_rd = gpd.GeoDataFrame(df['geometry'], crs = self.crs_rd) 
        self.mesh2d_map = self.mesh2d_rd.to_crs(self.crs_map)

        Mesh2d_node_x = ds.Mesh2d_node_x.values
        Mesh2d_node_y = ds.Mesh2d_node_y.values
        Mesh2d_node_z = ds.Mesh2d_node_z.values
        bedlevel = []
        geometry = []
        offset = np.min(ds.Mesh2d_face_nodes.values[ds.Mesh2d_face_nodes.values > -1])
        for face in ds.Mesh2d_face_nodes.values[mask]:
            coords = []
            node_heights = []
            for node in face:
                if np.isnan(node):
                    continue
                if node < 0:
                    continue
                node_index = node - offset
                node_heights.append(Mesh2d_node_z[int(node_index)])
                coords.append(
                    [
                        Mesh2d_node_x[int(node_index)], 
                        Mesh2d_node_y[int(node_index)]
                    ]
                )
            bedlevel.append(np.average(node_heights))
            geometry.append(Polygon(coords))
        ds.close()

        self.rectangles = gpd.GeoDataFrame({'bedlevel': bedlevel}, geometry=geometry, crs=self.crs_rd)

    def generate_geodfs(self):
        """Reads the netcdf file, and create geodataframes of the data tha tis expected to be used.

        Raises:
            Exception: _description_
        """
        ds = xr.open_dataset(self.location)
        network = [(ds.network_node_x.values[i], ds.network_node_y.values[i]) for i in range(len(ds.network_node_x.values))]
        if 'mesh1d_node_x' not in ds.variables.keys():
            raise Exception("Mesh coordinates are not in network.nc, try importing and exporting your model in the D-HYDRO GUI.")
        mesh1d = np.column_stack([ds.mesh1d_node_x.values, ds.mesh1d_node_y.values])
        mesh2d = np.column_stack([ds.Mesh2d_face_x.values, ds.Mesh2d_face_y.values])       

        # find the link ids, these have different names, so look for name that contains link and id
        for key in list(ds.keys()):
            if ('link' in key) and ('id' in key):
                links_ids = ds[key].values
                links_ids = [l.decode("utf-8") for l in links_ids]
                links = np.array([[int(l.split('_')[0]), int(l.split('_')[-1])] for l in links_ids])
                break

        df = pd.DataFrame()
        df['coors'] = network
        df['geometry'] = df['coors'].apply(Point)
        self.network_rd = gpd.GeoDataFrame(df['geometry'], crs = self.crs_rd)
        self.network_map = self.network_rd.to_crs(self.crs_map)
        self.network_node_id = ds.network_node_id.values  

        df = pd.DataFrame()
        df['coors'] = [tuple(row) for row in mesh1d]
        df['geometry'] = df['coors'].apply(Point)
        self.mesh1d_rd = gpd.GeoDataFrame(df['geometry'], crs = self.crs_rd) 
        self.mesh1d_map = self.mesh1d_rd.to_crs(self.crs_map)

        self.generate_mesh2d(Point([103895, 440297]), 1000)

        links_geo = list(np.stack((mesh1d[links[:, 0].astype(int)], mesh2d[links[:, 1].astype(int)]), axis=1)) # van link (index - index) naar linestring (point - point)
        line_strings = [LineString(line) for line in links_geo]
        self.links_rd = gpd.GeoDataFrame({'geometry':line_strings}, crs = self.crs_rd)
        self.links_map = self.links_rd.to_crs(self.crs_map)
        ds.close()

class DambreakWidget(WidgetStyling):
    def __init__(self, model_folder):
        self.network_loc = os.path.join(model_folder, r'dflowfm\network.nc')
        self.model_path = model_folder
        self.max_distance_around_breach = 2000
        self.bedlevel_layer = None
        
        self.set_default_layout_and_styling()
        
        self.crs_rd = 28992
        self.crs_map = 4326
        crs_rd = pyproj.CRS.from_epsg(self.crs_rd)
        crs_map = pyproj.CRS.from_epsg(self.crs_map)

        # Create a transformer object to convert between the two coordinate systems
        self.map_to_rd = pyproj.Transformer.from_crs(crs_map, crs_rd)
        self.rd_to_map = pyproj.Transformer.from_crs(crs_rd, crs_map)

        self.netcdf = NetcdfNetwork(self.network_loc, self.crs_rd, self.crs_map, self.rd_to_map, self.map_to_rd)
        self.keringen = gpd.read_file(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'data\keringen.shp'), crs = self.crs_rd)

        self.center = [51.90915, 4.74685]
        self.zoom = 12
        self.marker = Marker(location=self.center, draggable=True, name="marker")
        self.current_step = 0
        self.display_text = ""
        self.draw_control = ipl.DrawControl()
        self.draw_control.polyline =  {
            "shapeOptions": {
                "color": "grey",
                "weight": 3,
                "opacity": 1.0
            }
        }
        self.draw_control.polygon = {}
        self.draw_control.circlemarker = {}
        
        url_legend = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'data\legend.png')
        image = Image(url_legend, width = 300)

        self.bed_level_viewer = ipy.HTML("")

        self.image_widget = ipy.Image(
            value= image.data,
            format='png', 
            width=300,
            )

    def gdf_to_geojson(self, input_gdf:gpd.GeoDataFrame, layer_name:str) -> ipl.GeoJSON:
        json_data = json.loads(input_gdf.to_json())
        geojson = ipl.GeoJSON(data=json_data, name = layer_name)
        return geojson

    def draw_map(self):
        clear_output()
        if self.current_step == 0:
            instruction = "Step 0: Select the approximate location of your dambreach (0/4) "
            self.m = Map(center=self.center, zoom=self.zoom) 
            self.m.add(self.marker)
            display(self.m)
        elif self.current_step == 1:
            self.completed_line_input = False
            instruction = "Step 1: Draw dambreach location (1/4) "
            self.m = Map(center=self.center, zoom=self.zoom) 
            self.draw_control.on_draw(self.handle_draw)    
            self.m.add(self.draw_control)
            display(self.m)
            self.m.add(self.bedlevel_layer)
            self.add_keringen()
            self.add_links()
        elif self.current_step == 2:
            instruction = "Step 2: Select upstream 1D computational node (2/4) "
            self.m = Map(center=self.center, zoom=self.zoom)
            display(self.m)
            self.add_network()
            self.m.add(self.bedlevel_layer)
            self.m.add(self.marker)
            self.m.add(self.breach_point_map)
        elif self.current_step == 3:
            instruction = "Step 3: Select downstream 2D grid node (3/4)" 
            self.m = Map(center=self.center, zoom=self.zoom)
            display(self.m)
            self.add_mesh()
            self.m.add(self.bedlevel_layer)
            self.m.add(self.marker)
            self.m.add(self.upstream)
            self.m.add(self.breach_point_map)
        elif self.current_step == 4:
            instruction = "Done! Inspect your results (4/4) "
            self.m = Map(center=self.center, zoom=self.zoom)
            display(self.m)
            self.add_mesh()
            self.add_network()
            self.m.add(self.upstream)
            self.m.add(self.downstream)
            self.m.add(self.dambreach_map)
            self.m.add(self.breach_point_map)
            self.m.add(self.breach_line_map)
            self.add_keringen()
            self.add_links()
            
        
        if self.layer_exists('bedlevel'):
            self.bed_level_viewer.value = "<i>Select a grid cell to view the height...</i>"
            display(self.bed_level_viewer)

        html = ipy.HTML(value=f'<b style="color:black;font-size:18px;">{instruction}</b>')
        self.html_log = ipy.HTML(value='')
        display(html)  
        display(self.html_log)
        if self.current_step < 4:
            button_next_step = ipy.Button(description="Next step", style = self.button_style, layout = self.button_layout)
            output = ipy.Output()
            display(button_next_step, output)
            button_next_step.on_click(self.next_step)
        display(self.image_widget)

    def next_step(self, _):
        self.center = self.m.center 
        self.zoom = self.m.zoom
        if self.current_step == 0:
            self.current_step = 1
            self.focus_point_rd = Point(self.map_to_rd.transform(self.marker.location[0], self.marker.location[1]))
            self.keringen = self.filter_gdf(self.keringen, self.focus_point_rd, self.max_distance_around_breach)
            self.netcdf.links_rd = self.filter_gdf(self.netcdf.links_rd, self.focus_point_rd, self.max_distance_around_breach)
            self.netcdf.generate_mesh2d(self.focus_point_rd, self.max_distance_around_breach)
            self.create_bedlevel_layer()
            self.draw_map() 
        elif self.current_step == 1:
            if self.completed_line_input == True:
                self.current_step = 2
                self.draw_map()    
            else:
                self.html_log.value = '<b style="color:red;font-size:14px;">You do not have a valid dambreak yet, so you can not continue</b>'
        elif self.current_step == 2:
            self.current_step = 3
            self.index_closest_1d = self.snap_to_closest_point(self.marker.location, self.netcdf.mesh1d_map, self.netcdf.mesh1d_rd)
            self.upstream_point = self.netcdf.mesh1d_rd.iloc[self.index_closest_1d:self.index_closest_1d+1]
            self.upstream = GeoData(geo_dataframe = self.upstream_point.to_crs(self.crs_map),
                        point_style={'radius': 5, 'color': 'blue', 'fillOpacity': 1, 'fillColor': 'white', 'weight': 3},
                    name = 'UPSTREAM')
            self.draw_map() 
        elif self.current_step == 3:
            self.current_step = 4
            self.index_closest_2d = self.snap_to_closest_point(self.marker.location, self.netcdf.mesh2d_map, self.netcdf.mesh2d_rd)
            self.downstream_point = self.netcdf.mesh2d_rd.iloc[self.index_closest_2d:self.index_closest_2d+1]
            self.downstream = GeoData(geo_dataframe = self.downstream_point.to_crs(self.crs_map),
                    point_style={'radius': 5, 'color': 'green', 'fillOpacity': 1, 'fillColor': 'white', 'weight': 3},
                    name = 'DOWNSTREAM')
            self.calculate_breach()
            self.draw_map() 
            self.settings = {}
            self.settings['downstream'] = self.downstream_point
            self.settings['upstream'] = self.upstream_point
            self.settings['breach'] = self.breach_point
            self.settings['breach_line'] = self.breach_line_rd
        elif self.current_step == 4:
            self.draw_map()

    def handle_draw(self, target, action, geo_json:dict):
        self.completed_line_input = False
        self.html_log.value = ""
        self.draw_control.clear()
        self.remove_layer_from_map('new_line')
        self.remove_layer_from_map('BREACH')
        geo_json['properties']['style']['color'] = '#6bc2e5'
        geo_json['properties']['style']['weight'] = 4
        geo_json['geometry']['coordinates'] = geo_json['geometry']['coordinates'][0:2]
        new_line = ipl.GeoJSON(data=geo_json, name = 'new_line')
        self.m.add(new_line)
        out = self.snap_to_kering(geo_json['geometry']['coordinates'])
        if out is not None:
            self.m.add(self.breach_point_map)
            self.completed_line_input = True
            breach_row = self.breach_point.to_crs(self.crs_map).iloc[0]
            self.center = breach_row.geometry.y, breach_row.geometry.x
            self.marker.location = breach_row.geometry.y, breach_row.geometry.x
    
    def filter_gdf(self, gdf:gpd.GeoDataFrame, point:Point, buffer_distance:int) -> gpd.GeoDataFrame:
        buffer = point.buffer(buffer_distance)
        gdf = gdf[gdf.geometry.intersects(buffer)]
        return gdf

    def snap_to_kering(self, line):
        line_gdf = gpd.GeoDataFrame(geometry = [LineString(line)], crs = self.crs_map).to_crs(self.crs_rd)
        line_dambreak = line_gdf.iloc[0].geometry
        intersection_keringen = self.find_intersection(self.keringen, line_dambreak)
        if intersection_keringen == False:
            self.remove_layer_from_map('new_line')
            self.remove_layer_from_map('BREACH')
            self.html_log.value = '<b style="color:red;font-size:14px;">Your line does not intersect with a kering</b>'
            return None
    
        intersection_lins = self.find_intersection(self.netcdf.links_rd, line_dambreak)
        if intersection_lins is None:
            self.remove_layer_from_map('new_line')
            self.remove_layer_from_map('BREACH')
            self.html_log.value = '<b style="color:red;font-size:14px;">Your line does not intersect with a 1D2D link</b>'
            return None

        snap_rd = Point(intersection_keringen.x, intersection_keringen.y)
        self.breach_line_rd = line_gdf
        self.breach_line_map = GeoData(geo_dataframe = self.breach_line_rd.to_crs(self.crs_map),
                    style={'color': 'grey', 'radius':10, 'fillColor': 'grey', 'opacity':1, 'weight':2, 'dashArray':'1', 'fillOpacity':1},
                    name = 'BREACH_LINE')  
        self.breach_point = gpd.GeoDataFrame(geometry = [snap_rd], crs = self.crs_rd)
        self.breach_point_map = GeoData(geo_dataframe = self.breach_point.to_crs(self.crs_map),
                point_style={'radius': 5, 'color': 'red', 'fillOpacity': 1, 'fillColor': 'white', 'weight': 3},
                name = 'BREACH')
        return 0

    def find_intersection(self, lines:gpd.GeoDataFrame, line_dambreak:LineString) -> Point:
        center = Point(np.average(line_dambreak.coords.xy[0]), np.average(line_dambreak.coords.xy[1]))
        lines['distance_to_center'] = lines.distance(center)
        lines_ordered_by_distance = lines.sort_values(by='distance_to_center')
        found_intersection = False
        for index, row in lines_ordered_by_distance.iterrows():
            if row.geometry.intersects(line_dambreak):
                int_pt = row.geometry.intersection(line_dambreak)
                found_intersection = True
                break
        if found_intersection:
            return int_pt
        return None

    def calculate_breach(self):
        line_1 = LineString([[self.upstream_point.geometry.x, self.upstream_point.geometry.y], [self.breach_point.geometry.x, self.breach_point.geometry.y]])
        line_2 = LineString([[self.breach_point.geometry.x, self.breach_point.geometry.y], [self.downstream_point.geometry.x, self.downstream_point.geometry.y]])
        dambreach = gpd.GeoDataFrame(geometry = [line_1, line_2], crs = self.crs_rd)
        self.dambreach_map = GeoData(geo_dataframe = dambreach.to_crs(self.crs_map),
                    style={'color': 'black', 'radius':10, 'fillColor': 'black', 'opacity':1, 'weight':2, 'dashArray':'2', 'fillOpacity':1},
                    name = 'DAMBREACH_LINES')       

    def remove_layer_from_map(self, name:str):
        to_remove = []
        for index, layer in enumerate(self.m.layers):
            if layer.name == name:
                 to_remove.append(index)
        for index in reversed(to_remove):
            self.m.remove_layer(self.m.layers[index])

    def add_keringen(self):
        keringen_map = GeoData(geo_dataframe = self.keringen.to_crs(self.crs_map),
                style={'color': 'red', 'radius':10, 'fillColor': 'green', 'opacity':1, 'weight':3, 'dashArray':'2', 'fillOpacity':0.6},
                name = 'KERINGEN')
        self.m.add(keringen_map)
    
    def add_links(self):
        self.netcdf.links_map = self.netcdf.links_rd.to_crs(self.netcdf.crs_map)
        links_map = GeoData(geo_dataframe = self.netcdf.links_map,
                style={'color': 'green', 'radius':10, 'fillColor': 'green', 'opacity':1, 'weight':3, 'dashArray':'2', 'fillOpacity':0.6},
                name = 'LINKS')
        self.m.add(links_map)
        
    def add_network(self):
        mesh1d = self.filter_gdf(self.netcdf.mesh1d_rd, self.focus_point_rd, self.max_distance_around_breach)
        network_1d_points = GeoData(geo_dataframe = mesh1d.to_crs(self.crs_map),
            point_style={'radius': 5, 'color': 'blue', 'fillOpacity': 1, 'fillColor': 'black', 'weight': 3},
            name = 'NETWORK')
        self.m.add(network_1d_points)
    
    def add_mesh(self):
        mesh2d = self.filter_gdf(self.netcdf.mesh2d_rd, self.focus_point_rd, self.max_distance_around_breach)
        mesh2d_points = GeoData(geo_dataframe = mesh2d.to_crs(self.crs_map),
            point_style={'radius': 5, 'color': 'green', 'fillOpacity': 1, 'fillColor': 'black', 'weight': 3},
            name = 'MESH')
        self.m.add(mesh2d_points)
    
    def snap_to_closest_point(self, point:list, geodf_map:gpd.GeoDataFrame, geodf_rd:gpd.GeoDataFrame):
        point = self.map_to_rd.transform(point[0], point[1])
        index_closest = geodf_rd.distance(Point(point[0], point[1])).argmin()
        closest_point = geodf_map.iloc[index_closest]
        self.marker.location = (closest_point.geometry.y, closest_point.geometry.x)
        return index_closest

    def create_bedlevel_layer(self):
        def handle_click(feature=None,  **kwargs):
            try:
                bedlevel = round(feature.get('properties').get('bedlevel'),2)
                self.bed_level_viewer.value = f"Bedlevel: {bedlevel} m + NAP"
            except Exception as e:
                self.bed_level_viewer.value = f"Could not retrieve due to error: {e}"

        if self.bedlevel_layer is None:
            rectangles = self.netcdf.rectangles
            rectangles = rectangles.to_crs(self.crs_map)
            bedlevel =  dict(zip(rectangles.index, rectangles.bedlevel))
            bedlevel = {str(key): bedlevel[key] for key in bedlevel.keys()}
            data = json.loads(rectangles.to_json())
            self.bedlevel_layer = ipl.Choropleth(
                geo_data=data,
                choro_data=bedlevel,
                border_color='black',
                colormap=linear.YlOrBr_08,
                style={'fillOpacity': 0.5},
                hover_style={'fillOpacity': 1},
                name="bedlevel"
                )
            self.bedlevel_layer.on_click(handle_click)
        
    def layer_exists(self, layer_to_check:str) -> bool:
        for layer in self.m.layers:
            if layer.name == layer_to_check:
                return True
        return False

class ModifyDambreak(WidgetStyling):
    def __init__(self, model_path, dambreak_settings, keringen, use_widget_dambreak, dambreak_template):
        self.model_path = model_path
        self.structures_path = os.path.join(model_path, 'dflowfm\structures.ini')
    
        self.structure_textfile = self.remove_dambreaks_from_structures()

        self.dambreak_settings = {}
        self.wrote_backup = False
        
        if use_widget_dambreak == False:
            self.dambreak_settings_widget = dambreak_settings
            self.keringen = keringen
            self.add_dambreaks_from_widget(
                self.dambreak_settings_widget['breach_line'],
                self.keringen)
        else:
            self.use_template_dambreak(dambreak_template)
        
        self.set_default_layout_and_styling()

        self.settings_to_modify = ['crestLevelIni', 't0',
                                  'crestLevelMin', 'breachWidthIni', 'f1', 'f2', 'uCrit']
        self.settings_to_modify_names = {
            'crestLevelIni': "Initial crest level of dambreach (crestLevelIni)",
            't0': "Time of dikebreach relative to start of simulation (t0 in hours)", 
            'timeToBreachToMaximumDepth' : "Time to breach maximum depth (timeToBreachToMaximumDepth)",
            'crestLevelMin': "Minimum crest level (crestLevelMin)", 
            'breachWidthIni': "Initial breach width (breachWidthIni)", 
            'f1': "f1 paramter in Verheij–van der Knaap (2002) formula", 
            'f2': "f2 paramter in Verheij–van der Knaap (2002) formula", 
            'uCrit': "uCrit paramter in Verheij–van der Knaap (2002) formula"
        }
        self.settings_in_hours = [ 't0']
        self.widgets = {}
        for setting in self.settings_to_modify:
            self.widgets[setting] = ipy.FloatText(
                value = self.convert_to_sas(setting, float(self.dambreak_settings[setting])),
                description=f'{self.settings_to_modify_names[setting]}:',
                disabled=False
                )
            self.widgets[setting].layout = self.item_layout
            self.widgets[setting].style = self.item_style_wide_description
        
        self.widgets_to_display = [widget for widget in self.widgets.values()]
    
    def convert_to_sas(self, key:str, value):
        if key in self.settings_in_hours:
            return value / 60 / 60
        return value

    def convert_to_mdu(self, key:str, value):
        if key in self.settings_in_hours:
            return value * 60 * 60
        return value

    def remove_dambreaks_from_structures(self):
        # find occurences of dambreak structure type
        pattern = r"type\s*=\s*dambreak"
        matches = {}
        with open(self.structures_path, "r") as f:
            structures_textfile = f.readlines()  

        # find all line that contain structure type = dambreak    
        for i, line in enumerate(structures_textfile):
            if re.search(pattern, line):
                matches[i] = line

        # if there are no structures, return text file as it is
        if len(matches) == 0: 
            return structures_textfile
        
        to_remove_start = [] # keep track of which parts need to be removed
        to_remove_end = [] # keep track of which parts need to be removed
        for match in list(matches.keys()): # for each match of type = dambreak
            found = False
            start_index = match
            pattern = r"\[Structure\]"
            for i in reversed(range(start_index - 6, start_index)):
                if re.search(pattern, structures_textfile[i]):
                    found = True
                    dambreak_settings = structures_textfile[i:]
                    to_remove_start.append(i)
                    break
            if not found:
                raise Exception(f"Error, did not find dambreak structure in structures.ini")
                
            # find end of dambreak strucutre
            for i in range(len(dambreak_settings)):
                if dambreak_settings[i] == '\n' or i == len(dambreak_settings)-1:
                    to_remove_end.append(i+1)
                    break
        
        if len(to_remove_end) != len(to_remove_start):
            raise Exception('error')
        
        for i in reversed(range(len(to_remove_start))):
            structures_textfile = structures_textfile[:to_remove_start[i]] + structures_textfile[to_remove_end[i] + to_remove_start[i]:]

        return structures_textfile

    def should_line_be_flipped(self, coords_line:np.array, point:Point) -> bool:
        """Checks on which side the point lies, returns True if the line has to be flipped, False if not. Needed for D-HYDRO otherwise flow trough dambreak is negative.

        Args:
            coords_line (np.array): coordinates of line (2d array)
            point (Point): point on 2D grid

        Returns:
            bool: Boolean to indicat if the line should be flipped.
        """
        A = coords_line[0]
        B = coords_line[-1]
        P = (point.x, point.y)

        AB = (B[0] - A[0], B[1] - A[1]) 
        AP = (P[0] - A[0], P[1] - A[1])  
        cross_product = AB[0] * AP[1] - AB[1] * AP[0]
        if cross_product > 0:
            return True
        return False
    
    def add_dambreaks_from_widget(self, db_gdf:gpd.GeoDataFrame, fw_gdf:gpd.GeoDataFrame, max_dist:int=100, z_default:int=0):
        fw_point_list = []
        for ix, branch in fw_gdf.iterrows():
            coords = branch.geometry.coords[:]

            for x, y, z in coords:
                fw_point_list.append([x, y, z])

        fw_kdtree = KDTree(data=fw_point_list)

        if len(db_gdf) != 1:
            raise Exception('error')
        
        db = db_gdf.iloc[0]
    
        geom = db.geometry
        if not isinstance(geom, LineString):
            raise ValueError("I need a LineString...")

        coords = np.array(geom.coords[:])
        centroid = geom.interpolate(0.5, normalized=True)
        if len(coords[0]) == 2:
            mid_point = [centroid.x, centroid.y, z_default]
        elif len(coords[0]) == 3:
            mid_point = [centroid.x, centroid.y, centroid.z]
        else:
            raise ValueError("can't cope with {} coordinates".format(len(coords[0])))

        _, ix_point_in_kering = fw_kdtree.query(mid_point, distance_upper_bound=max_dist)
        point_in_kering = fw_point_list[ix_point_in_kering]

        # D-HYDRO needs the hinterland to be at a certain side of the dike breach. Here, check if the downstream point is on correct side, and flip if it is not.
        flip_line = False
        try:
            flip_line = self.should_line_be_flipped(coords, self.dambreak_settings_widget['downstream'].iloc[0].geometry) 
        except:
            print("Flip line function failed")
            pass
        if flip_line:
            coords = np.flip(coords, axis=0)

        self.dambreak_settings['id'] = 'comb_0.0'
        self.dambreak_settings['name'] = 'd3d9f584-086e-4d19-ae01-e75bdca23d21'
        self.dambreak_settings['type'] = 'dambreak'
        self.dambreak_settings['numCoordinates'] = 2
        self.dambreak_settings['xCoordinates'] = ' '.join(str(x) for x in coords[:, 0].tolist())
        self.dambreak_settings['yCoordinates'] = ' '.join(str(x) for x in coords[:, 1].tolist())
        self.dambreak_settings['startLocationX'] = point_in_kering[0]
        self.dambreak_settings['startLocationY'] = point_in_kering[1]
        self.dambreak_settings['algorithm'] = 2
        self.dambreak_settings['crestLevelIni'] = round(point_in_kering[2] - 4, 2) # - db.crestlevelini
        self.dambreak_settings['crestLevelMin'] = -2
        self.dambreak_settings['breachWidthIni'] = 5
        self.dambreak_settings['t0'] = 0
        self.dambreak_settings['timeToBreachToMaximumDepth'] = 360.0
        self.dambreak_settings['f1'] = 1.3
        self.dambreak_settings['f2'] = 0.04
        self.dambreak_settings['uCrit'] = 0.2
        self.dambreak_settings['breachWidthIni'] = 5
        self.dambreak_settings['waterLevelDownstreamLocationX'] = self.dambreak_settings_widget['downstream'].iloc[0].geometry.x
        self.dambreak_settings['waterLevelDownstreamLocationY'] = self.dambreak_settings_widget['downstream'].iloc[0].geometry.y
        # self.dambreak_settings['waterLevelUpstreamLocationX'] = self.dambreak_settings_widget['upstream'].iloc[0].geometry.x
        # self.dambreak_settings['waterLevelUpstreamLocationY'] = self.dambreak_settings_widget['upstream'].iloc[0].geometry.y
    
    def use_template_dambreak(self, dambreak_template:dict):
        self.dambreak_settings['id'] = 'comb_0.0'
        self.dambreak_settings['name'] = 'd3d9f584-086e-4d19-ae01-e75bdca23d21'
        self.dambreak_settings['type'] = 'dambreak'
        self.dambreak_settings['algorithm'] = 2
        self.dambreak_settings['numCoordinates'] = 2
        self.dambreak_settings['crestLevelIni'] = 0 
        self.dambreak_settings['crestLevelMin'] = -2
        self.dambreak_settings['breachWidthIni'] = 5
        self.dambreak_settings['t0'] = 0
        self.dambreak_settings['timeToBreachToMaximumDepth'] = 360.0
        self.dambreak_settings['f1'] = 1.3
        self.dambreak_settings['f2'] = 0.04
        self.dambreak_settings['uCrit'] = 0.2
        self.dambreak_settings['breachWidthIni'] = 5

        for setting in dambreak_template:
            self.dambreak_settings[setting] = dambreak_template[setting]

    def write_to_structures(self, write_output:bool=True, backup_original:bool=True):       
        if write_output:
            if backup_original and self.wrote_backup == False:
                self.wrote_backup = True
                backup_file_loc = os.path.join(os.path.dirname(self.structures_path), 'stuctures_backup.ini')
                shutil.copyfile(self.structures_path, backup_file_loc)

            with open(self.structures_path ,'w') as f:
                for line in self.structure_textfile:
                    f.write(line)
                f.write('\n')
                f.write('[Structure]\n')
                for key, value in self.dambreak_settings.items():
                    f.write(f"    {key.ljust(35)}= {value}\n")
    
    def update_widget_values(self):
        for setting in self.settings_to_modify:
            self.widgets[setting].value = self.convert_to_sas(setting, float(self.dambreak_settings[setting]))
        
    def display_widgets(self):
        self.update_widget_values()
        items = ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style)
        display(items)
        
        button = ipy.Button(description="Update settings", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.update_settings_widget)

    def update_settings_widget(self, _):
        for setting in self.settings_to_modify:
            self.dambreak_settings[setting] = self.convert_to_mdu(setting, self.widgets[setting].value)
        
        self.write_to_structures(write_output = True)

        clear_output(wait=True)
        self.display_widgets()
        display("Dambreak settings are:")
        print_settings_dict = {k: self.dambreak_settings[k] for k in self.settings_to_modify}
        for key in print_settings_dict.keys():
            print_settings_dict[key] = self.convert_to_sas(key, print_settings_dict[key])
        print(json.dumps(print_settings_dict, indent=4))

