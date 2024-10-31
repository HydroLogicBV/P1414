import os
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
# from widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import numpy as np
import shutil
from shapely.geometry import Point
from notebooks.background_scripting.v1.boundaryconditions import BoundaryConditionsFile, BoundaryConditionBlock
import geopandas as gpd
import matplotlib.pyplot as plt
from typing import Tuple
import pandas as pd

class ModifyBoundaries(WidgetStyling):
    """
    Class for modifying the boundary conditions of the model
    """
    def __init__(self, model_folder:str, mdu_settings:dict, dambreak_settings:dict):
        self.set_default_layout_and_styling()
        self.crs = 28992
        self.model_folder = model_folder
        files = os.listdir(os.path.join(self.model_folder, 'dflowfm'))
        self.bnd_file, self.bnd_backup_file = None, None
        for file in files:
            if file.endswith('.bc.backup') and 'boundar' in file:
                self.bnd_backup_file = os.path.join(self.model_folder, 'dflowfm', file)
            if file.endswith('.bc') and 'boundar' in file:
                self.bnd_file = os.path.join(self.model_folder, 'dflowfm', file)   
        if self.bnd_file == None:
            raise Exception("Boundary conditions file could not be found in model directory")
        if self.bnd_backup_file == None:
            self.bnd_backup_file = f"{self.bnd_file}.backup"
            shutil.copyfile(self.bnd_file, self.bnd_backup_file)
            
        
        self.bc = BoundaryConditionsFile(self.bnd_backup_file)

        self.home_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.boundary_areas = gpd.read_file(os.path.join(self.home_dir, 'data', 'BoundaryConditions.shp'))

        self.assign_groups_to_boundaries()
        self.create_modify_groups()

        self.t_start =  int(mdu_settings['tStart'])
        self.t_stop = int(mdu_settings['tStop'])
        self.t_breach = int(dambreak_settings['t0'])

    
    def assign_groups_to_boundaries(self):
        for block in self.bc.blocks:
            if block.type.lower() != "[forcing]":
                continue
            if block.location is None:
                continue
            intersection = self.boundary_areas.intersects(block.location)
            names = self.boundary_areas[intersection]['Name']
            if len(names) > 0:
                block.boundary_group = names.iloc[0]
    
    def create_modify_groups(self):
        self.modify_groups = {}
        for block in self.bc.blocks:
            if block.type.lower() != "[forcing]":
                continue
            if block.get_property('function').lower() == 'qhtable':
                continue

            if block.location is not None and block.boundary_group is None:
                continue # if it is a location boundary, but not assigned to a group, do not add them.

            if block.get_property('function').lower() == 'constant':
                block.transform_from_constant_to_timeseries()
            
            if block.boundary_group is not None:
                group_name = block.boundary_group
            else:
                group_name = block.get_property('name')
            
            if group_name in self.modify_groups.keys():
                self.modify_groups[group_name].add_boundary(block)
            else:
                self.modify_groups[group_name] = BoundaryModificationGroup(group_name, block.get_property('quantity', get_match_nr=2), self)
                self.modify_groups[group_name].add_boundary(block)
        

    def generate_discharge_timeseries(self, basic_discharge:int, peak_discharge:int, peak_offset:int, peak_duration:int, **kwargs) -> Tuple[list, list]:
        peak_offset = peak_offset * 3600
        peak_duration = peak_duration * 3600
        if self.t_start % 300 != 0 or self.t_stop % 300 != 0:
            raise Exception("tstart or tstop is not valid")
        times = []
        values = []
        for t in range(self.t_start, self.t_stop + 300, 300):
            times.append(t/3600)
            Q_basic = basic_discharge
            Q_peak = 0
            if t < peak_offset + peak_duration and t > peak_offset:
                Q_peak =  (peak_discharge-basic_discharge) * (1 + np.cos(np.pi + 2*np.pi*((t-peak_offset)/peak_duration))) * 1/2
            Q_total = Q_basic + Q_peak
            values.append(Q_total)
        return times, values
    
    def generate_waterlevel_timeseries(self, basic_waterlevel:float, tidal_offset:float, tidal_range:float, peak_waterlevel_addition:float, peak_offset:float, peak_duration:float, **kwargs) -> Tuple[list, list]:
        tidal_offset = tidal_offset * 3600
        peak_offset = peak_offset * 3600
        peak_duration = peak_duration * 3600
        if self.t_start % 300 != 0 or self.t_stop % 300 != 0:
            raise Exception("tstart or tstop is not valid")
        tidal_duration = 12 * 60 * 60 + 25 * 60
        times = []
        values = []
        for t in range(self.t_start, self.t_stop + 300, 300):
            times.append(t/3600)
            wl_tide = basic_waterlevel + (tidal_range * np.cos(np.pi + 2 * np.pi * ((t - tidal_offset) / tidal_duration)))
            wl_peak = 0
            if t < peak_offset + peak_duration and t > peak_offset:
                wl_peak = peak_waterlevel_addition * (1 + np.cos(np.pi + 2 * np.pi * ((t - peak_offset) / peak_duration))) * 1/2
            wl = wl_tide + wl_peak
            values.append(wl)
        
        return times, values

    def plot_timeseries(self, times:list, values:list, title:str, quantity: str, unit:str, location:str):
        """
        Generate a plot of the discharge timeseries, together with tstart, tstop and tbreach.
        """
        fig, ax = plt.subplots(figsize = (13, 5))
        ax.plot(times, values, color = '#3587A4', lw = 3, label = f'{quantity} {location}')
        range = max(times) - min(times) 
        ax.set_xlim([min(times) - range * 0.01, max(times) + range * 0.01])
        ax.set_ylim(ax.get_ylim())
        ax.set_title(title)
        ax.set_xlabel("T (hours)")
        ax.set_ylabel(f"{quantity} ({unit})")
        ax.plot([self.t_start/3600, self.t_start/3600], [-10000, 100000], color = 'orange', label = "Start of simulation")
        ax.plot([self.t_stop/3600, self.t_stop/3600], [-10000, 100000], color = 'orange', label = "End of simulation")
        ax.plot([self.t_breach/3600, self.t_breach/3600], [-10000, 100000], color = 'red', ls='--', label = "Timestep dikebreach")
        ax.legend()
    
    def display_widgets(self):
        self.boxes = []
        self.widgets = {}
        for group in self.modify_groups.keys():        
            self.widgets[group] = {} 

            self.widgets[group]['csv'] = ipy.Text(
                value = self.modify_groups[group].csv_path,
                description=f'CSV file for boundary condition',
                disabled=False,
                style = self.item_style,
                layout = self.item_layout
                )
        
            if len(self.widgets[group]['csv'].value) <= 0:
                for key, value in self.modify_groups[group].parameters.items():
                    if key == 'csv':
                        continue
                    unit = self.modify_groups[group].get_unit(key)
                    self.widgets[group][key] = ipy.FloatText(
                        value =  value,
                        description= f"{key} ({unit})",
                        disabled=False
                        )
                    self.widgets[group][key].layout = self.item_layout
                    self.widgets[group][key].style = self.item_style
            else:
                if self.modify_groups[group].csv_read_error != None:
                    self.widgets[group]['error'] = ipy.HTML(value=f'<b style="color:red;font-size:16px;">{self.modify_groups[group].csv_read_error}</b>')
            
            header_widget = ipy.HTML(f'<h3 style="margin:0px;padding:0px;text-align:center;">Settings for {group}</3>')
            header_widget.layout = self.item_layout
            header_widget.style = self.item_style

            self.widgets_to_display = [header_widget] + [widget for widget in self.widgets[group].values()]
            
            self.boxes.append(ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style))
        display_button = True
        if len(self.boxes) == 0:
            self.boxes.append(ipy.HTML("No boundary conditions found that can be modified."))
            display_button = False
        box = ipy.VBox(children=self.boxes, layout=self.box_layout, style = self.box_style)
        display(box)
        
        button = ipy.Button(description="Update settings", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        
        if display_button:
            display(button, output)

        button.on_click(self.update_settings_widget)
        

    def update_settings_widget(self, _):
        clear_output(wait=True)
        for group_name in self.widgets:
            for parameter_name in self.widgets[group_name].keys():
                group = self.modify_groups[group_name]
                value = self.widgets[group_name][parameter_name].value
                group.parameters[parameter_name] = value
            times, values = group.update_values()
            if times != None and values != None:
                self.plot_timeseries(times, values, group.title, group.quantity, group.get_unit(group.quantity), group.title)
        self.bc.write(self.bnd_file)

        self.display_widgets()



class BoundaryModificationGroup(WidgetStyling):
    def __init__(self, title:str, parameter:str, parent_class:ModifyBoundaries):
        self.parent = parent_class
        self.blocks_to_modify = []
        self.title = title
        self.quantity = parameter
        self.parameter_to_unit = {
            "discharge": "m3/s",
            "waterlevel": "m + NAP",
            "range": "m",
            "offset": "hours",
            "duration": "hours"
        }

        if self.quantity not in ['waterlevelbnd', 'dischargebnd']:
            raise Exception(f"Boundary paramter found ({self.parameter}) not in supported parameters")
        
        self.csv_path = None
        self.csv_read_error = None

        self.set_inputs()

    def add_boundary(self, boundary:BoundaryConditionBlock):
        self.blocks_to_modify.append(boundary)

    def set_inputs(self):
        if self.quantity == 'dischargebnd':
            self.parameters = {
                "csv": "",
                "basic_discharge": 1000,
                "peak_discharge": 2000,
                "peak_offset": 0,
                "peak_duration": 12,
                }
        elif self.quantity == 'waterlevelbnd':
            self.parameters = {
                "csv": "",
                "basic_waterlevel": 0, 
                "tidal_offset": 0,
                "tidal_range": 0,
                "peak_waterlevel_addition": 0,
                "peak_offset": 0,
                "peak_duration": 12
                }
            
    def get_unit(self, parameter:str) -> str:
        for p in self.parameter_to_unit:
            if p in parameter:
                return self.parameter_to_unit[p]
        return ''

    def update_values(self) -> Tuple[list, list]:
        times, values = None, None
        parameters = self.parameters
        if 'csv' in self.parameters.keys():
            self.csv_path = parameters.pop('csv')
        if len(self.csv_path) > 0:
            times, values = self.read_csv()
        elif self.quantity == 'dischargebnd':
            times, values = self.parent.generate_discharge_timeseries(**parameters)
        elif self.quantity == 'waterlevelbnd':
            times, values = self.parent.generate_waterlevel_timeseries(**parameters)
        if times != None and values != None:
            for block in self.blocks_to_modify:
                block.timeseries.set_timeseries(times, values)
        return times, values

    def read_csv(self) -> Tuple[list, list]:
        self.csv_read_error = None
        if os.path.exists(self.csv_path) == False:
            self.csv_read_error = f"The inputted path to csv file ({self.csv_path}) can not be found"
            return None, None
        
        if self.csv_path.endswith('.csv') == False:
            self.csv_read_error = f"The path to the csv file should end with .csv"
            return None, None
        
        try:
            csv_data = pd.read_csv(self.csv_path, delimiter=';', dtype=float)
        except Exception as e:
            self.csv_read_error = e
            return None, None
        
        if any([x not in csv_data.columns for x in ['time', 'value']]):
            self.csv_read_error = f"csv file should have a 'time' and a 'value' column. Columns found are: {list(csv_data.columns)}"
            return None, None
        
        times = csv_data['time']
        values = csv_data['value']

        if min(times) > 0:
            self.csv_read_error = f"Earliest timestep found in csv is hour={min(times)}, but times should start at hour=0"
            return None, None
        
        times_in_seconds = times * 3600
        if max(times_in_seconds) < self.parent.t_stop:
            self.csv_read_error = f"Latest timestep found in csv is hour={max(times)}, but models ends at hour={self.parent.t_stop/3600}"
            return None, None
        return list(times), list(values)          
