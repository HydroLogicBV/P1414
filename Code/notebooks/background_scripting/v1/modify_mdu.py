import os
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import json
from datetime import datetime

class ModifyMDU(WidgetStyling):
    def __init__(self, model_folder):
        files = os.listdir(os.path.join(model_folder, 'dflowfm'))
        mdu_file = [file for file in files if file.endswith('.mdu')][0]
        self.mdu_path = os.path.join(model_folder, 'dflowfm', mdu_file)
        self.run_bat_file = os.path.join(model_folder, 'run.bat')
        self.mdu_lines = self.read_mdu()
        
        self.modify_parameter_mdu(parameter = "statsInterval", new_value = 1) # always do this

        self.settings_to_modify = ['tStart', 'tStop', 'mapInterval']
        self.settings = {}
        self.settings['refDate'] = datetime.strptime(self.read_parameter_mdu('refDate'), '%Y%m%d')
        
        self.set_default_layout_and_styling()


        self.widgets = {}
        for setting in self.settings_to_modify:
            self.settings[setting] =  float(self.read_parameter_mdu(setting))
            self.widgets[setting] = ipy.FloatText(
                value = self.settings[setting],
                description=f'{setting}:',
                disabled=False
                )
            
            self.widgets[setting].layout = self.item_layout
            self.widgets[setting].style = self.item_style
        
        self.settings['DHYDRO location'] = self.read_run_bat()
        self.widgets['DHYDRO location'] = ipy.Text(
            value = self.settings['DHYDRO location'],
            description=f'Location of dhydro run_dimr:',
            disabled=False,
            style = self.item_style,
            layout = self.item_layout
            )

        self.widgets_to_display = [self.widgets[setting] for setting in self.widgets.keys()]
        
    def read_run_bat(self):
        with open(self.run_bat_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("call "):
                path =  line.split('"')[1]
        return path

    def modify_run_bat(self):
        with open(self.run_bat_file, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("call "):
                lines[i] = 'call \"{}\"\n'.format(self.settings['DHYDRO location'])
        with open(self.run_bat_file, 'w') as f:
            for line in lines:
                f.write(line)


    def search_string_in_file(self, string, lines):
        line_numbers, strings = [], []
        for i, line in enumerate(lines):
            if string in line:
                line_numbers.append(i)
                strings.append(line)
        if len(strings) == 1:
            return line_numbers[0], strings[0]
        else:
            return None
        
    def read_mdu(self):
        with open(os.path.join(self.mdu_path), 'r') as f:
            lines = f.readlines()
        return lines
    
    def save_mdu(self):
        with open(os.path.join(self.mdu_path), 'w') as f:
            for line in self.mdu_lines:
                f.write(line)
    
    def modify_parameter_mdu(self, parameter, new_value):
        search_string = f" {parameter} " # add spaces to prevent accidentaly having wrong value
        index, line = self.search_string_in_file(search_string, self.mdu_lines)
        key = line.split("=")[0]
        value = f"= {new_value}"
        self.mdu_lines[index] = f"{key}{value}\n"
    
    def read_parameter_mdu(self, parameter):
        search_string = f" {parameter} " # add spaces to prevent accidentaly having wrong value
        index, line = self.search_string_in_file(search_string, self.mdu_lines)
        value = line.split("=")[-1].strip(' ')
        if '\n' in value:
            value = value.replace('\n', '')
        return value

    def display_widgets(self):
        items = ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style)
        display(items)
        
        button = ipy.Button(description="Update settings", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.update_settings_widget)
    
    def update_settings_widget(self, b):
        for setting in self.settings_to_modify:
            widget_value = self.widgets[setting].value
            self.settings[setting] = widget_value
            self.modify_parameter_mdu(parameter = setting, new_value = widget_value)

        self.settings['DHYDRO location'] = self.widgets['DHYDRO location'].value
        self.modify_run_bat()    
        
        self.save_mdu()
        clear_output(wait=True)
        self.display_widgets()
        display("MDU settings are:")
        print_settings_dict = {k: self.settings[k] for k in self.widgets.keys()}
        print(json.dumps(print_settings_dict, indent=4))
