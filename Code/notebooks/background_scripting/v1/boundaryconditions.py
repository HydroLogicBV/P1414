import os
import shutil
from shapely.geometry import Point

class TimeSeries():
    def __init__(self, lines:list=None, times:list=None, values:list=None):
        self.times = []
        self.values = []

        if lines is not None:
            self.read(lines)
        else:
            self.times = times
            self.values = values
            
    @staticmethod
    def check_if_number(s:str) -> bool:
        return all(c.isdigit() or c in '-.' for c in s)

    def read(self, lines:list):
        for line in lines:
            parts = line.strip().split()
            
            if len(parts) != 2:
                continue

            if not all(self.check_if_number(p) for p in parts):
                continue
        
            self.times.append(float(parts[0]))
            self.values.append(float(parts[1]))
        
    def set_timeseries(self, times:list, values:list):
        if len(times) != len(values):
            raise Exception(f"Times (len: {len(times)}) and values (len: {len(values)}) are not of the same length")
        self.times = times
        self.values = values            

class BoundaryConditionBlock():
    def __init__(self, lines:list):
        self.type = None
        self.property_keys = []
        self.property_values = []
        self.timeseries = None
        self.tail = None
        self.location = None
        self.boundary_group = None

        self.valid = True

        if len(lines) <= 2:
            self.valid = False
            return None
        
        self.read(lines)

        if self.type.lower() == "[forcing]":
            self.location = self.get_location()
    
    def get_location(self):
        id = self.get_property('name')
        parts = id.split('_')
        if len(parts) != 2:
            return None
        try:
            x = float(parts[0])
            y = float(parts[1])
            location = Point(x, y)
            return location
        except:
            return None

    def get_property(self, key:str, get_match_nr:int=1):
        # watch out, it always returns the first match of the property if get_last is not specified.
        get_match_nr = get_match_nr - 1
        matches = []
        for index, k in enumerate(self.property_keys):
            if k.lower() == key.lower():
                matches.append(self.property_values[index])
        if get_match_nr + 1 > len(matches):
            return None
        return matches[get_match_nr]
    
    def overwite_properties(self, key_value_pairs:list):
        self.property_keys = []
        self.property_values = []
        for key, value in key_value_pairs:
            self.property_keys.append(key)
            self.property_values.append(value)
    
    def transform_from_constant_to_timeseries(self):
        if not self.get_property('function') == 'constant':
            raise Exception(f"Trying to transform boundary {self.get_property('name')} from constant to timeseries, but boundary function is not set to constant, but to {self.get_property('function')}")
        new_key_value_pairs = [
            ['name', self.get_property('name')],
            ['function', 'timeseries'],
            ['timeInterpolation', self.get_property('timeInterpolation')],
            ['quantity', 'time'],
            ['unit', 'hours since 2000-01-01 00:00:00'],
            ['quantity', self.get_property('quantity')],
            ['unit', self.get_property('unit')],
        ]
        self.overwite_properties(new_key_value_pairs)
        self.timeseries = TimeSeries(times = [0, 9999], values = [0, 0])

    def read(self, lines:list):
        self.type = lines[0].strip()

        last_property_index = None
        for index, line in enumerate(lines):
            if "=" in line:
                key = line.split('=', maxsplit=1)[0].strip()
                value = line.split('=', maxsplit=1)[-1].strip()
                self.property_keys.append(key)
                self.property_values.append(value)
                last_property_index = index
        
        if last_property_index < len(lines):
            self.tail = lines[last_property_index+1:]
        
        if self.get_property('function') == 'timeseries':
            self.timeseries = TimeSeries(lines = self.tail)
    
    def __str__(self):
        return f"type = {self.type}, name = {self.get_property('name')}"

class BoundaryConditionsFile():
    def __init__(self, file:str, additional_block_headers:list=[]):
        self.block_headers = ["[general]", "[forcing]"]
        self.block_headers += additional_block_headers

        self.blocks = []

        self.read_file(file)
    
    def add_block(self, block_lines:list):
        block = BoundaryConditionBlock(block_lines)
        if block.valid == True:
            self.blocks.append(block)

    def read_file(self, file):
        with open(file, 'r') as f:
            lines = []
            for line in f.readlines():
                if "#" in line:
                    lines.append(line.split("#")[0])
                else:
                    lines.append(line)
        
        block_lines = []
        for line in lines:
            if line.lower().strip() in self.block_headers:
                self.add_block(block_lines)
                block_lines = []
            block_lines.append(line)
        self.add_block(block_lines)
    
    def write(self, path:str):
        max_key_length = 9
        with open(path, 'w') as f:
            for block in self.blocks:
                f.write(f"{block.type}\n")

                for i in range(len(block.property_keys)):
                    key = block.property_keys[i]
                    value = block.property_values[i]
                    spaces = ' ' * (max_key_length - len(key))
                    f.write(f"{key}{spaces}= {value}\n")
                
                if block.timeseries is not None:
                    for t in range(len(block.timeseries.times)):
                        f.write(f"{round(block.timeseries.times[t],3)}   {round(block.timeseries.values[t], 3)}\n")
                    f.write('\n')
                elif block.tail is not None:
                    for line in block.tail:
                        f.write(line)
                else:
                    f.write('\n')


