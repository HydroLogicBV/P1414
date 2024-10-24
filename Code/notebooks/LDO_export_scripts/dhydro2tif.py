import sys
import os
import numpy as np
from flowmeshreader import load_meta_data, load_map_data, mesh_to_tiff, load_mesh, load_data
import os
from typing import Tuple

class DhydroToTif():
    def __init__(self, model_dir:str, file_locations:dict, resolution:float=None, silent:bool=False):
        self.silent = silent
        self.model_dir = model_dir
        self.file_locations = file_locations

        self.fm_dir = os.path.join(self.model_dir, self.file_locations['fm_dir'])
        self.output_dir =  os.path.join(self.fm_dir, self.file_locations['output_dir'])
        self.ldo_export_dir = os.path.join(self.output_dir, 'LDO_export')
        self.create_export_dir(self.ldo_export_dir)
        self.map_file = os.path.join(self.output_dir, self.file_locations['map_file'])

        if resolution == None:
            self.resolution = self.get_resolution(self.map_file)
            if not self.silent:
                print(f"Using a resolution of {self.resolution} meters, interpreted by reading the output")
        else:
            if not self.silent:
                print(f"Using a resolution of {self.resolution} meters, as specified by the user")

        self.interpolation = 'nearest'
        self.distance_tolerance = self.resolution * 1.5

        self.times, self.timestep = self.get_time()
        self.show_times()
        
    def get_resolution(self, map_file:str) -> float:
        """Get resolution of a rectangular mesh by reading the map file.
        Input:
            map_file (str): location of the map netcdf file.

        Returns:
            float: resolution of 2D grid
        """
        self.edges = load_map_data(map_file, 'Mesh2d_edge_nodes')
        self.node_x = load_data(map_file, 'Mesh2d_node_x')
        self.node_y = load_data(map_file, 'Mesh2d_node_y')
        self.node_id_offset = np.min(self.edges)
        distances = []
        for edge in self.edges[:50]:
            from_node, to_node = [int(e - self.node_id_offset) for e in edge]
            edge_length = ((self.node_x[from_node] - self.node_x[to_node]) ** 2 + (self.node_y[from_node] - self.node_y[to_node]) ** 2) ** 0.5
            distances.append(edge_length)
        unique_distances = set(distances)
        if len(unique_distances) != 1:
            raise Exception(f"Error occured when interpreting model resolution. Found these distances: {unique_distances}. Note: this script is meant to be run with rectangular meshes. You can also manually specify the resolution in the intialisation of {self.__name__}")
        return distances[0]

    def show_times(self):
        times = load_data(self.map_file, 'time')
        if not self.silent:
            print(
                    (f"netcdf start time at t={min(times)} seconds\n"
                    f"netcdf end time at t={max(times)} seconds\n"
                    f"Timestep of output data is {times[1]-times[0]} seconds")
            )

    def get_time(self):
        times = load_data(self.map_file, 'time')
        return times, times[1]-times[0]

    def get_start_index(self, t_breach:int) -> Tuple[int, float]:
        start_index = max(np.where(self.times > t_breach)[0][0] - 1, 0)
        return start_index, self.times[start_index]
    
    def create_export_dir(self, export_dir:str):
        if not os.path.exists(export_dir):
            if os.path.exists(os.path.dirname(export_dir)):
                os.mkdir(export_dir)
            else:
                raise Exception(f"This directory does not exist: {os.path.dirname(export_dir)}")
    
    def create_waterdepth_animation(self, t_breach:int, netcdf_variable:str="Mesh2d_waterdepth", output_dir:str=None):
        map_data = load_map_data(self.map_file, netcdf_variable)
        start_index, start_time= self.get_start_index(t_breach)
        if output_dir == None:
            output_dir = os.path.join(self.ldo_export_dir, 'waterdepth_animation')
        self.create_export_dir(output_dir)
        if not self.silent:
            print("=========================================================")
            print(f"Generating waterdepth animation tif files, saving them to {output_dir}")
        for t in range(start_index, map_data.shape[0]):
            if not self.silent:
                print(f" ({t+1}/{map_data.shape[0]}) Generated tif file for timestep {t} (t={int(self.times[t])}, t since breach={int(self.times[t] - start_time)})")
            _, _, grid_data = self.mesh_to_tiff_wrapper(
                data=map_data[t, :],
                output_file_path=os.path.join(output_dir, f"{int(self.times[t] - start_time)}.tiff")
                )
        if not self.silent:
            print(f"Finished generating tiff files\n")
    
    def create_maximum_waterdepth_tiff(self, netcdf_variable:str='Mesh2d_waterdepth', output_file:str=None):
        map_data = load_map_data(self.map_file, netcdf_variable)
        map_data = np.max(map_data, axis=0)
        if output_file == None:
            output_file = os.path.join(self.ldo_export_dir, 'depth-maximum.tiff')
        if not self.silent:
            print("=========================================================")
            print(f"Generating maximum waterdepth tiff file, saving it to {output_file}")
        _, _, grid_data = self.mesh_to_tiff_wrapper(data=map_data,output_file_path=output_file)
        if not self.silent:
            print(f"Finished generating tiff file\n")
    
    def create_maximum_velocity_tiff(self, netcdf_variable:str='Mesh2d_ucmag', output_file:str=None):
        map_data = load_map_data(self.map_file, netcdf_variable)
        map_data = np.max(map_data, axis=0)
        if output_file is None:
            output_file = os.path.join(self.ldo_export_dir, 'velocity-maximum.tiff')
        if not self.silent:
            print("=========================================================")
            print(f"Generating maximum velocity tiff file, saving it to {output_file}")
        _, _, grid_data = self.mesh_to_tiff_wrapper(data=map_data,output_file_path=output_file)
        if not self.silent:
            print(f"Finished generating tiff file\n")
    
    def create_bathymetry_tiff(self, netcdf_variable:str='Mesh2d_flowelem_bl', output_file:str=None):
        map_data = load_data(self.map_file, netcdf_variable)
        if output_file is None:
            output_file = os.path.join(self.ldo_export_dir, 'bathymetry.tiff')
        if not self.silent:
            print("=========================================================")
            print(f"Generating bathymetry tiff file, saving it to {output_file}")
        _, _, grid_data = self.mesh_to_tiff_wrapper(data=map_data,output_file_path=output_file)
        if not self.silent:
            print(f"Finished generating tiff file\n")

    def mesh_to_tiff_wrapper(self, data:np.ndarray, output_file_path:str):
        return mesh_to_tiff(
            data=data,
            input_file_path=self.map_file,
            output_file_path=output_file_path,
            resolution=self.resolution,
            distance_tol=self.distance_tolerance,
            interpolation=self.interpolation
        )

    



