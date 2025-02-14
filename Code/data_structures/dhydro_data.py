import importlib

import geopandas as gpd
import pandas as pd
import os
import uuid
from geo_tools.clip_tools import _clip_structures_by_branches
from geo_tools.merge_networks import merge_networks

from data_structures.fixedweirs_helpers import (create_dambreak_data,
                                                create_fixed_weir_data)
from data_structures.hydamo_helpers import (check_and_fix_duplicate_code,
                                            convert_to_dhydamo_data)
from data_structures.roi_data_model import ROIDataModel as DataModel
from data_structures.to_dhydro_helpers import to_dhydro, write_dimr
import os


class DHydroData:
    """
    Helper class to convert raw data to:
    1. load raw data into a DHydamoDataModel
    2. load a DHydamo geopackage into a DHydamoDataModel
    3. save a DHydamoDataModel into a DHydamo geopackage
    4. save a DHYdamoDataModel into a D-HYDRO model
    5. Clip structures by (buffered) branches extent

    Attributes:
        ddm (DHydamoDataModel): instance of DHydamoDataModel that is used by this clas
        features (list): list of features that are present in self.ddm and do not have None values
    """

    def __init__(self):
        pass

    def clip_structures_by_branches(self, buffer: float = 1, min_overlap: float = 0.5) -> None:
        """
        Class method to drop structures from the data that are farther from a branch than "buffer"
        Structures that are linestring are required to overlap with a branch for at least "min_overlap"

        Args:
            buffer (float): distance from branch that structures are kept (default = 1m)
            min_overlap (float): percentage of branch and LineString structure that should overlap in order to keep it.
                                 this drop, for example, culverts that are perpendicular to a branch

        Returns:
            None
        """
        print(f'Clipping structures with a buffer of {buffer} meter')
        self.ddm = _clip_structures_by_branches(self, buffer=buffer, min_overlap=min_overlap)

    def dambreaks_from_config(self, config: str, defaults: str):
        dm = DataModel()
        dm = create_dambreak_data(config=config, defaults=defaults, dm=dm)

        self._set_ddm(ddm=dm)

    def fixed_weirs_from_raw_data(self, config: str, defaults: str, min_length: float = None):
        dm = DataModel()
        dm = create_fixed_weir_data(config=config, defaults=defaults, dm=dm, min_length=min_length)
        print('Fixed weirs added')
        self._set_ddm(ddm=dm)

    def hydamo_from_raw_data(
        self, config: str, defaults: str, GIS_folder: str, dhydro_config: str, branch_snap_dist: float = 10
    ) -> None:
        """
        Class method to load raw_data into a DHydamoDataModel. This datamodel validates data against expected values

        Args:
            defaults (str): default settings to use (should be in ./dataset_configs)
            config (str): configuration file to use (should be in ./dataset_configs)

        Returns:
            None
        """
        # load features and add to DHydamoDataModel
        ddm = DataModel()
        ddm = convert_to_dhydamo_data(ddm=ddm, defaults=defaults, config=config, GIS_folder=GIS_folder, dhydro_config=dhydro_config, branch_snap_dist=branch_snap_dist)

        self._set_ddm(branch_snap_dist=branch_snap_dist, ddm=ddm)

    def hydamo_from_gpkg(self, gpkg_path: str, branch_snap_dist: float = 10) -> None:
        """
        Class method to load data from geopackage and validate data against DHydamoDataModel

        Args:
            gpkg_path (str): file location of the geopackage to load

        Returns
            None
        """

        # set geopackage path to self
        # self.gpkg_path = gpkg_path

        # initialize datamodel
        ddm = DataModel()

        # loop over datamodel attributes and check if they are presentin the geopackage
        # if so, set them to the datamodel
        attributes = ddm.__dict__.keys()
        for attribute in attributes:
            if not os.path.exists(gpkg_path):
                print(f"WARNING: gpkg path '{gpkg_path}' not found, skipping this GPKG layer")
                continue
            
            try:
                data = gpd.read_file(gpkg_path, layer=attribute)
                print("succesfully loaded {}".format(attribute))
            except ValueError:
                print("failed to load {}".format(attribute))
                continue
            
            if len(data) > 0:
                setattr(ddm, attribute, data)


        # set datamodel to self
        self._set_ddm(branch_snap_dist=branch_snap_dist, ddm=ddm)

    def hydamo_to_gpkg(self, output_gpkg: str) -> None:
        """
        Class method that saves DHydamoDataModel to geopackage

        Args:
            output_gpkg (str): file location to save geopackage to

        Returns:
            None
        """
        # self.gpkg_path = output_gpkg
        self.ddm.to_gpkg(output_gpkg=output_gpkg)

    def relocate_comments_in_roughness_file(self, textfile: str):
        # Rename the old textfile such that the new textfile can have the old name
        textfile_old = textfile[:-4] + r'_oud.ini'
        os.rename(textfile, textfile_old)        
        with open(textfile_old, 'r') as f:
            lines = f.readlines()
        
        better_lines = [sub.replace('# numLevels lines containing space separated lists', '\n# numLevels lines containing space separated lists') for sub in lines]
        
        with open(textfile, 'w') as f:
            for line in better_lines:
                f.write(line)

    def to_dhydro(self, config: str, output_folder: str, load_mesh2d_path:str = None,defaults: str = "defaults", write=True):
        """
        Class method that converts a DHydamoDataModel to a D-HYDRO Model and saves unless write=False

        Args:
            config (str): configuration file to use (should be in ./dataset_configs)
            output_folder (str): folder to save D-HYDRO model to

        Returns:
            None
        """
        # check if data has been loaded and correct attributes are set
        if (not hasattr(self, "ddm")) | (not hasattr(self, "features")):
            raise AttributeError("Modeldatabase not loaded")

        if hasattr(importlib.import_module("dataset_configs." + config), "Dambreak"):
            self.dambreaks_from_config(config=config, defaults=defaults)

        if not load_mesh2d_path == None:
            if not os.path.exists(load_mesh2d_path):
                raise AttributeError("Existing model mesh not found, check the path")

        to_dhydro(self=self, config=config, load_mesh2d_path=load_mesh2d_path, output_folder=output_folder)

        if write:
            self.write_dimr(output_folder=output_folder)

        # Relocate the comments in the roughness files if those files exist
        if write:
            roughness_files_to_check = [output_folder + r'\dflowfm\roughness_FloodPlain1.ini',
                                        output_folder + r'\dflowfm\roughness_Main.ini']
            for rfile in roughness_files_to_check:
                if os.path.exists(rfile):
                    self.relocate_comments_in_roughness_file(rfile)

    def write_dimr(self, output_folder: str):
        return write_dimr(fm=self.fm, output_folder=output_folder)

    def _set_ddm(self, ddm: DataModel, branch_snap_dist: float = None) -> None:
        """
        Class method to add a DHydamoDataModel to self while checking for pre-existing data.
        This method simply assigns the data if there is no pre-existing data.
        Otherwise, it appends to the pre-existing data

        Args:
            ddm (DHydamoDataModel): a DHydamoDataModel that is to be added to self.ddm

        Returns:
            None
        """

        # check if a datamodel already exists
        if not hasattr(self, "ddm"):
            # if not, just assign
            self.ddm = ddm

            # and fill list of features present
            self.features = []
            for key, value in self.ddm.__dict__.items():
                # check if attribute is not None
                if value is not None:
                    self.features.append(key)

        # if it does exist, check per attribute in datamodel if it already exists
        else:
            for key, value in ddm.__dict__.items():
                # check if attribute is not None
                if value is not None:
                    # if an attribute does not exist, simply assign and update features list
                    if getattr(self.ddm, key) is None:
                        setattr(self.ddm, key, getattr(ddm, key))
                        self.features.append(key)

                    # else, concatonate the old and new geodataframe and assign to the datamodel
                    else:
                        in_gdf = getattr(ddm, key)

                        if key == "waterloop":
                            in_gdf = merge_networks(
                                data_base_input=in_gdf,
                                data_match_input=getattr(self.ddm, key),
                                max_dist=branch_snap_dist,
                            )

                        new_gdf = gpd.GeoDataFrame(
                            pd.concat([getattr(self.ddm, key), in_gdf], ignore_index=True)
                        )
                        new_gdf = check_and_fix_duplicate_code(new_gdf)

                        setattr(self.ddm, key, new_gdf)
