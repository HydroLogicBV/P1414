import importlib

import geopandas as gpd
from hydrolib.dhydamo.core.hydamo import HyDAMO

from data_structures.dhydamo_data_model import DHydamoDataModel
from data_structures.dhydamo_helpers import (add_branches, add_bridges,
                                             add_culverts, add_pumps,
                                             add_weirs, build_1D_model,
                                             build_2D_model, simple_fm_model,
                                             write_model)
from data_structures.hydamo_helpers import convert_to_dhydamo_data


class DHydamoData:
    def __init__(self):
        pass

    def from_raw_data(self, defaults: str, config: str) -> None:
        # load features and add to DHydamoDataModel
        self.ddm = convert_to_dhydamo_data(defaults=defaults, config=config)

        # add succesfully loaded features to features list
        self.features = []
        for key, value in self.ddm.__dict__.items():
            if value is not None:
                self.features.append(key)

    def from_dhydamo_gpkg(self, gpkg_path):
        """
        Class method to load data from geopackage and check it against DHydamoDataModel
        """

        self.ddm = DHydamoDataModel()
        self.gpkg_path = gpkg_path
        self.features = []
        attributes = self.ddm.__dict__.keys()

        for attribute in attributes:
            try:
                data = gpd.read_file(self.gpkg_path, layer=attribute)
            except ValueError:
                continue

            setattr(self.ddm, attribute, data)
            self.features.append(attribute)

    def to_dhydro(self, config: str, output_folder: str):
        # check if data has been loaded and correct attributes are set
        if (
            (not hasattr(self, "ddm"))
            | (not hasattr(self, "gpkg_path"))
            | (not hasattr(self, "features"))
        ):
            raise AttributeError("Modeldatabase not loaded")

        # load configuration file
        model_config = getattr(importlib.import_module("dataset_configs." + config), "Models")

        # build FM model if set in config file
        if hasattr(model_config, "FM"):
            # initialize a simple FM model
            self.fm = simple_fm_model(
                start_time=model_config.FM.start_time, stop_time=model_config.FM.stop_time
            )
            self.hydamo = HyDAMO()

            # build 1D model
            if model_config.FM.one_d:
                # check if there are branches in the data
                if "waterloop" not in self.features:
                    raise ValueError("Missing branches")

                # initialize hydamo object and add branches

                self.hydamo = add_branches(
                    features=self.features, gpkg_path=self.gpkg_path, hydamo=self.hydamo
                )
                # Loop over structures and add when present in the data
                struct_functions = {
                    "brug": add_bridges,
                    "duiker": add_culverts,
                    "gemaal": add_pumps,
                    "stuw": add_weirs,
                }

                for structure, function in struct_functions.items():
                    if structure in self.features:
                        print("\nworking on {}\n".format(structure))
                        self.hydamo = function(
                            gpkg_path=self.gpkg_path,
                            hydamo=self.hydamo,
                            max_snap_dist=model_config.FM.max_snap_dist,
                        )

                # add 1D model
                self.fm, self.hydamo = build_1D_model(
                    fm=self.fm, features=self.features, hydamo=self.hydamo
                )

            # build 2D model
            if model_config.FM.two_d:
                # compute extent of
                extent = self.ddm.waterloop.dissolve(by=None).convex_hull.buffer(
                    model_config.FM.two_d_buffer
                )

                # add 2D model
                self.fm = build_2D_model(
                    dx=model_config.FM.dx,
                    dy=model_config.FM.dy,
                    elevation_raster_path=model_config.FM.elevation_raster_path,
                    extent=extent,
                    fm=self.fm,
                    one_d=model_config.FM.one_d,
                )

        # save D-HYDRO model
        write_model(self.fm, self.hydamo, output_folder=output_folder, one_d=model_config.FM.one_d)

    def to_dhydamo_gpkg(self, output_gpkg: str) -> None:
        self.gpkg_path = output_gpkg
        self.ddm.to_gpkg(output_gpkg=output_gpkg)
