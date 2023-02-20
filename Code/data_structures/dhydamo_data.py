import geopandas as gpd
import pandas as pd
from geo_tools.clip_tools import _clip_structures_by_branches
from geo_tools.merge_networks import merge_networks

from data_structures.dhydamo_data_model import DHydamoDataModel
from data_structures.hydamo_helpers import check_and_fix_duplicate_code, convert_to_dhydamo_data
from data_structures.to_dhydro_helpers import to_dhydro, write_dimr


class DHydamoData:
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
        gpkg_path (str): file location of geopackage (depreciated)
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

        self.ddm = _clip_structures_by_branches(self, buffer=buffer, min_overlap=min_overlap)

    def from_raw_data(self, defaults: str, config: str, branch_snap_dist: float = 10) -> None:
        """
        Class method to load raw_data into a DHydamoDataModel. This datamodel validates data against expected values

        Args:
            defaults (str): default settings to use (should be in ./dataset_configs)
            config (str): configuration file to use (should be in ./dataset_configs)

        Returns:
            None
        """
        # load features and add to DHydamoDataModel
        ddm = convert_to_dhydamo_data(defaults=defaults, config=config)

        self._set_data(ddm=ddm, branch_snap_dist=branch_snap_dist)

    def from_dhydamo_gpkg(self, gpkg_path: str, branch_snap_dist: float = 10) -> None:
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
        ddm = DHydamoDataModel()

        # loop over datamodel attributes and check if they are presentin the geopackage
        # if so, set them to the datamodel
        attributes = ddm.__dict__.keys()
        for attribute in attributes:
            try:
                # print("succesfully loaded {}".format(attribute))
                data = gpd.read_file(gpkg_path, layer=attribute)
            except ValueError:
                # print("failed to load {}".format(attribute))
                continue

            setattr(ddm, attribute, data)

        # set datamodel to self
        self._set_data(ddm=ddm, branch_snap_dist=branch_snap_dist)

    def to_dhydamo_gpkg(self, output_gpkg: str) -> None:
        """
        Class method that saves DHydamoDataModel to geopackage

        Args:
            output_gpkg (str): file location to save geopackage to

        Returns:
            None
        """
        # self.gpkg_path = output_gpkg
        self.ddm.to_gpkg(output_gpkg=output_gpkg)

    def to_dhydro(self, config: str, output_folder: str, write=True):
        """
        Class method that converts a DHydamoDataModel to a D-HYDRO Model and saves unless write=False

        Args:
            config (str): configuration file to use (should be in ./dataset_configs)
            output_folder (str): folder to save D-HYDRO model to

        Returns:
            None
        """
        to_dhydro(self=self, config=config)

        if write:
            self.write_dimr(output_folder=output_folder)

    def write_dimr(self, output_folder: str):
        return write_dimr(fm=self.fm, output_folder=output_folder)

    def _set_data(self, ddm: DHydamoDataModel, branch_snap_dist: float) -> None:
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
