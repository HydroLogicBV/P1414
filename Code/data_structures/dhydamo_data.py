import geopandas as gpd
from clip_tools.clip_tools import _clip_structures_by_branches

from data_structures.dhydamo_data_model import DHydamoDataModel
from data_structures.hydamo_helpers import convert_to_dhydamo_data
from data_structures.to_dhydro_helpers import to_dhydro


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

    def from_raw_data(self, defaults: str, config: str) -> None:
        """
        Class method to load raw_data into a DHydamoDataModel. This datamodel validates data against expected values

        Args:
            defaults (str): default settings to use (should be in ./dataset_configs)
            config (str): configuration file to use (should be in ./dataset_configs)

        Returns:
            None
        """
        # load features and add to DHydamoDataModel
        self.ddm = convert_to_dhydamo_data(defaults=defaults, config=config)

        # add succesfully loaded features to features list
        self.features = []
        for key, value in self.ddm.__dict__.items():
            if value is not None:
                self.features.append(key)

    def from_dhydamo_gpkg(self, gpkg_path: str) -> None:
        """
        Class method to load data from geopackage and validate data against DHydamoDataModel

        Args:
            gpkg_path (str): file location of the geopackage to load

        Returns
            None
        """

        self.ddm = DHydamoDataModel()
        self.gpkg_path = gpkg_path
        self.features = []
        attributes = self.ddm.__dict__.keys()

        for attribute in attributes:
            try:
                # print("succesfully loaded {}".format(attribute))
                data = gpd.read_file(self.gpkg_path, layer=attribute)
            except ValueError:
                # print("failed to load {}".format(attribute))
                continue

            setattr(self.ddm, attribute, data)
            self.features.append(attribute)

    def to_dhydro(self, config: str, output_folder: str, extent: gpd.GeoDataFrame = None):
        """
        Class method that save a DHydamoDataModel to a D-HYDRO Model

        Args:
            config (str): configuration file to use (should be in ./dataset_configs)
            output_folder (str): folder to save D-HYDRO model to
            extent (gpd.GeoDataFrame): optional extent in a one-row geodataframe that is used to create a 2D model domain

        Returns:
            None
        """
        return to_dhydro(self=self, config=config, output_folder=output_folder, extent=extent)

    def to_dhydamo_gpkg(self, output_gpkg: str) -> None:
        """
        Class method that saves DHydamoDataModel to geopackage

        Args:
            output_gpkg (str): file location to save geopackage to

        Returns:
            None
        """
        self.gpkg_path = output_gpkg
        self.ddm.to_gpkg(output_gpkg=output_gpkg)
