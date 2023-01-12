import geopandas as gpd
from clip_tools.clip_tools import _clip_structures_by_branches

from data_structures.dhydamo_data_model import DHydamoDataModel
from data_structures.hydamo_helpers import convert_to_dhydamo_data
from data_structures.to_dhydro_helpers import to_dhydro


class DHydamoData:
    def __init__(self):
        pass

    def clip_structures_by_branches(self, buffer: float = 1, min_overlap: float = 0.5):
        self.ddm = _clip_structures_by_branches(self, buffer=buffer, min_overlap=min_overlap)

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
        Class method to load data from geopackage and validate data against DHydamoDataModel
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

    def to_dhydro(self, config: str, output_folder: str, extent: gpd.GeoDataFrame = None):
        return to_dhydro(self=self, config=config, output_folder=output_folder, extent=extent)

    def to_dhydamo_gpkg(self, output_gpkg: str) -> None:
        self.gpkg_path = output_gpkg
        self.ddm.to_gpkg(output_gpkg=output_gpkg)
