from typing import Optional

import geopandas as gpd
import pandera as pa
from hydrolib.dhydamo.core.hydamo import HyDAMO
from pandera.typing import DataFrame, Series
from pandera.typing.geopandas import GeoSeries
from pydantic import BaseModel
from tqdm import tqdm

from data_structures.dhydamo_helpers import (
    add_branches,
    add_bridges,
    add_culverts,
    add_pumps,
    add_weirs,
    build_1D_model,
    simple_fm_model,
    write_model,
)
from data_structures.dhydamo_helpers_old import build_dhydro

## TODO: gemeten dwarsprofielen


class BasicSchema(pa.SchemaModel):
    class Config:
        name = "BasicSchema"
        strict = True  # no additional columns allowed


class PDBasicShema(BasicSchema):
    globalid: Series[str]
    geometry: GeoSeries = pa.Field(nullable=True)  # by design, to add to gpkg


class GPDBasicShema(BasicSchema):
    code: Series[str] = pa.Field(coerce=True)
    globalid: Series[str]
    geometry: GeoSeries


class BrugSchema(GPDBasicShema):
    intreeverlies: Series[float] = pa.Field(nullable=True)
    lengte: Series[float] = pa.Field(nullable=True)
    ruwheid: Series[float]
    typeruwheid: Series[str]
    uittreeverlies: Series[float] = pa.Field(nullable=True)


class DuikerSchema(GPDBasicShema):
    breedteopening: Series[float] = pa.Field(nullable=True)
    hoogtebinnenonderkantbene: Series[float] = pa.Field(nullable=True)
    hoogtebinnenonderkantbov: Series[float] = pa.Field(nullable=True)
    hoogteopening: Series[float] = pa.Field(nullable=True)
    intreeverlies: Series[float] = pa.Field(nullable=True)
    lengte: Series[float] = pa.Field(nullable=True)
    ruwheid: Series[float]
    typeruwheid: Series[str]
    uittreeverlies: Series[float] = pa.Field(nullable=True)
    vormkoker: Series[float] = pa.Field(nullable=True)  # accepteer enkele waarden


class GemaalSchema(GPDBasicShema):
    pass


class Hydroobject_normgpSchema(PDBasicShema):
    hydroobjectid: Series[str]
    normgeparamprofielid: Series[str]


class KunstwerkopeningSchema(PDBasicShema):
    afvoercoefficient: Series[float]
    hoogstedoorstroombreedte: Series[float] = pa.Field(nullable=True)
    hoogstedoorstroomhoogte: Series[float] = pa.Field(nullable=True)
    laagstedoorstroombreedte: Series[float] = pa.Field(nullable=True)
    laagstedoorstroomhoogte: Series[float] = pa.Field(nullable=True)
    stuwid: Series[str]
    vormopening: Series[float] = pa.Field(nullable=True, coerce=True)


class NormgeparamprofielwaardeSchema(BasicSchema):
    geometry: GeoSeries = pa.Field(nullable=True)  # by design, to add to gpkg
    normgeparamprofielid: Series[str]
    ruwheidhoog: Series[float] = pa.Field(coerce=True)
    ruwheidlaag: Series[float] = pa.Field(coerce=True)
    soortparameter: Series[str]
    typeruwheid: Series[str]
    waarde: Series[float]


class PompSchema(PDBasicShema):
    code: Series[str]  # addition to confluence
    gemaalid: Series[str]
    maximalecapaciteit: Series[float] = pa.Field(nullable=True)


class RegelmiddelSchema(PDBasicShema):
    code: Series[str]  # addition to confluence
    kunstwerkopeningid: Series[str]
    overlaatonderlaat: Series[str]
    soortregelbaarheid: Series[float] = pa.Field(nullable=True, coerce=True)
    stuwid: Series[str]


class SturingSchema(PDBasicShema):
    bovengrens: Series[float] = pa.Field(nullable=True)
    code: Series[str]  # addition to confluence
    doelvariabele: Series[str]  # addition to confluence
    ondergrens: Series[float] = pa.Field(nullable=True)
    pompid: Series[str]
    streefwaarde: Series[float] = pa.Field(nullable=True)


class StuwSchema(GPDBasicShema):
    afvoercoefficient: Series[float]
    code: Series[str]
    soortstuw: Series[float] = pa.Field(nullable=True, coerce=True)  # addition to confluence


class WaterloopSchema(GPDBasicShema):
    # ruwheid: Series[float]  # addition to confluence
    typeruwheid: Series[str]  # addition to confluence


class DHydamoDataModel(BaseModel):
    brug: Optional[DataFrame[BrugSchema]]
    duiker: Optional[DataFrame[DuikerSchema]]
    gemaal: Optional[DataFrame[GemaalSchema]]
    hydroobject_normgp: Optional[DataFrame[Hydroobject_normgpSchema]]
    kunstwerkopening: Optional[DataFrame[KunstwerkopeningSchema]]
    normgeparamprofielwaarde: Optional[DataFrame[NormgeparamprofielwaardeSchema]]
    pomp: Optional[DataFrame[PompSchema]]
    regelmiddel: Optional[DataFrame[RegelmiddelSchema]]
    sturing: Optional[DataFrame[SturingSchema]]
    stuw: Optional[DataFrame[StuwSchema]]
    waterloop: Optional[DataFrame[WaterloopSchema]]

    class Config:
        validate_assignment = True  # validates new attributes that are assigned

    def to_gpkg(self, output_gpkg: str) -> None:
        # fields = list(self.__dict__.keys())
        for key, value in self.__dict__.items():
            if value is not None:
                getattr(self, key).to_file(output_gpkg, layer=key, driver="GPKG")


class DHydamoData:
    def __init__(self):
        pass

    # TODO: convert raw data to HYDAMO

    def from_gpkg(self, gpkg_path):
        """
        Class method to load data from geopackage and check it against DHydamoDataModel
        """
        self.gpkg_path = gpkg_path

        brug = gpd.read_file(self.gpkg_path, layer="brug")
        duiker = gpd.read_file(self.gpkg_path, layer="duiker")
        gemaal = gpd.read_file(self.gpkg_path, layer="gemaal")
        hydroobject_normgp = gpd.read_file(self.gpkg_path, layer="hydroobject_normgp")
        kunstwerkopening = gpd.read_file(self.gpkg_path, layer="kunstwerkopening")
        normgeparamprofielwaarde = gpd.read_file(self.gpkg_path, layer="normgeparamprofielwaarde")
        pomp = gpd.read_file(self.gpkg_path, layer="pomp")
        regelmiddel = gpd.read_file(self.gpkg_path, layer="regelmiddel")
        sturing = gpd.read_file(self.gpkg_path, layer="sturing")
        stuw = gpd.read_file(self.gpkg_path, layer="stuw")
        waterloop = gpd.read_file(self.gpkg_path, layer="waterloop")

        self.ddm = DHydamoDataModel(
            brug=brug,
            duiker=duiker,
            gemaal=gemaal,
            hydroobject_normgp=hydroobject_normgp,
            kunstwerkopening=kunstwerkopening,
            normgeparamprofielwaarde=normgeparamprofielwaarde,
            pomp=pomp,
            regelmiddel=regelmiddel,
            sturing=sturing,
            stuw=stuw,
            waterloop=waterloop,
        )
        self.features = list(self.ddm.__dict__.keys())

    def to_dhydro(self, output_folder, **kwargs):
        if (not hasattr(self, "ddm")) | (not hasattr(self, "gpkg_path")):
            raise AttributeError("Modeldatabase not loaded")

        if "waterloop" not in self.features:
            raise ValueError("Missing branches")

        self.hydamo = HyDAMO()
        self.hydamo = add_branches(
            features=self.features, gpkg_path=self.gpkg_path, hydamo=self.hydamo
        )

        struct_functions = {
            "brug": add_bridges,
            "duiker": add_culverts,
            "gemaal": add_pumps,
            "stuw": add_weirs,
        }

        for structure, function in struct_functions.items():
            if structure in self.features:
                print("working on {}".format(structure))
                self.hydamo = function(
                    gpkg_path=self.gpkg_path, hydamo=self.hydamo, max_snap_dist=1
                )

        self.fm = simple_fm_model()
        self.fm, self.hydamo = build_1D_model(
            fm=self.fm, features=self.features, hydamo=self.hydamo
        )
        write_model(self.fm, self.hydamo, output_folder=output_folder)

    def to_dhydro_old(self, output_folder, **kwargs):
        if (not hasattr(self, "ddm")) | (not hasattr(self, "gpkg_path")):
            raise AttributeError("Modeldatabase not loaded")

        build_dhydro(
            features=self.features, gpkg_file=self.gpkg_path, output_folder=output_folder, **kwargs
        )
