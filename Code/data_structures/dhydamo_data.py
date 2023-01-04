from typing import Optional

import geopandas as gpd
import pandera as pa
from pandera.typing import DataFrame, Series
from pandera.typing.geopandas import GeoSeries
from pydantic import BaseModel

from data_structures.dhydamo_helpers import build_dhydro

## TODO: gemeten dwarsprofielen


class BasicSchema(pa.SchemaModel):
    class Config:
        name = "BasicSchema"
        strict = True  # no additional columns allowed


class PDBasicShema(BasicSchema):
    globalid: Series[str]
    geometry: GeoSeries = pa.Field(nullable=True)  # by design, to add to gpkg


class GPDBasicShema(BasicSchema):
    code: Series[str]
    globalid: Series[str]
    geometry: GeoSeries


class BrugSchema(GPDBasicShema):
    intreeverlies: Series[float] = pa.Field(nullable=True)
    lengte: Series[float] = pa.Field(nullable=True)
    ruwheid: Series[float]
    typeruwheid: Series[int]
    uittreeverlies: Series[float] = pa.Field(nullable=True)


class DuikerSchema(GPDBasicShema):
    breedteopening: Series[float] = pa.Field(nullable=True)
    hoogtebinnenonderkantbene: Series[float] = pa.Field(nullable=True)
    hoogtebinnenonderkantbov: Series[float] = pa.Field(nullable=True)
    hoogteopening: Series[float] = pa.Field(nullable=True)
    intreeverlies: Series[float] = pa.Field(nullable=True)
    lengte: Series[float] = pa.Field(nullable=True)
    ruwheid: Series[float]
    typeruwheid: Series[int]
    uittreeverlies: Series[float] = pa.Field(nullable=True)
    vormkoker: Series[float] = pa.Field(nullable=True)


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
    ruwheidhoog: Series[float]
    ruwheidlaag: Series[float]
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
    ruwheid: Series[float]  # addition to confluence
    typeruwheid: Series[int]  # addition to confluence


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


class DHydamoData:
    def __init__(self):
        pass

    # TODO: convert raw data to HYDAMO
    # TODO: save HYDAMO to GPKG

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

    def to_dhydro(self, output_folder, **kwargs):
        if (not hasattr(self, "ddm")) | (not hasattr(self, "gpkg_path")):
            raise AttributeError("Modeldatabase not loaded")

        build_dhydro(gpkg_file=self.gpkg_path, output_folder=output_folder, **kwargs)
