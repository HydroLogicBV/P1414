from typing import Optional

import pandera as pa
from pandera.typing import DataFrame, Series
from pandera.typing.geopandas import GeoSeries
from pydantic import BaseModel

roughness_mapping = {
    "Chezy": "Chezy",
    "Manning": "Manning",
    "StricklerKn": "StricklerNikuradse",
    "StricklerKs": "Strickler",
    "White Colebrook": "WhiteColebrook",
    "Bos en Bijkerk": "deBosBijkerk",
    "Onbekend": "Strickler",
    "Overig": "Strickler",
}
ROUGHNESS_LIST = list(roughness_mapping)
HYDAMO_SHAPE_NUMS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 98, 99]
HYDAMO_WEIR_TYPES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24, 25, 26, 98, 99]
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
    intreeverlies: Series[float] = pa.Field(ge=0, le=1)
    lengte: Series[float] = pa.Field(gt=0)
    ruwheid: Series[float] = pa.Field(gt=0)
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_LIST)
    uittreeverlies: Series[float] = pa.Field(ge=0, le=1)


class DuikerSchema(GPDBasicShema):
    breedteopening: Series[float] = pa.Field(gt=0)
    doorstroomopening: Series[str]
    hoogtebinnenonderkantbene: Series[float]
    hoogtebinnenonderkantbov: Series[float]
    hoogteopening: Series[float] = pa.Field(gt=0)
    intreeverlies: Series[float] = pa.Field(ge=0, le=1)
    lengte: Series[float] = pa.Field(gt=0)
    ruwheid: Series[float]
    typeruwheid: Series[str]
    uittreeverlies: Series[float] = pa.Field(ge=0, le=1)
    vormkoker: Series[int] = pa.Field(
        isin=HYDAMO_SHAPE_NUMS, coerce=True
    )  # accepteer enkel waarden volgend hydamo standaard


class GemaalSchema(GPDBasicShema):
    pass


class Hydroobject_normgpSchema(PDBasicShema):
    hydroobjectid: Series[str]
    normgeparamprofielid: Series[str]


class KunstwerkopeningSchema(PDBasicShema):
    afvoercoefficient: Series[float]
    hoogstedoorstroombreedte: Series[float] = pa.Field(gt=0)
    hoogstedoorstroomhoogte: Series[float]
    laagstedoorstroombreedte: Series[float] = pa.Field(gt=0)
    laagstedoorstroomhoogte: Series[float]
    stuwid: Series[str]
    vormopening: Series[int] = pa.Field(
        isin=HYDAMO_SHAPE_NUMS, coerce=True
    )  # accepteer enkel waarden volgend hydamo standaard


class NormgeparamprofielwaardeSchema(BasicSchema):
    geometry: GeoSeries = pa.Field(nullable=True)  # by design, to add to gpkg
    normgeparamprofielid: Series[str]
    ruwheidhoog: Series[float] = pa.Field(gt=0, coerce=True)
    ruwheidlaag: Series[float] = pa.Field(gt=0, coerce=True)
    soortparameter: Series[str]
    typeruwheid: Series[str]
    waarde: Series[float]


class PompSchema(PDBasicShema):
    code: Series[str] = pa.Field(coerce=True)  # addition to confluence
    gemaalid: Series[str]
    maximalecapaciteit: Series[float] = pa.Field(ge=0)


class RegelmiddelSchema(PDBasicShema):
    code: Series[str]  # addition to confluence
    kunstwerkopeningid: Series[str]
    overlaatonderlaat: Series[str]
    soortregelbaarheid: Series[int] = pa.Field(isin=[1, 2, 3, 4, 98, 99], coerce=True)
    stuwid: Series[str]


class SturingSchema(PDBasicShema):
    bovengrens: Series[float] = pa.Field(nullable=True)
    code: Series[str] = pa.Field(coerce=True)  # addition to confluence
    doelvariabele: Series[str]  # addition to confluence
    ondergrens: Series[float] = pa.Field(nullable=True)
    pompid: Series[str]
    streefwaarde: Series[float] = pa.Field(nullable=True)


class StuwSchema(GPDBasicShema):
    afvoercoefficient: Series[float] = pa.Field(coerce=True)
    code: Series[str]
    soortstuw: Series[float] = pa.Field(
        isin=HYDAMO_WEIR_TYPES, coerce=True
    )  # addition to confluence


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
        """save fields that are not None to gpkg"""
        for key, value in self.__dict__.items():
            if value is not None:
                getattr(self, key).to_file(output_gpkg, layer=key, driver="GPKG")