from typing import Any, Optional

import pandera as pa
from pandera.typing import DataFrame, Series
from pandera.typing.geopandas import GeoSeries
from pydantic import BaseModel

from data_structures.dhydamo_data_model_checks import (
    geometry_check,
    globalid_check,
    none_geometry_check,
)
from data_structures.hydamo_globals import (
    HYDAMO_SHAPE_NUMS,
    HYDAMO_WEIR_TYPES,
    MANAGEMENT_DEVICE_TYPES,
    ROUGHNESS_MAPPING_LIST,
)


class BasicSchema(pa.SchemaModel):
    class Config:
        coerce = True  # allow datatype conversion
        name = "BasicSchema"
        strict = True  # no additional columns allowed


class PDBasicShema(BasicSchema):
    globalid: Series[str] = pa.Field(unique=True)
    geometry: GeoSeries = pa.Field(nullable=True)

    @pa.check("globalid", element_wise=True, name="globalidcheck")
    def _globalid_check(cls, globalid: str) -> bool:
        return globalid_check(globalid=globalid)

    @pa.check("geometry", element_wise=True, name="None geometry check")
    def _none_geometry_check(cls, geometry: None) -> bool:
        return none_geometry_check(geometry=geometry)


class GPDBasicShema(BasicSchema):
    globalid: Series[str] = pa.Field(unique=True)
    geometry: GeoSeries

    @pa.check("globalid", element_wise=True, name="globalidcheck")
    def _globalid_check(cls, globalid: str) -> bool:
        return globalid_check(globalid=globalid)

    @pa.check("geometry", element_wise=True, name="Point or line geometry check")
    def _geometry_check(cls, geometry: None) -> bool:
        return geometry_check(geometry=geometry)


class BrugSchema(GPDBasicShema):
    code: Series[str] = pa.Field(coerce=True)
    intreeverlies: Series[float] = pa.Field(ge=0, le=1)
    lengte: Series[float] = pa.Field(gt=0)
    ruwheid: Series[float] = pa.Field(gt=0)
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)
    uittreeverlies: Series[float] = pa.Field(ge=0, le=1)


class DuikerSchema(GPDBasicShema):
    breedteopening: Series[float] = pa.Field(gt=0)
    code: Series[str] = pa.Field(coerce=True)
    doorstroomopening: Series[str]
    hoogtebinnenonderkantbene: Series[float]
    hoogtebinnenonderkantbov: Series[float]
    hoogteopening: Series[float] = pa.Field(gt=0)
    intreeverlies: Series[float] = pa.Field(ge=0, le=1)
    lengte: Series[float] = pa.Field(gt=0)
    ruwheid: Series[float]
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)
    uittreeverlies: Series[float] = pa.Field(ge=0, le=1)
    vormkoker: Series[int] = pa.Field(
        isin=HYDAMO_SHAPE_NUMS, coerce=True
    )  # accepteer enkel waarden volgend hydamo standaard


class GemaalSchema(GPDBasicShema):
    code: Series[str] = pa.Field(coerce=True)


class Hydroobject_normgpSchema(PDBasicShema):
    hydroobjectid: Series[str] = pa.Field(unique=True)
    normgeparamprofielid: Series[str] = pa.Field(unique=True)


class KunstwerkopeningSchema(PDBasicShema):
    afvoercoefficient: Series[float]
    hoogstedoorstroombreedte: Series[float] = pa.Field(gt=0)
    hoogstedoorstroomhoogte: Series[float]
    laagstedoorstroombreedte: Series[float] = pa.Field(gt=0)
    laagstedoorstroomhoogte: Series[float]
    stuwid: Series[str] = pa.Field(unique=True)
    vormopening: Series[int] = pa.Field(
        isin=HYDAMO_SHAPE_NUMS, coerce=True
    )  # accepteer enkel waarden volgend hydamo standaard


class NormgeparamprofielwaardeSchema(BasicSchema):
    geometry: GeoSeries = pa.Field(nullable=True)  # by design, to add to gpkg
    normgeparamprofielid: Series[str]
    ruwheidhoog: Series[float] = pa.Field(gt=0, coerce=True)
    ruwheidlaag: Series[float] = pa.Field(gt=0, coerce=True)
    soortparameter: Series[str]
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)
    waarde: Series[float]


class PompSchema(PDBasicShema):
    code: Series[str] = pa.Field(coerce=True)  # addition to confluence
    gemaalid: Series[str] = pa.Field(unique=True)
    maximalecapaciteit: Series[float] = pa.Field(ge=0)


class ProfielgroepSchema(PDBasicShema):
    brugid: Optional[Series[str]] = pa.Field(nullable=True)
    stuwid: Optional[Series[str]] = pa.Field(nullable=True)


class ProfiellijnSchema(GPDBasicShema):
    profielgroepid: Series[str] = pa.Field(coerce=True, unique=True)


class ProfielpuntSchema(GPDBasicShema):
    code: Series[str] = pa.Field(coerce=True)
    profiellijnid: Series[str] = pa.Field(coerce=True)
    codevolgnummer: Series[int] = pa.Field(coerce=True)


class RegelmiddelSchema(GPDBasicShema):
    code: Series[str]  # addition to confluence
    kunstwerkopeningid: Series[str] = pa.Field(unique=True)
    overlaatonderlaat: Series[str]
    soortregelbaarheid: Series[int] = pa.Field(
        isin=list(MANAGEMENT_DEVICE_TYPES.values()), coerce=True
    )
    stuwid: Series[str] = pa.Field(unique=True)


class RuwheidsprofielSchema(BasicSchema):
    code: Series[str]  # addition to confluence
    geometry: GeoSeries = pa.Field(nullable=True)  # by design, to add to gpkg
    profielpuntid: Series[str] = pa.Field(coerce=True, unique=True)
    ruwheidhoog: Series[float] = pa.Field(gt=0, coerce=True)
    ruwheidlaag: Series[float] = pa.Field(gt=0, coerce=True)
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)


class SturingSchema(PDBasicShema):
    bovengrens: Series[float] = pa.Field(nullable=True)
    code: Series[str] = pa.Field(coerce=True)  # addition to confluence
    doelvariabele: Series[str]  # addition to confluence
    ondergrens: Series[float] = pa.Field(nullable=True)
    pompid: Series[str] = pa.Field(unique=True)
    streefwaarde: Series[float] = pa.Field(nullable=True)


class StuwSchema(GPDBasicShema):
    afvoercoefficient: Series[float] = pa.Field(coerce=True)
    code: Series[str] = pa.Field(coerce=True)
    soortstuw: Series[float] = pa.Field(
        isin=HYDAMO_WEIR_TYPES, coerce=True
    )  # addition to confluence


class WaterloopSchema(GPDBasicShema):
    code: Series[str] = pa.Field(coerce=True, unique=True)
    typeruwheid: Series[str] = pa.Field(isin=ROUGHNESS_MAPPING_LIST)  # addition to confluence


class DHydamoDataModel(BaseModel):
    """
    Datamodel that validates attributes using pydantic and pandera, using pandera SchemaModels.
    The former only checks for type (i.e. DataFrame), the latter checks the content of the DataFrames
    Ensures that data is compatible with DHydamo

    Attributes:
        brug (gpd.GeoDataFrame): geodataframe containing bridge data
        duiker (gpd.GeoDataFrame): geodataframe containing culvert data
        gemaal (gpd.GeoDataFrame): geodataframe containing pump station data
        hydroobject_normgp (gpd.GeoDataFrame): geodataframe containing norm profile id data
        kunstwerkopening (gpd.GeoDataFrame): geodataframe containing weir opening data
        normgeparamprofielwaarde (gpd.GeoDataFrame): geodataframe containing norm profile data
        pomp (gpd.GeoDataFrame): geodataframe containing pump data
        profielgroep (gpd.GeoDataFrame): geodataframe containing measured profile groups
        profiellijn (gpd.GeoDataFrame): geodataframe containing measured profile lines
        profielpunt (gpd.GeoDataFrame): geodataframe containing measured profile points
        regelmiddel (gpd.GeoDataFrame): geodataframe containing management device data for weirs
        ruwheidsprofiel (gpd.GeoDataFrame): geodataframe containing roughness data for measured profiles
        sturing (gpd.GeoDataFrame): geodataframe containing management data for pumps
        stuw (gpd.GeoDataFrame): geodataframe containing weir data
        waterloop (gpd.GeoDataFrame): geodataframe containing branch data
    """

    brug: Optional[DataFrame[BrugSchema]]
    duiker: Optional[DataFrame[DuikerSchema]]
    gemaal: Optional[DataFrame[GemaalSchema]]
    hydroobject_normgp: Optional[DataFrame[Hydroobject_normgpSchema]]
    kunstwerkopening: Optional[DataFrame[KunstwerkopeningSchema]]
    normgeparamprofielwaarde: Optional[DataFrame[NormgeparamprofielwaardeSchema]]
    pomp: Optional[DataFrame[PompSchema]]
    profielgroep: Optional[DataFrame[ProfielgroepSchema]]
    profiellijn: Optional[DataFrame[ProfiellijnSchema]]
    profielpunt: Optional[DataFrame[ProfielpuntSchema]]
    regelmiddel: Optional[DataFrame[RegelmiddelSchema]]
    rivier_profielen: Optional[Any]
    rivier_profielen_data: Optional[Any]
    rivier_profielen_ruwheid: Optional[Any]
    ruwheidsprofiel: Optional[DataFrame[RuwheidsprofielSchema]]
    sturing: Optional[DataFrame[SturingSchema]]
    stuw: Optional[DataFrame[StuwSchema]]
    waterloop: Optional[DataFrame[WaterloopSchema]]

    class Config:
        validate_assignment = True  # validates new attributes that are assigned

    def to_gpkg(self, output_gpkg: str) -> None:
        """
        Class method that saves fields that are not None to gpkg

        Args:
            output_gpkg (str): file location to write geopackge to

        Returns:
            None
        """
        for key, value in self.__dict__.items():
            if value is not None:
                try:
                    getattr(self, key).to_file(output_gpkg, layer=key, driver="GPKG")
                except:
                    print("\nfailed to save {}\n".format(key))
