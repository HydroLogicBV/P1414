# generated by datamodel-codegen:
#   filename:  rtcDataConfig.json
#   timestamp: 2022-09-23T08:50:37+00:00

from __future__ import annotations

from enum import Enum
from typing import Any, List, Optional, Union

from pydantic import BaseModel, Extra, Field, conint, constr

# from . import _


class RtcAggregationTypeEnumStringType(Enum):
    BLOCK = "BLOCK"
    LINEAR = "LINEAR"


class RtcEnsembleModeEnumStringType(Enum):
    JOINT = "JOINT"
    TREE = "TREE"
    INDEPENDENT = "INDEPENDENT"


class RtcExternalBooleanSimpleType(BaseModel):
    __root__: Union[bool, constr(regex=r"^([\$][\(-_a-z]+[\$])$")]


class RtcExternalIntegerSimpleType(BaseModel):
    __root__: Union[int, constr(regex=r"^([\$][\(-_a-z]+[\$])$")]


class RtcExternalParameterSimpleType(BaseModel):
    __root__: Union[float, constr(regex=r"^([#-\$][\(-_a-z]+[#-\$])$")]


class RtcPIExtrapolationOptionEnumStringType(Enum):
    BLOCK = "BLOCK"
    PERIODIC = "PERIODIC"


class RtcPIInterpolationOptionEnumStringType(Enum):
    BLOCK = "BLOCK"
    LINEAR = "LINEAR"


class RtcExtrapolationOption(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[RtcPIExtrapolationOptionEnumStringType] = Field(None, alias="$")


class RtcInterpolationOption(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[RtcPIInterpolationOptionEnumStringType] = Field(None, alias="$")


class _Validation(Enum):
    NO = "NO"
    STATE = "STATE"
    UPDATE = "UPDATE"
    UPDATE_EXCEPT_STATE = "UPDATE_EXCEPT_STATE"
    FORECAST = "FORECAST"
    FORECAST_EXCEPT_T0 = "FORECAST_EXCEPT_T0"
    ALL = "ALL"
    ALL_EXCEPT_STATE = "ALL_EXCEPT_STATE"


class RtcSeparatorEnumStringType(Enum):
    _ = "."
    __1 = ","
    __2 = ";"


class RtcTimeSeriesSimpleType(BaseModel):
    __root__: constr(min_length=1)


class RtcTimeZoneSimpleType(BaseModel):
    __root__: float = Field(
        ...,
        description="The timeZone (in decimal hours shift from GMT)\n            e.g. -1.0 or 3.5. If not present GMT is assumed",
    )


class RtcUnitEnumStringType(Enum):
    m = "m"
    m_2 = "m^2"
    m_3 = "m^3"
    m_3_s = "m^3/s"
    s = "s"


class RtcVariableTypeEnumStringType(Enum):
    CONTINUOUS = "CONTINUOUS"
    INTEGER = "INTEGER"
    TIMEINSTANCE = "TIMEINSTANCE"


class RtcDateType(BaseModel):
    __root__: constr(regex=r"^([\d][\d][\d][\d]\-[\d][\d]\-[\d][\d])$")


class RtcTimeSeriesType(Enum):
    accumulative = "accumulative"
    instantaneous = "instantaneous"


class RtcTimeStepUnitEnumStringType(Enum):
    second = "second"
    minute = "minute"
    hour = "hour"
    day = "day"
    week = "week"


class RtcTimeType(BaseModel):
    __root__: constr(regex=r"^([\d][\d]\:[\d][\d]\:[\d][\d])$")


class XsBoolean(BaseModel):
    __root__: bool


class XsPositiveInteger(BaseModel):
    __root__: conint(ge=1)


class XsString(BaseModel):
    __root__: str


class RtcCSVTimeSeriesFileComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    _adjointOutput: Optional[XsBoolean] = Field(None, alias="@adjointOutput")
    _decimalSeparator: Optional[RtcSeparatorEnumStringType] = Field(
        None, alias="@decimalSeparator"
    )
    _delimiter: Optional[RtcSeparatorEnumStringType] = Field(None, alias="@delimiter")


class RtcDateTimeComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    _date: RtcDateType = Field(..., alias="@date")
    _time: RtcTimeType = Field(..., alias="@time")


class RtcElementId(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcQuantityId(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcUnit(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[RtcUnitEnumStringType] = Field(None, alias="$")


class RtcOpenMIExchangeItemComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_elementId: RtcElementId = Field(
        ...,
        alias="rtc:elementId",
        description="OpenMI element ID, corresponds to the locationId",
    )
    rtc_quantityId: RtcQuantityId = Field(
        ...,
        alias="rtc:quantityId",
        description="OpenMI quantity ID, corresponds to the parameterId",
    )
    rtc_unit: RtcUnit = Field(
        ..., alias="rtc:unit", description="Selection of supported units"
    )


class RtcLocationId(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcParameterId(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcQualifierIdItem(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcUnit1(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcAdjointOutput(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsBoolean] = Field(None, alias="$")


class RtcTimeSeriesFile(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcUseBinFile(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsBoolean] = Field(None, alias="$")


class RtcPITimeSeriesExportFileComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_adjointOutput: Optional[RtcAdjointOutput] = Field(
        None, alias="rtc:adjointOutput"
    )
    rtc_timeSeriesFile: RtcTimeSeriesFile = Field(
        ...,
        alias="rtc:timeSeriesFile",
        description="Name of the file containing timeseries data. ",
    )
    rtc_useBinFile: Optional[RtcUseBinFile] = Field(
        None,
        alias="rtc:useBinFile",
        description='When true the events in the PI time series file are read from / written into a binairy file instead of the xml file.\nThe xml file only contains the time series headers and optionally a time zone.\nThe binairy file has the same name as the xml file only the extension is "bin" instead of "xml". The byte order in the bin file is always Intel x86.\n                ',
    )


class RtcTimeSeriesFile1(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsString] = Field(None, alias="$")


class RtcUseBinFile1(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    __1: Optional[XsBoolean] = Field(None, alias="$")


class RtcPITimeSeriesImportFileComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_timeSeriesFile: RtcTimeSeriesFile1 = Field(
        ...,
        alias="rtc:timeSeriesFile",
        description="Name of the file containing timeseries data. ",
    )
    rtc_useBinFile: Optional[RtcUseBinFile1] = Field(
        None,
        alias="rtc:useBinFile",
        description="OBSOLETE. Still here for backwards compatibility. Remove after next release.",
    )


class RtcTimeStepComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    _divider: Optional[XsPositiveInteger] = Field(None, alias="@divider")
    _multiplier: Optional[XsPositiveInteger] = Field(None, alias="@multiplier")
    _unit: RtcTimeStepUnitEnumStringType = Field(..., alias="@unit")


class RtcPITimeSeriesComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_extrapolationOption: Optional[RtcExtrapolationOption] = Field(
        None,
        alias="rtc:extrapolationOption",
        description="Extrapolation option in data import",
    )
    rtc_interpolationOption: Optional[RtcInterpolationOption] = Field(
        None,
        alias="rtc:interpolationOption",
        description="Interpolation option in data import",
    )
    rtc_locationId: RtcLocationId = Field(
        ..., alias="rtc:locationId", description="Location ID in Delft-FEWS PI-XML file"
    )
    rtc_parameterId: RtcParameterId = Field(
        ...,
        alias="rtc:parameterId",
        description="Parameter ID in Delft-FEWS PI-XML file",
    )
    rtc_qualifierId: Optional[List[RtcQualifierIdItem]] = Field(
        None, alias="rtc:qualifierId"
    )
    rtc_timeStep: Optional[RtcTimeStepComplexType] = Field(
        None,
        alias="rtc:timeStep",
        description="Equidistant time step of time series with optional multiplier of divider",
    )
    rtc_unit: Optional[RtcUnit1] = Field(
        None,
        alias="rtc:unit",
        description="Optional check for this unit during import, write this unit optionally when export the time series",
    )


class RtcRTCTimeSeriesComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    _id: constr(min_length=1) = Field(..., alias="@id")
    _validation: Optional[_Validation] = Field(None, alias="@validation")
    _vectorLength: Optional[conint(ge=1, le=2147483647)] = Field(
        None, alias="@vectorLength"
    )
    rtc_OpenMIExchangeItem: Optional[RtcOpenMIExchangeItemComplexType] = Field(
        None,
        alias="rtc:OpenMIExchangeItem",
        description="Time series definition of the OpenMI format for the online coupling of models during runtime",
    )
    rtc_PITimeSeries: Optional[RtcPITimeSeriesComplexType] = Field(
        None,
        alias="rtc:PITimeSeries",
        description="Time series definition of the PI XML time series format of Delft-FEWS",
    )


class RtcRTCSeriesExportComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_CSVTimeSeriesFile: Optional[RtcCSVTimeSeriesFileComplexType] = Field(
        None,
        alias="rtc:CSVTimeSeriesFile",
        description="Comma-separated file for data exports. Note that this option is only used in the exportSeries element. If selected, all available time series will be exported.",
    )
    rtc_PITimeSeriesFile: Optional[RtcPITimeSeriesExportFileComplexType] = Field(
        None, alias="rtc:PITimeSeriesFile"
    )
    rtc_timeSeries: List[RtcRTCTimeSeriesComplexType] = Field(
        ..., alias="rtc:timeSeries", min_items=1
    )


class RtcRTCSeriesImportComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_PITimeSeriesFile: Optional[RtcPITimeSeriesImportFileComplexType] = Field(
        None, alias="rtc:PITimeSeriesFile"
    )
    rtc_timeSeries: List[RtcRTCTimeSeriesComplexType] = Field(
        ..., alias="rtc:timeSeries", min_items=1
    )


class RtcRTCDataConfigComplexType(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    rtc_exportSeries: RtcRTCSeriesExportComplexType = Field(
        ...,
        alias="rtc:exportSeries",
        description="Export time series RTC-Tools genenerates and exports to XML or csv files or supplies to other applications via other interfaces",
    )
    rtc_importSeries: RtcRTCSeriesImportComplexType = Field(
        ...,
        alias="rtc:importSeries",
        description="Import time series RTC-Tools imports from XML files or other interfaces",
    )


class Model(BaseModel):
    class Config:
        extra = Extra.forbid

    _: Optional[str] = Field(None, alias="#")
    _xmlns_rtc: Optional[Any] = Field("http://www.wldelft.nl/fews", alias="@xmlns:rtc")
    _xmlns_xs: Optional[Any] = Field(
        "http://www.w3.org/2001/XMLSchema", alias="@xmlns:xs"
    )
    rtc_rtcDataConfig: Optional[_.RtcRtcDataConfig] = Field(
        None, alias="rtc:rtcDataConfig"
    )
