from pathlib import Path
from typing import List

# import data_functions as daf
import geopandas as gpd
import numpy as np
from hydrolib.core.io.crosssection.models import CrossDefModel, CrossLocModel
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.core.io.friction.models import FrictionModel
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.structure.models import StructureModel
from hydrolib.dhydamo.converters.df2hydrolibmodel import Df2HydrolibModel
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh
from hydrolib.dhydamo.io.dimrwriter import DIMRWriter


def add_branches(features: List[str], gpkg_path: str, hydamo: HyDAMO) -> HyDAMO:
    hydamo.branches.read_gpkg_layer(gpkg_path, layer_name="waterloop", index_col="code")
    # Check for circular features in branches

    hydamo.branches_popped = hydamo.branches.copy()
    for _, branch in hydamo.branches.iterrows():

        start = branch.geometry.coords[0]
        end = branch.geometry.coords[-1]
        if start == end:
            code = branch.code
            hydamo.branches_popped = hydamo.branches_popped.drop(code)

    hydamo.branches = hydamo.branches_popped.set_geometry("geometry")

    if ("hydroobject_normgp" in features) and ("normgeparamprofielwaarde" in features):
        hydamo.param_profile.read_gpkg_layer(
            gpkg_path,
            layer_name="hydroobject_normgp",
            index_col="globalid",
        )
        hydamo.param_profile_values.read_gpkg_layer(
            gpkg_path,
            layer_name="normgeparamprofielwaarde",
            index_col="normgeparamprofielid",
        )
    return hydamo


def add_bridges(
    gpkg_path: str, hydamo: HyDAMO, default_height=10, max_snap_dist: float = 5
) -> HyDAMO:
    hydamo.bridges.read_gpkg_layer(gpkg_path, layer_name="brug", index_col="code")
    hydamo.bridges.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=max_snap_dist)
    hydamo.bridges.dropna(axis=0, inplace=True, subset=["branch_offset"])
    for i, bridge in hydamo.bridges.iterrows():
        hydamo.structures.add_bridge(
            id=bridge.code,
            name=bridge.code,
            length=bridge.lengte,
            branchid=bridge.branch_id,
            chainage=bridge.branch_offset,
            frictiontype=bridge.typeruwheid,
            csdefid=bridge.code,  # TODO Check influence of csdefid
            shift=default_height,  # TODO Validate the use of offset
            friction=bridge.ruwheid,
            inletlosscoeff=bridge.intreeverlies,
            outletlosscoeff=bridge.uittreeverlies,
        )

    return hydamo


def add_culverts(gpkg_path: str, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
    hydamo.culverts.read_gpkg_layer(gpkg_path, layer_name="duiker")
    hydamo.culverts.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=max_snap_dist)
    hydamo.culverts.dropna(axis=0, inplace=True, subset=["branch_offset"])
    for i, culvert in hydamo.culverts.iterrows():
        hydamo.structures.add_culvert(
            id=culvert.code,
            name=culvert.code,
            branchid=culvert.branch_id,
            chainage=culvert.branch_offset,
            leftlevel=culvert.hoogtebinnenonderkantbov,
            rightlevel=culvert.hoogtebinnenonderkantbene,
            length=culvert.lengte,
            inletlosscoeff=culvert.intreeverlies,
            outletlosscoeff=culvert.uittreeverlies,
            # crosssection=daf.get_crosssection_culvert_AGV(
            #     shape=culvert.vormkoker,
            #     height=culvert.hoogteopening,
            #     width=culvert.breedteopening,
            #     closed=1,
            # ),
            crosssection=eval(culvert.doorstroomopening),
            bedfrictiontype=culvert.typeruwheid,
            bedfriction=culvert.ruwheid,
        )

    return hydamo


def add_pumps(gpkg_path: str, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
    hydamo.pumpstations.read_gpkg_layer(gpkg_path, layer_name="gemaal")
    hydamo.pumps.read_gpkg_layer(gpkg_path, layer_name="pomp")
    hydamo.management.read_gpkg_layer(gpkg_path, layer_name="sturing")

    hydamo.pumpstations.snap_to_branch(
        hydamo.branches, snap_method="overal", maxdist=max_snap_dist
    )
    hydamo.pumps.drop(
        axis=0,
        index=hydamo.pumpstations[np.isnan(hydamo.pumpstations.branch_offset)].index,
        inplace=True,
    )
    hydamo.pumpstations.dropna(axis=0, inplace=True, subset=["branch_offset"])

    hydamo.structures.convert.pumps(
        hydamo.pumpstations, pumps=hydamo.pumps, management=hydamo.management
    )

    return hydamo


def add_weirs(gpkg_path: str, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
    hydamo.weirs.read_gpkg_layer(gpkg_path, layer_name="stuw")
    hydamo.opening.read_gpkg_layer(gpkg_path, layer_name="kunstwerkopening")
    hydamo.management_device.read_gpkg_layer(gpkg_path, layer_name="regelmiddel")

    hydamo.weirs.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=max_snap_dist)
    hydamo.weirs.dropna(axis=0, inplace=True, subset=["branch_offset"])

    hydamo.structures.convert.weirs(
        weirs=hydamo.weirs, opening=hydamo.opening, management_device=hydamo.management_device
    )

    return hydamo


def build_1D_model(
    features: List, fm: FMModel, hydamo: HyDAMO, max_snap_dist: float = 5
) -> HyDAMO:
    # Add observation points (empty for now)
    hydamo.observationpoints.add_points(
        [], [], locationTypes=["1d", "1d"], snap_distance=max_snap_dist
    )

    if "brug" in features:
        bbrug = True
    else:
        bbrug = False

    if "duiker" in features:
        bduiker = True
    else:
        bduiker = False

    if "gemaal" in features:
        bgemaal = True
    else:
        bgemaal = False

    if "stuw" in features:
        brstuw = True
        bustuw = True
    else:
        brstuw = False
        bustuw = False

    # Collect all structures in a Structures dataframe
    structures = hydamo.structures.as_dataframe(
        bridges=bbrug, culverts=bduiker, pumps=bgemaal, rweirs=brstuw, uweirs=bustuw
    )

    mesh.mesh1d_add_branches_from_gdf(
        fm.geometry.netfile.network,
        branches=hydamo.branches,
        branch_name_col="code",
        node_distance=20,
        max_dist_to_struc=3,
        structures=structures,
    )

    # Set the crosssections for the branches

    hydamo.crosssections.convert.profiles(
        crosssections=hydamo.profile,
        crosssection_roughness=hydamo.profile_roughness,
        # profile_groups=hydamo.profile_group, # skip as no bridges or weirs with profile
        profile_lines=hydamo.profile_line,
        param_profile=hydamo.param_profile,
        param_profile_values=hydamo.param_profile_values,
        branches=hydamo.branches,
        roughness_variant="High",
    )

    # Set a default cross section
    default = hydamo.crosssections.add_rectangle_definition(
        height=5.0,
        width=5.0,
        closed=False,
        roughnesstype="StricklerKs",
        roughnessvalue=30,
        name="default",
    )
    hydamo.crosssections.set_default_definition(definition=default, shift=10.0)
    return fm, hydamo


def build_2D_model(
    dx: float,
    dy: float,
    extent: gpd.GeoDataFrame,
    fm: FMModel,
    elevation_raster_path: str = None,
    one_d=False,
) -> FMModel:

    network = fm.geometry.netfile.network

    mesh.mesh2d_add_rectilinear(network=network, polygon=extent.geometry.values[0], dx=dx, dy=dy)

    if elevation_raster_path is not None:
        mesh.mesh2d_altitude_from_raster(
            network=network,
            rasterpath=elevation_raster_path,
            where="node",  # Face does not work
            stat="nanmean",
            fill_option="fill_value",
            fill_value=-10,
        )

    if one_d:
        mesh.links1d2d_add_links_1d_to_2d(network=network)
        # mesh.links1d2d_add_links_2d_to_1d_embedded(network=network)
    return fm


def simple_fm_model(start_time: int = 20160601, stop_time: int = 2 * 86400) -> FMModel:
    fm = FMModel()
    fm.time.refdate = start_time
    fm.time.tstop = stop_time

    return fm


def write_model(fm: FMModel, hydamo: HyDAMO, output_folder: str, one_d=True):
    if one_d:
        models = Df2HydrolibModel(hydamo)
        # Export to DIMR configuration
        fm.geometry.structurefile = [StructureModel(structure=models.structures)]
        fm.geometry.crosslocfile = CrossLocModel(crosssection=models.crosslocs)
        fm.geometry.crossdeffile = CrossDefModel(definition=models.crossdefs)

        fm.geometry.frictfile = []
        for i, fric_def in enumerate(models.friction_defs):
            fric_model = FrictionModel(global_=fric_def)
            fric_model.filepath = f"roughness_{i}.ini"
            fm.geometry.frictfile.append(fric_model)

    # fm.output.obsfile = [ObservationPointModel(observationpoint=models.obspoints)]

    # extmodel = ExtModel()
    # extmodel.boundary = models.boundaries_ext
    # extmodel.lateral = models.laterals_ext
    # fm.external_forcing.extforcefilenew = extmodel

    # fm.geometry.inifieldfile = IniFieldModel(initial=models.inifields)
    # for ifield, onedfield in enumerate(models.onedfieldmodels):
    #     fm.geometry.inifieldfile.initial[ifield].datafile = OneDFieldModel(
    #         global_= onedfield
    #     )
    # Now we write the file structure:
    output_path = Path(output_folder)
    output_path.mkdir(exist_ok=True, parents=True)
    fm.filepath = Path(output_folder) / "fm" / "test.mdu"
    dimr = DIMR()
    dimr.component.append(
        FMComponent(name="DFM", workingDir=output_path / "fm", model=fm, inputfile=fm.filepath)
    )
    dimr.save(recurse=True)
    # import shutil
    # shutil.copy(data_path / "initialWaterDepth.ini", folder / "fm")

    dimr = DIMRWriter(output_path=output_path)
    dimr.write_dimrconfig(fm=fm)  # , rr_model=drrmodel, rtc_model=drtcmodel)
    dimr.write_runbat()
