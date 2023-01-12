import importlib
from pathlib import Path
from typing import List, Union

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
from shapely.geometry import Polygon

from data_structures.dhydamo_data_model import DHydamoDataModel


def add_branches(
    ddm: DHydamoDataModel, features: List[str], gpkg_path: str, hydamo: HyDAMO
) -> HyDAMO:
    hydamo.branches.set_data(ddm.waterloop, index_col="code")

    # Check for circular features in branches
    hydamo.branches_popped = hydamo.branches.copy()
    for _, branch in hydamo.branches.iterrows():

        start = branch.geometry.coords[0]
        end = branch.geometry.coords[-1]
        if start == end:
            code = branch.code
            hydamo.branches_popped = hydamo.branches_popped.drop(code)

    hydamo.branches = hydamo.branches_popped.set_geometry("geometry")

    if (
        ("profielgroep" in features)
        and ("profielpunt" in features)
        and ("profiellijn" in features)
        and ("ruwheidsprofiel" in features)
    ):
        hydamo.profile.read_gpkg_layer(
            gpkg_path,
            layer_name="profielpunt",
            groupby_column="profiellijnid",
            order_column="codevolgnummer",
            index_col="code",
        )

        # Snap profiles to branch
        # hydamo.profile_roughness.read_gpkg_layer(gpkg_path, layer_name="ruwheidsprofiel")
        hydamo.profile_roughness.set_data(ddm.ruwheidsprofiel)
        hydamo.profile.snap_to_branch(hydamo.branches, snap_method="intersecting")
        hydamo.profile.dropna(axis=0, inplace=True, subset=["branch_offset"])
        # hydamo.profile_line.read_gpkg_layer(gpkg_path, layer_name="profiellijn")
        # hydamo.profile_group.read_gpkg_layer(gpkg_path, layer_name="profielgroep")
        hydamo.profile_line.set_data(ddm.profiellijn)
        hydamo.profile_group.set_data(ddm.profielgroep)

        hydamo.profile.drop("code", axis=1, inplace=True)
        hydamo.profile["code"] = hydamo.profile["profiellijnid"]

    if ("hydroobject_normgp" in features) and ("normgeparamprofielwaarde" in features):
        hydamo.param_profile.set_data(ddm.hydroobject_normgp, index_col="globalid")
        hydamo.param_profile_values.set_data(
            ddm.normgeparamprofielwaarde, index_col="normgeparamprofielid"
        )
    return hydamo


def add_bridges(
    ddm: DHydamoDataModel, hydamo: HyDAMO, default_height=10, max_snap_dist: float = 5
) -> HyDAMO:

    hydamo.bridges.set_data(ddm.brug, index_col="code")
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


def add_culverts(ddm: DHydamoDataModel, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
    hydamo.culverts.set_data(ddm.duiker)
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
            crosssection=eval(
                culvert.doorstroomopening
            ),  # dict saved as string, so eval turns it into dict
            bedfrictiontype=culvert.typeruwheid,
            bedfriction=culvert.ruwheid,
        )

    return hydamo


def add_pumps(ddm: DHydamoDataModel, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
    hydamo.pumpstations.set_data(ddm.gemaal)
    # hydamo.pumps.set_data(ddm.pomp)
    hydamo.pumpstations.snap_to_branch(
        hydamo.branches, snap_method="overal", maxdist=max_snap_dist
    )
    hydamo.pumpstations.dropna(axis=0, inplace=True, subset=["branch_offset"])

    hydamo.pumps.set_data(ddm.pomp[ddm.pomp["gemaalid"].isin(hydamo.pumpstations["globalid"])])
    hydamo.management.set_data(ddm.sturing)

    hydamo.structures.convert.pumps(
        hydamo.pumpstations, pumps=hydamo.pumps, management=hydamo.management
    )

    return hydamo


def add_weirs(ddm: DHydamoDataModel, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
    hydamo.weirs.set_data(ddm.stuw)
    hydamo.opening.set_data(ddm.kunstwerkopening)
    hydamo.management_device.set_data(ddm.regelmiddel)

    hydamo.weirs.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=max_snap_dist)
    hydamo.weirs.dropna(axis=0, inplace=True, subset=["branch_offset"])

    hydamo.structures.convert.weirs(
        weirs=hydamo.weirs, opening=hydamo.opening, management_device=hydamo.management_device
    )

    return hydamo


def build_1D_model(
    features: List,
    fm: FMModel,
    hydamo: HyDAMO,
    max_snap_dist: float = 5,
    node_distance=20,
    max_dist_to_struct=3,
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
        node_distance=node_distance,
        max_dist_to_struc=max_dist_to_struct,
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
    coupling: str = "1Dto2D",
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
        if coupling == "1Dto2D":
            mesh.links1d2d_add_links_1d_to_2d(network=network)
        elif coupling == "2Dto1D":
            mesh.links1d2d_add_links_2d_to_1d_embedded(network=network)
    return fm


def simple_fm_model(start_time: int = 20160601, stop_time: int = 2 * 86400) -> FMModel:
    fm = FMModel()
    fm.time.refdate = start_time
    fm.time.tstop = stop_time

    return fm


def to_dhydro(
    self, config: str, output_folder: str, extent: Union[gpd.GeoDataFrame, Polygon] = None
):
    """ """
    # check if data has been loaded and correct attributes are set
    if (not hasattr(self, "ddm")) | (not hasattr(self, "features")):
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
        if model_config.FM.one_d_bool:
            # check if there are branches in the data
            if "waterloop" not in self.features:
                raise ValueError("Missing branches")

            # add branches, profiles and norm profiles
            print("\nworking on 1D branches\n")
            self.hydamo = add_branches(
                ddm=self.ddm, features=self.features, gpkg_path=self.gpkg_path, hydamo=self.hydamo
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
                        ddm=self.ddm,
                        hydamo=self.hydamo,
                        max_snap_dist=model_config.FM.one_d.max_snap_dist,
                    )

            # add 1D model
            print("\nCompletig 1D model\n")
            self.fm, self.hydamo = build_1D_model(
                fm=self.fm,
                features=self.features,
                hydamo=self.hydamo,
                node_distance=model_config.FM.one_d.node_distance,
                max_dist_to_struct=model_config.FM.one_d.max_dist_to_struct,
            )

        # build 2D model
        if model_config.FM.two_d_bool:
            # compute extent of 1D network
            # Add to gpgk?
            if extent is None:
                extent = self.ddm.waterloop.dissolve(by=None).convex_hull.buffer(
                    model_config.FM.two_d.two_d_buffer
                )

            # add 2D model
            print("\nBuilding 2D model grid\n")
            self.fm = build_2D_model(
                coupling=model_config.FM.two_d.coupling_type,
                dx=model_config.FM.two_d.dx,
                dy=model_config.FM.two_d.dy,
                elevation_raster_path=model_config.FM.two_d.elevation_raster_path,
                extent=extent,
                fm=self.fm,
                one_d=model_config.FM.one_d_bool,
            )

    # save D-HYDRO model
    write_model(self.fm, self.hydamo, output_folder=output_folder, one_d=model_config.FM.one_d)


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
