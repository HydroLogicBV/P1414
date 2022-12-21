#%%
import importlib
from pathlib import Path

import contextily as ctx
import matplotlib.pyplot as plt
import numpy as np
from hydrolib.core.io.crosssection.models import CrossDefModel, CrossLocModel
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.core.io.friction.models import FrictionModel
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.structure.models import StructureModel
from hydrolib.dhydamo.converters.df2hydrolibmodel import Df2HydrolibModel
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh
from hydrolib.dhydamo.io.common import ExtendedDataFrame, ExtendedGeoDataFrame
from hydrolib.dhydamo.io.dimrwriter import DIMRWriter
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from tqdm import tqdm

import data_functions as daf
from data_functions import *

# importlib.reload(data_functions)
#%% ############################################
# Data laden

# folder = r"D:\work\P1414_ROI"
folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HDSR\HDSR_hydamo.gpkg"
norm_profile_gpkg = folder + r"\GIS\HDSR\norm_profielen_test.gpkg"

extent_shp_path = folder + "\GIS\WAGV\AGV_mask.shp"
twod_depth_path = "D:\Work\Project\P1414\GIS\WAGV\AGV_dummy_depth_v2.tif"

hydamo = HyDAMO()
hydamo.branches.read_gpkg_layer(gpkg_file, layer_name="waterloop", index_col="code")
hydamo.bridges.read_gpkg_layer(gpkg_file, layer_name="brug", index_col="code")
hydamo.culverts.read_gpkg_layer(gpkg_file, layer_name="duiker")
hydamo.weirs.read_gpkg_layer(gpkg_file, layer_name="stuw")
hydamo.opening.read_gpkg_layer(gpkg_file, layer_name="kunstwerkopening")
hydamo.management_device.read_gpkg_layer(gpkg_file, layer_name="regelmiddel")
hydamo.pumpstations.read_gpkg_layer(gpkg_file, layer_name="gemaal")
hydamo.pumps.read_gpkg_layer(gpkg_file, layer_name="pomp")
hydamo.management.read_gpkg_layer(gpkg_file, layer_name="sturing")

# %% ##################################################
# Check for circular features in branches

circular_branches = []
hydamo.branches_popped = hydamo.branches.copy()
j = 0
for i, branch in hydamo.branches.iterrows():

    start = branch.geometry.coords[0]
    end = branch.geometry.coords[-1]
    if start == end:
        code = branch.code
        hydamo.branches_popped = hydamo.branches_popped.drop(code)

hydamo.branches = hydamo.branches_popped.set_geometry("geometry")


#%% #######################################################
# Snap profiles to branch
hydamo.param_profile.read_gpkg_layer(
    norm_profile_gpkg,
    layer_name="hydroobject_normgp",
    index_col="globalid",
)
hydamo.param_profile_values.read_gpkg_layer(
    norm_profile_gpkg,
    layer_name="normgeparamprofielwaarde",
    index_col="normgeparamprofielid",
)

# hydamo.profile_roughness.read_gpkg_layer(profile_gpkg, layer_name="ruwheidsprofiel")
# hydamo.profile.snap_to_branch(hydamo.branches, snap_method="intersecting")
# hydamo.profile.dropna(axis=0, inplace=True, subset=["branch_offset"])
# hydamo.profile_line.read_gpkg_layer(profile_gpkg, layer_name="profiellijn")
# hydamo.profile_group.read_gpkg_layer(profile_gpkg, layer_name="profielgroep")
# hydamo.profile.drop("code", axis=1, inplace=True)
# hydamo.profile["code"] = hydamo.profile["profiellijnid"]


# Snap structures to branches
hydamo.bridges.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=10)
hydamo.bridges.dropna(axis=0, inplace=True, subset=["branch_offset"])

hydamo.culverts.snap_to_branch(hydamo.branches, snap_method="ends", maxdist=15)
hydamo.culverts.dropna(axis=0, inplace=True, subset=["branch_offset"])  # Without branch
hydamo.culverts.dropna(axis=0, inplace=True, subset=["hoogteopening"])  # Without dimensions
hydamo.culverts.drop(
    hydamo.culverts[hydamo.culverts["vormkoker"] != (1 or 3)].index, inplace=True
)  # Without unknown shapes (i.e. only circle and rectangle)

hydamo.weirs.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=10)
hydamo.weirs.dropna(axis=0, inplace=True, subset=["branch_offset"])

hydamo.pumpstations.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=10)
hydamo.pumps.drop(
    axis=0,
    index=hydamo.pumpstations[np.isnan(hydamo.pumpstations.branch_offset)].index,
    inplace=True,
)
hydamo.pumpstations.dropna(axis=0, inplace=True, subset=["branch_offset"])

#%% ######################################################
# Add structures
hydamo.structures.convert.weirs(
    weirs=hydamo.weirs, opening=hydamo.opening, management_device=hydamo.management_device
)

hydamo.structures.convert.pumps(
    hydamo.pumpstations, pumps=hydamo.pumps, management=hydamo.management
)


# Add structures in a loop to provide the right format
for i, bridge in tqdm(hydamo.bridges.iterrows()):
    hydamo.structures.add_bridge(
        id=bridge.code,
        name=bridge.code,
        length=bridge.lengte,
        branchid=bridge.branch_id,
        chainage=bridge.branch_offset,
        frictiontype=daf.get_roughness(bridge.typeruwheid),
        csdefid=bridge.code,  # TODO Check influence of csdefid
        shift=bridge.branch_offset,  # TODO Validate the use of offset
        friction=bridge.ruwheid,
        inletlosscoeff=bridge.intreeverlies,
        outletlosscoeff=bridge.uittreeverlies,
    )

for i, culvert in tqdm(hydamo.culverts.iterrows()):
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
        crosssection=daf.get_crosssection_culvert_AGV(
            shape=culvert.vormkoker,
            height=culvert.hoogteopening,
            width=culvert.breedteopening,
            closed=1,
        ),
        bedfrictiontype=get_roughness(culvert.typeruwheid),
        bedfriction=culvert.ruwheid,
    )


#%% #####################################################
# Finalize the input to the dhydro model
# Add observation points (empty for now)
hydamo.observationpoints.add_points([], [], locationTypes=["1d", "1d"], snap_distance=10.0)

# Collect all structures in a Structures dataframe
structures = hydamo.structures.as_dataframe(
    rweirs=True,
    bridges=True,
    uweirs=True,
    culverts=True,
    orifices=False,
    pumps=False,
)

#%% ####################################################
# Create the FM model
fm = FMModel()
# Set start and stop time
fm.time.refdate = 20160601
fm.time.tstop = 2 * 3600 * 24

mesh.mesh1d_add_branches_from_gdf(
    fm.geometry.netfile.network,
    branches=hydamo.branches,
    branch_name_col="code",
    node_distance=20,
    max_dist_to_struc=3,
    structures=structures,
)


#%% #################################################################
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

# %% 2D
twod = False
if twod:
    print("doing 2D")
    extent = gpd.read_file(extent_shp_path)
    network = fm.geometry.netfile.network

    mesh.mesh2d_add_rectilinear(network=network, polygon=extent.geometry.values[0], dx=100, dy=100)

    mesh.links1d2d_add_links_1d_to_2d(network=network)
    mesh.mesh2d_altitude_from_raster(
        network=network,
        rasterpath=twod_depth_path,
        where="node",  # Face does not work
        stat="mean",
        fill_option="fill_value",
        fill_value=100,
    )

    print("2D done")
# %% ####################################################

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


fm.filepath = Path(folder) / "fm" / "test.mdu"
dimr = DIMR()
dimr.component.append(
    FMComponent(name="DFM", workingDir=Path(folder) / "fm", model=fm, inputfile=fm.filepath)
)
dimr.save(recurse=True)
# import shutil
# shutil.copy(data_path / "initialWaterDepth.ini", folder / "fm")

dimr = DIMRWriter(output_path=folder)
dimr.write_dimrconfig(fm)  # , rr_model=drrmodel, rtc_model=drtcmodel)
dimr.write_runbat()

# %%
