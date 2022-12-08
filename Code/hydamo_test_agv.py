import importlib
from pathlib import Path

import contextily as ctx
import matplotlib.pyplot as plt
import numpy as np
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.structure.models import *
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

hydamo = HyDAMO(extent_file=folder + "\GIS\WAGV\AGV_mask.shp")
hydamo.branches.read_shp(
    path=folder + r"\GIS\Uitgesneden watergangen\AGV.shp",
    column_mapping={
        "ruwheidsty": "typeruwheid",
    },
    index_col="code",
)

hydamo.bridges.read_shp(
    path=folder + r"\GIS\WAGV\brug_v13\brug_v13_clipped.shp",
    column_mapping={
        "naam": "globalid",
        "ruwheidsty": "typeruwheid",
        "intreeverl": "intreeverlies",
        "uittreever": "uittreeverlies",
    },
    index_col="code",
)

hydamo.culverts.read_shp(
    path=folder + r"\GIS\WAGV\duikersifonhevel_v13\duikersifonhevel_v13_clipped.shp",
    column_mapping={
        "naam": "globalid",
        "hoogtebinn": "hoogtebinnenonderkantbene",
        "hoogtebin0": "hoogtebinnenonderkantbov",
        "vormkokeri": "vormkoker",
        "hoogteopen": "hoogteopening",
        "breedteope": "breedteopening",
        "intreeverl": "intreeverlies",
        "uittreever": "uittreeverlies",
        "ruwheidsty": "typeruwheid",
    },
    index_col="code",
)

hydamo.weirs.read_shp(
    path=folder + r"\GIS\WAGV\stuw_v13\stuw_v13_clipped.shp",
    column_mapping={
        "code": "globalid",
        "naam": "code",
        "afvoercoef": "afvoercoefficient",
        "soortstuwi": "soortstuw",
    },
    index_col="globalid",
)
hydamo.opening.read_gpkg_layer(
    gpkg_path=folder + r"\GIS\WAGV\doorstroomopening_stuw_v13_clipped.gpkg",
    layer_name="doorstroomopening_stuw_v13_clipped",
    column_mapping={
        "laagstedoo": "laagstedoorstroombreedte",
        "laagstedo0": "laagstedoorstroomhoogte",
        "stuwcode": "stuwid",
        "code": "globalid",
        "vormopenin": "vormopening",
        "afvoercoef": "afvoercoefficient",
        "hoogstedoo": "hoogstedoorstroombreedte",
        "hoogstedo0": "hoogstedoorstroomhoogte",
    },
)

hydamo.pumpstations.read_shp(
    path=folder + r"\GIS\WAGV\pomp_gemaal_v13\pomp_gemaal_v13_clipped_streefpeil.shp",
    column_mapping={"gemaalcode": "globalid"},
)

# hydamo.profile.read_shp(
#     path = folder + r"\GIS\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped.shp",
#     column_mapping =   {#"code":"globalid",
#                         "metingprof":"profiellijnid",
#                         "codevolgnu":"codevolgnummer",
#                         },
#     index_col = "code",
# )

# hydamo.profile_roughness.read_shp(
#     path = folder + r"\GIS\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped.shp",
#     column_mapping =   {#"code":"profielpuntid",
#                         "ruwheidsty":"typeruwheid",
#                         "ruwheidswa":"ruwheidlaag",
#                         "ruwheidsw0":"ruwheidhoog",
#                         },
#     index_col = "code",
# )


#%%  ################################################
# Test plot the features

fig = plt.figure()
ax = plt.gca()
hydamo.branches.geometry.plot(ax=ax, label="Channels", linewidth=0.5, color="black")
hydamo.culverts.geometry.plot(ax=ax, label="Culverts", linewidth=2, color="blue")
hydamo.bridges.geometry.plot(ax=ax, label="Bridges", markersize=2, color="red")
hydamo.weirs.geometry.plot(ax=ax, label="Weirs", markersize=2, color="orange")
ctx.add_basemap(ax, crs=28992, source=ctx.providers.OpenStreetMap.Mapnik)
plt.show()

# print(hydamo)


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

    # if branch.code is None:
    #     # hydamo.branches_popped.loc[i, "code"] = "opgevuld {}".format(j)
    #     j += 1

# remove branches with none as index
# try:
#     hydamo.branches_popped = hydamo.branches_popped.drop([None], axis=0)
# except KeyError:
#     pass
# print(hydamo.branches_popped)
# hydamo.branches_popped.dropna(axis="geometry", inplace=True)
hydamo.branches = hydamo.branches_popped.set_geometry("geometry")


#%% #######################################################
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
hydamo.pumpstations.dropna(axis=0, inplace=True, subset=["branch_offset"])

# Define pumpinformation in a separate ExtendedDataFrame
pumps = ExtendedDataFrame()
pumps["gemaalid"] = hydamo.pumpstations["globalid"]
pumps["maximalecapaciteit"] = hydamo.pumpstations["maximaleca"]
pumps["globalid"] = hydamo.pumpstations["code"]
pumps["code"] = daf.getuniquecode("Pumps", len(hydamo.pumpstations["globalid"]))

# Define managementinformation in a separate ExtendedDataFrame
sturing = ExtendedDataFrame()
sturing["pompid"] = hydamo.pumpstations["code"]
sturing["streefwaarde"] = hydamo.pumpstations["peil1"]
sturing["ondergrens"] = hydamo.pumpstations["peil1"] - 0.05  # Very crude #TODO fix this conversion
sturing["bovengrens"] = hydamo.pumpstations["peil1"] + 0.05  # Very crude #TODO fix this conversion
sturing["doelvariabele"] = "waterstand"

# Add information on kunstwerkopening to opening
hydamo.opening["kunstwerkopeningid"] = hydamo.opening["globalid"]
hydamo.opening["overlaatonderlaat"] = "Overlaat"

# Drop duplicates in globalid of pumpstations and covert back to ExtendedGeoDataFrame
subset_pumps = hydamo.pumpstations.drop_duplicates(subset=["globalid"])
pumpgeotype = hydamo.pumpstations.geotype
hydamo.pumpstations = ExtendedGeoDataFrame(geotype=pumpgeotype)
hydamo.pumpstations.set_data(subset_pumps)

#%% ######################################################
# Add structures
hydamo.structures.convert.weirs(
    weirs=hydamo.weirs, opening=hydamo.opening, management_device=hydamo.opening
)

hydamo.structures.convert.pumps(hydamo.pumpstations, pumps=pumps, management=sturing)

# Define a list of possible roughness types
roughness_list = [
    "Bos en Bijkerk",
    "Chezy",
    "Manning",
    "StricklerKn",
    "StricklerKs",
    "White Colebrook",
]

# Add structures in a loop to provide the right format
for i, bridge in enumerate(tqdm(hydamo.bridges.code)):
    hydamo.structures.add_bridge(
        id=hydamo.bridges.code[i],
        name=hydamo.bridges.code[i],
        length=hydamo.bridges.lengte[i],
        branchid=hydamo.bridges.branch_id[i],
        chainage=hydamo.bridges.branch_offset[i],  # TODO Validate the use of offset
        frictiontype=roughness_list[hydamo.bridges.typeruwheid[i]],
        csdefid=hydamo.bridges.code[i],  # TODO Check influence of csdefid
        shift=hydamo.bridges.branch_offset[i],  # TODO Validate the use of offset
        friction=hydamo.bridges.ruwheid[i],
        inletlosscoeff=hydamo.bridges.intreeverlies[i],
        outletlosscoeff=hydamo.bridges.uittreeverlies[i],
    )

for i, culvert in tqdm(hydamo.culverts.iterrows()):
    hydamo.structures.add_culvert(
        id=culvert.code,
        name=culvert.code,
        branchid=culvert.branch_id,
        chainage=culvert.branch_offset,  # TODO Valdiate the use of offset
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
        ),  # TODO Check for automatic cross section based on shape
        numlosscoeff=None,  # TODO Check definition
        relopening=None,  # TODO Check definition
        losscoeff=None,  # TODO Check definition
        bedfrictiontype=roughness_list[culvert.typeruwheid],
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
#%% #################################################################
# Set the crosssections for the branches

# hydamo.crosssections.convert.profiles(
#     crosssections=hydamo.profile,
#     crosssection_roughness=hydamo.profile_roughness,
#     profile_groups=hydamo.profile_group,
#     profile_lines=hydamo.profile_line,
#     param_profile=hydamo.param_profile,
#     param_profile_values=hydamo.param_profile_values,
#     branches=hydamo.branches,
#     roughness_variant="High",
# )

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
    # max_dist_to_struc=None,
    structures=structures,
)

models = Df2HydrolibModel(hydamo)


# %% ####################################################
# Export to DIMR configuration
fm.geometry.structurefile = [StructureModel(structure=models.structures)]
# fm.geometry.crosslocfile = CrossLocModel(crosssection=models.crosslocs)
# fm.geometry.crossdeffile = CrossDefModel(definition=models.crossdefs)

# fm.geometry.frictfile = []
# for i, fric_def in enumerate(models.friction_defs):
#    fric_model = FrictionModel(global_=fric_def)
#    fric_model.filepath = f"roughness_{i}.ini"
#    fm.geometry.frictfile.append(fric_model)

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
