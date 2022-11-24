from pathlib import Path

import contextily as ctx
import matplotlib.pyplot as plt
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.structure.models import *
from hydrolib.dhydamo.converters.df2hydrolibmodel import Df2HydrolibModel
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh
from hydrolib.dhydamo.io.dimrwriter import DIMRWriter
from tqdm import tqdm

#%%
folder = r"D:\work\P1414"
#folder = r"D:\Work\Project\P1414"

hydamo = HyDAMO(extent_file=folder + "\GIS\WAGV\AGV_mask.shp")
hydamo.branches.read_shp(

    path = folder + r"\GIS\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp",
    column_mapping ={
        "ruwheidsty":"typeruwheid",
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
        "naam": "globalid",
        "afvoercoef": "afvoercoefficient",
        "soortstuwi": "soortstuw",
    },
    index_col="code",
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

#hydamo.profile_roughness.set_data
#%% Test plot the features

fig = plt.figure()
ax = plt.gca()
hydamo.branches.geometry.plot(ax=ax, label='Channels',linewidth=0.5,color='black')
hydamo.culverts.geometry.plot(ax=ax, label="Culverts", linewidth=2, color="blue")
hydamo.bridges.geometry.plot(ax=ax, label="Bridges", markersize=2, color="red")
hydamo.weirs.geometry.plot(ax=ax, label="Weirs", markersize=2, color="orange")
ctx.add_basemap(ax, crs=28992, source=ctx.providers.OpenStreetMap.Mapnik)
plt.show()

#print(hydamo)

# %% Check for circular features in branches

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
hydamo.branches_popped = hydamo.branches_popped.drop([None], axis=0)

hydamo.branches = hydamo.branches_popped

#%% Snap structures to branches
hydamo.bridges.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=10)
hydamo.bridges.dropna(axis=0, inplace=True, subset=["branch_offset"])

hydamo.culverts.snap_to_branch(hydamo.branches, snap_method="ends", maxdist=15)
hydamo.culverts.dropna(axis=0, inplace=True, subset=["branch_offset"]) # Without branch
hydamo.culverts.dropna(axis=0,inplace=True,subset=["hoogteopening"]) # Without dimensions
hydamo.culverts.drop(hydamo.culverts[hydamo.culverts['vormkoker'] != (1 or 3)].index, inplace=True) # Without unknown shapes

# hydamo.weirs.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=10)
# hydamo.weirs.dropna(axis=0, inplace=True, subset=["branch_offset"])

# hydamo.structures.convert.weirs(
#     weirs=hydamo.weirs,
# )

# hydamo.structures.convert.culverts(hydamo.culverts, management_device=None)
# hydamo.structures.convert.bridges(hydamo.bridges)

# Add structures in a loop to provide the right format 
def get_crosssection_culvert_AGV(shape: int = 1, height: float = None, 
                                width: float = None, closed:int = 1):
    shapedict = {1:"circle",3: "rectangle"} 
    shape_str = shapedict[shape]

    # Include the diameter when the culvert is a circle
    if shape == 1: diameter = width
    else: diameter = None

    crosssection = {"shape"    : shape_str,
                    "diameter" : diameter,
                    "height"   : height,
                    "width"    : width,
                    "closed"   : closed,
                    }
    return crosssection


# Define a list of possible roughness types to interpret the underlying data 
roughness_list=["Bos en Bijkerk", "Chezy", "Manning", 
                "StricklerKn", "StricklerKs", "White Colebrook",
    ]


for i, bridge in enumerate(tqdm(hydamo.bridges.code)):
    hydamo.structures.add_bridge(
        id = hydamo.bridges.code[i],
        name = hydamo.bridges.code[i],
        length = hydamo.bridges.lengte[i],
        branchid = hydamo.bridges.branch_id[i],
        chainage = hydamo.bridges.branch_offset[i],  #TODO Validate the use of offset
        frictiontype = roughness_list[hydamo.bridges.typeruwheid[i]],
        csdefid = hydamo.bridges.code[i],  #TODO Check influence of csdefid
        shift = hydamo.bridges.branch_offset[i],  #TODO Validate the use of offset
        friction = hydamo.bridges.ruwheid[i],
        inletlosscoeff = hydamo.bridges.intreeverlies[i],
        outletlosscoeff = hydamo.bridges.uittreeverlies[i],
    )

for i, culvert in enumerate(tqdm(hydamo.culverts.iterrows())):
    hydamo.structures.add_culvert(
        id = hydamo.culverts.code[i],
        name = hydamo.culverts.code[i],
        branchid = hydamo.culverts.branch_id[i],
        chainage = hydamo.culverts.branch_offset[i], #TODO Valdiate the use of offset 
        leftlevel = hydamo.culverts.hoogtebinnenonderkantbov[i],
        rightlevel = hydamo.culverts.hoogtebinnenonderkantbene[i],
        length = hydamo.culverts.lengte[i],
        inletlosscoeff = hydamo.culverts.intreeverlies[i],
        outletlosscoeff = hydamo.culverts.uittreeverlies[i],
        crosssection = get_crosssection_culvert_AGV(shape=hydamo.culverts.vormkoker[i],
                                                    height=hydamo.culverts.hoogteopening[i],
                                                    width=hydamo.culverts.breedteopening[i],
                                                    closed=1), #TODO Check for automatic cross section based on shape
        numlosscoeff = None , #TODO Check definition
        relopening = None, #TODO Check definition
        losscoeff = None, #TODO Check definition
        bedfrictiontype = roughness_list[hydamo.culverts.typeruwheid[i]],
        bedfriction = hydamo.culverts.ruwheid[i],
    )


# Add observation points (empty for now)
hydamo.observationpoints.add_points(
    [],
    [],
    locationTypes=["1d","1d"],
    snap_distance=10.0
)

# Collect all structures in a Structures dataframe
structures = hydamo.structures.as_dataframe(
    rweirs=False,
    bridges=True,
    uweirs=False,
    culverts=True,
    orifices=False,
    pumps=False,
)
#%% Set the profiles for the branches
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
#%% Create the FM model
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


# %% Export to DIMR configuration
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
