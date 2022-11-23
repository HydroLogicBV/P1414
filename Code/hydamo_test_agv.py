import contextily as ctx
import matplotlib.pyplot as plt
from hydrolib.core.io.structure.models import *
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.dhydamo.io.dimrwriter import DIMRWriter
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh
from hydrolib.dhydamo.converters.df2hydrolibmodel import Df2HydrolibModel
from pathlib import Path


#%%
folder = r"D:\work\P1414"

hydamo = HyDAMO(extent_file= folder +"\GIS\WAGV\AGV_mask.shp")
hydamo.branches.read_shp(
    path = folder + r"\GIS\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp",
    column_mapping={"ruwheidsty": "typeruwheid"},
    index_col="code",
)

hydamo.bridges.read_shp(
    path = folder + r"\GIS\WAGV\brug_v13\brug_v13_clipped.shp",
    column_mapping={"naam": "globalid", "ruwheidsty":"typeruwheid", 
                    "intreeverl": "intreeverlies", "uittreever":"uittreeverlies"},
    index_col="code",
)

hydamo.culverts.read_shp(
    path = folder + r"\GIS\WAGV\duikersifonhevel_v13\duikersifonhevel_v13_clipped.shp",
    column_mapping = {"naam": "globalid","hoogtebinn": "hoogtebinnenonderkantbene",
                         "hoogtebin0":"hoogtebinnenonderkantbov", "vormkokeri":"vormkoker",
                         "hoogteopen":"hoogteopening", "breedteope":"breedteopening",
                         "intreeverl": "intreeverlies", "uittreever":"uittreeverlies",
                         "ruwheidsty":"typeruwheid" },
    index_col="code",
)

hydamo.weirs.read_shp(
    path = folder + r"\GIS\WAGV\stuw_v13\stuw_v13_clipped.shp",
    column_mapping = {"naam": "globalid", "afvoercoef": "afvoercoefficient", 
                      "soortstuwi":"soortstuw"},
    index_col = "code",
)


#%% Test plot the features

# fig = plt.figure()
# ax = plt.gca()
# hydamo.bridges.geometry.plot(ax=ax, label="Channel", linewidth=2, color="blue")
# ctx.add_basemap(ax, crs=28992, source=ctx.providers.OpenStreetMap.Mapnik)
# plt.show()

# print(hydamo)

# %% Check for circular features in branches

circular_branches = []
hydamo.branches_popped = hydamo.branches.copy()

for i, branch in enumerate(hydamo.branches.geometry):
    
    try:
        start = branch.coords[0]
        end = branch.coords[-1]
        if (start == end):
            code = hydamo.branches.code[i]
            hydamo.branches_popped = hydamo.branches_popped.drop(code)
    except IndexError: print('List index out of range for index',i)

hydamo.branches = hydamo.branches_popped     

#%% Snap structures to branches
hydamo.bridges.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=1100)
hydamo.bridges.dropna(axis=0, inplace=True, subset=["branch_offset"])

hydamo.culverts.snap_to_branch(hydamo.branches, snap_method='ends',maxdist=15)
hydamo.culverts.dropna(axis=0,inplace=True,subset=["branch_offset"])

hydamo.weirs.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=10)
hydamo.weirs.dropna(axis=0, inplace=True, subset=["branch_offset"])

#hydamo.structures.convert.culverts(hydamo.culverts, management_device=None)
#hydamo.structures.convert.bridges(hydamo.bridges)

# Adding the bridges in a loop does work, however, two fields are still missing:
# - csDefId: no clue what that does yet
# - shift: might be branch_offset or chainage as I defined it right now

for i, bridge in enumerate(hydamo.bridges.code):
    hydamo.structures.add_bridge(
        id = hydamo.bridges.code[i],
        name = hydamo.bridges.code[i],
        length = hydamo.bridges.lengte[i], 
        branchid = hydamo.bridges.branch_id[i],
        chainage = hydamo.bridges.branch_offset[i], # Unsure
        frictiontype='StricklerKs',
        csdefid = hydamo.bridges.code[i], # Unsure  
        shift = hydamo.bridges.branch_offset[i], # Unsure
        friction = hydamo.bridges.ruwheid[i],
        inletlosscoeff = hydamo.bridges.intreeverlies[i],
        outletlosscoeff = hydamo.bridges.uittreeverlies[i],
    )

# Collect all structures in a Structures dataframe
structures = hydamo.structures.as_dataframe(
    rweirs=False,
    bridges=True,
    uweirs=False,
    culverts=False,
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
    height=5.0, width=5.0, closed=False, roughnesstype="StricklerKs", roughnessvalue=30, name="default"
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
#fm.geometry.crosslocfile = CrossLocModel(crosssection=models.crosslocs)
#fm.geometry.crossdeffile = CrossDefModel(definition=models.crossdefs)

#fm.geometry.frictfile = []
#for i, fric_def in enumerate(models.friction_defs):
#    fric_model = FrictionModel(global_=fric_def)
#    fric_model.filepath = f"roughness_{i}.ini"
#    fm.geometry.frictfile.append(fric_model)

#fm.output.obsfile = [ObservationPointModel(observationpoint=models.obspoints)]

#extmodel = ExtModel()
#extmodel.boundary = models.boundaries_ext
#extmodel.lateral = models.laterals_ext
#fm.external_forcing.extforcefilenew = extmodel

#fm.geometry.inifieldfile = IniFieldModel(initial=models.inifields)
# for ifield, onedfield in enumerate(models.onedfieldmodels):
#     fm.geometry.inifieldfile.initial[ifield].datafile = OneDFieldModel(
#         global_= onedfield
#     )
#Now we write the file structure:

fm.filepath = Path(folder) / "fm" / "test.mdu"
dimr = DIMR()
dimr.component.append(
    FMComponent(name="DFM", workingDir=Path(folder) / "fm", model=fm, inputfile=fm.filepath)    
)
dimr.save(recurse=True)
#import shutil
#shutil.copy(data_path / "initialWaterDepth.ini", folder / "fm")

dimr = DIMRWriter(output_path=folder)
dimr.write_dimrconfig(fm)#, rr_model=drrmodel, rtc_model=drtcmodel)
dimr.write_runbat()
# %%
