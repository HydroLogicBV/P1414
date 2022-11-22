import contextily as ctx
import matplotlib.pyplot as plt
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.dhydamo.io.dimrwriter import DIMRWriter
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh

from pathlib import Path
#%%
folder = r"D:\work\P1414"

hydamo = HyDAMO(extent_file= folder +"\GIS\WAGV\AGV_mask.shp")
hydamo.branches.read_shp(
    path = folder +"\GIS\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp",
    column_mapping={"ruwheidsty": "typeruwheid"},
    index_col="code",
)

#%% Test plot the features

# fig = plt.figure()
# ax = plt.gca()
# hydamo.branches.geometry.plot(ax=ax, label="Channel", linewidth=2, color="blue")
# ctx.add_basemap(ax, crs=28992, source=ctx.providers.OpenStreetMap.Mapnik)
# plt.show()

# print(hydamo)

# %% Check for circular features in branches

circular_branches = []
hydamo_branches_popped = hydamo.branches.copy()

for i, branch in enumerate(hydamo.branches.geometry):
    
    try:
        start = branch.coords[0]
        end = branch.coords[-1]
        if (start == end):
            code = hydamo.branches.code[i]
            hydamo_branches_popped = hydamo_branches_popped.drop(code)
    except IndexError: print('List index out of range for index',i)
    
print(len(hydamo.branches))
print(len(hydamo_branches_popped))   


#%% Create the FM model
fm = FMModel()
# Set start and stop time
fm.time.refdate = 20160601
fm.time.tstop = 2 * 3600 * 24

mesh.mesh1d_add_branches_from_gdf(
    fm.geometry.netfile.network,
    branches=hydamo_branches_popped,
    branch_name_col="code",
    node_distance=20,
    # max_dist_to_struc=None,
    # structures=structures,
)


# %% Export to DIMR configuration
#fm.geometry.structurefile = [StructureModel(structure=models.structures)]
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