import contextily as ctx
import matplotlib.pyplot as plt
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh

#%%
hydamo = HyDAMO(extent_file=r"D:\Work\Project\P1414\GIS\WAGV\AGV_mask.shp")
hydamo.branches.read_shp(
    path="D:\Work\Project\P1414\GIS\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp",
    column_mapping={"ruwheidsty": "typeruwheid"},
    index_col="naam",
)

#%%

# fig = plt.figure()
# ax = plt.gca()
# hydamo.branches.geometry.plot(ax=ax, label="Channel", linewidth=2, color="blue")
# ctx.add_basemap(ax, crs=28992, source=ctx.providers.OpenStreetMap.Mapnik)
# plt.show()

# print(hydamo)


#%%
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
    # structures=structures,
)

# %%
# Check for circular features in branches
circular_branches = []
for i, branch in enumerate(hydamo.branches.geometry):
    
    try:
        start = branch.coords[0]
        end = branch.coords[-1]
        if (start == end):
            circular_branches.append([branch,i])
    except IndexError: print('List index out of range for index',i)
    
for circle in circular_branches:
    print('ID circles: ',circle[1])   
