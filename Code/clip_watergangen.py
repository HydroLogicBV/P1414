# %% Init
print("Initializing")
import uuid
from copy import copy

import geopandas as gpd
import numpy as np
from scipy.spatial import KDTree

BUFFER_DIST = 10
EPSG = 28992

AGV = True
HHD = True
HDSR = True
HHR = True
HHSK = True


def clip_branches(
    in_branches_path: str,
    buffered_branches: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:

    print("Intersecting buffered branches with input branches")

    # Load input branches, explode multilinestrings, and assign unique ids
    in_branches = (
        gpd.read_file(in_branches_path).to_crs(buffered_branches.crs).explode(index_parts=False)
    )
    in_branches["new_id"] = [str(uuid.uuid4()) for _ in range(in_branches.shape[0])]

    # intersect branches with buffered branches from old model
    out_branches = gpd.overlay(
        in_branches, buffered_branches, how="intersection", keep_geom_type=True
    )

    print("Checking for split branches and fixing them")

    # setting index to new_id for easier indexing later
    in_branches.set_index("new_id", inplace=True)
    out_branches.set_index("new_id", inplace=True)

    # find branches that are MultiLineString due to partial intersection
    MLS_branches_bool = out_branches.geometry.type == "MultiLineString"
    MLS_branches = out_branches.loc[MLS_branches_bool]

    print("Number of split branches: {}".format(np.sum(MLS_branches_bool)))

    # replace branches that are MultiLineString with original branch data
    for ix, branch in MLS_branches.iterrows():
        out_branches.loc[branch.name] = in_branches.loc[branch.name]

    n_split_branches = np.sum(out_branches.geometry.type == "MultiLineString")
    print("Remaining number of split branches: {}".format(n_split_branches))

    if n_split_branches > 0:
        print(out_branches.loc[out_branches.geometry.type == "MultiLineString"])
        raise ValueError("still found split branches...")

    # return out_branches
    return skip_branches_dist(in_branches=out_branches)


def skip_branches_dist(
    in_branches: gpd.GeoDataFrame, max_distance: float = 0.5, min_connectivity: int = 2
) -> gpd.GeoDataFrame:
    """
    Function that removes branches from GeoDataFrame if they are not sufficiently connected to other branches
    The criterion for this is min_connectivity for the number of branches that should coincide in a point
    Max_distance is the distance tolerance for branches starting and ending in a point

    Args:
        in_branches (gpd.GeoDataFrame): GeoDataFrame with branches containging a geometry column with LineStrings
        max_distance (float): maximum distance at which branches are considered connected
        min_connectivity (int): minumum number of branches that should meet in a point

    Returns:
        out_branhces (gpd.GeoDataFrame): GeoDataFrame with branches that fullfill the connectivity criterion
    """

    print("Checking for isolated branches")
    n_in = in_branches.shape[0]

    # skip branches further than max_dist away from another branch
    point_list = []
    for ix, branch in in_branches.iterrows():
        coords = branch.geometry.coords
        point_list.append(coords[0])
        point_list.append(coords[-1])

    # build spatial KDTree from start and end points
    kdtree = KDTree(point_list)
    sdm = kdtree.sparse_distance_matrix(
        kdtree, max_distance=max_distance, output_type="coo_matrix"
    )

    # matrix is symetrical, so picking one of the two axis to count number of non empty entries (explicit zeros allowed)
    nnz = sdm.getnnz(axis=0)

    # compare number of points within max_distance to min_connectivity
    bool_array = nnz < min_connectivity

    # Drop branches with insufficient connectivity on both sides
    out_branches = copy(in_branches)
    for ix, (name, branch) in enumerate(in_branches.iterrows()):
        if (bool_array[ix * 2]) & (bool_array[ix * 2 + 1]):
            out_branches.drop(index=name, inplace=True)
        # elif bool_array[ix * 2 + 1]:
        #     out_branches.drop(index=name, inplace=True)

    # out_branches = in_branches
    n_out = out_branches.shape[0]
    print("Dropped {} isolated branches".format(n_in - n_out))
    return out_branches


def read_rm_branches(rm_branches_path: str) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    rm_branches = gpd.read_file(rm_branches_path).to_crs(crs=EPSG)
    ondd_bool = (rm_branches["Source"].str.contains("OND")) | (
        rm_branches["Target"].str.contains("OND")
    )
    onderdoorgangen = rm_branches.loc[ondd_bool, :]
    rm_branches = rm_branches.loc[~ondd_bool, :]
    return rm_branches, onderdoorgangen


p_folder = r"D:\Work\Project\P1414\GIS"
old_rm_branches_path = p_folder + r"\Randstadmodel_oud\rm_Branches_28992.shp"
buffered_branches_path = p_folder + r"\Randstadmodel_oud\buffered_branches.shp"

old_rm_branches, onderdoorgangen = read_rm_branches(old_rm_branches_path)
buffered_old_rm_branches = gpd.GeoDataFrame(
    geometry=old_rm_branches.buffer(distance=BUFFER_DIST)
).dissolve(by=None)
buffered_old_rm_branches.to_file(buffered_branches_path)

onderdoorgangen_path = p_folder + r"\Randstadmodel_oud\onderdoorgangen.shp"
onderdoorgangen.to_file(onderdoorgangen_path)

# %% AGV
if AGV:
    print("AGV")

    agv_branches_path = p_folder + r"\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp"
    clipped_agv_branches_path = p_folder + r"\Uitgesneden watergangen\AGV_v3_test.shp"

    intersected_branches = clip_branches(
        in_branches_path=agv_branches_path, buffered_branches=buffered_old_rm_branches
    )
    intersected_branches.to_file(clipped_agv_branches_path)
# %% HHD
if HHD:
    print("HHD")

    HHD_branches_path = (
        p_folder + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water.shp"
    )
    clipped_HHD_branches_path = p_folder + r"\Uitgesneden watergangen\HHD_v3_test.shp"

    intersected_branches = clip_branches(
        in_branches_path=HHD_branches_path, buffered_branches=buffered_old_rm_branches
    )
    intersected_branches.to_file(clipped_HHD_branches_path)

# %% HDSR
if HDSR:
    print("HDSR")

    HDSR_branches_path = p_folder + r"\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
    clipped_HDSR_branches_path = p_folder + r"\Uitgesneden watergangen\HDSR_v3_test.shp"

    intersected_branches = clip_branches(
        in_branches_path=HDSR_branches_path, buffered_branches=buffered_old_rm_branches
    )
    intersected_branches.to_file(clipped_HDSR_branches_path)

# %% HHR
if HHR:
    print("HHR")

    HHR_branches_path = p_folder + r"\HHRijnland\Legger\Watergang\Watergang_as.shp"
    clipped_HHR_branches_path = p_folder + r"\Uitgesneden watergangen\HHR_v3_test.shp"

    intersected_branches = clip_branches(
        in_branches_path=HHR_branches_path, buffered_branches=buffered_old_rm_branches
    )
    intersected_branches.to_file(clipped_HHR_branches_path)

# %% HHSK
if HHSK:
    print("HHSK")

    HHSK_branches_path = p_folder + r"\HHSK\Legger\Hoofdwatergang.shp"
    clipped_HHSK_branches_path = p_folder + r"\Uitgesneden watergangen\HHSK_v3_test.shp"

    intersected_branches = clip_branches(
        in_branches_path=HHSK_branches_path, buffered_branches=buffered_old_rm_branches
    )
    intersected_branches.to_file(clipped_HHSK_branches_path)
