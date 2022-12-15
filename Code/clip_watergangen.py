# %% Init
print("Initializing")
import sys

import geopandas as gpd

BUFFER_DIST = 50
EPSG = 28992


def clip_branches(
    in_branches_path: str,
    buffered_branches: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    in_branches = gpd.read_file(in_branches_path).to_crs(buffered_branches.crs)
    out_branches = gpd.overlay(
        in_branches, buffered_branches, how="intersection", keep_geom_type=True
    )
    out_branches = out_branches.explode()
    return out_branches


def read_rm_branches(rm_branches_path: str) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    rm_branches = gpd.read_file(rm_branches_path).to_crs(crs=EPSG)
    ondd_bool = (rm_branches["Source"].str.contains("OND")) | (
        rm_branches["Target"].str.contains("OND")
    )
    onderdoorgangen = rm_branches.loc[ondd_bool, :]
    rm_branches = rm_branches.loc[~ondd_bool, :]
    return rm_branches, onderdoorgangen


old_rm_branches_path = r"D:\Work\Project\P1414\GIS\Randstadmodel_oud\rm_Branches_28992.shp"

old_rm_branches, onderdoorgangen = read_rm_branches(old_rm_branches_path)
buffered_old_rm_branches = gpd.GeoDataFrame(
    geometry=old_rm_branches.buffer(distance=BUFFER_DIST)
).dissolve(by=None)

onderdoorgangen_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\onderdoorgangen.shp"
onderdoorgangen.to_file(onderdoorgangen_path)

# %% AGV
print("AGV")

agv_branches_path = r"D:\Work\Project\P1414\GIS\WAGV\hydroobject_v13\hydroobject_v13_clipped.shp"
clipped_agv_branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\AGV.shp"

intersected_branches = clip_branches(
    in_branches_path=agv_branches_path, buffered_branches=buffered_old_rm_branches
)
intersected_branches.to_file(clipped_agv_branches_path)
# %% HHD
print("HHD")

HHD_branches_path = r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water.shp"
clipped_HHD_branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHD.shp"

intersected_branches = clip_branches(
    in_branches_path=HHD_branches_path, buffered_branches=buffered_old_rm_branches
)
intersected_branches.to_file(clipped_HHD_branches_path)

# %% HDSR
print("HDSR")

HDSR_branches_path = r"D:\Work\Project\P1414\GIS\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
clipped_HDSR_branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HDSR.shp"

intersected_branches = clip_branches(
    in_branches_path=HDSR_branches_path, buffered_branches=buffered_old_rm_branches
)
intersected_branches.to_file(clipped_HDSR_branches_path)

# %% HHR
print("HHR")

HHR_branches_path = r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Watergang\Watergang_as.shp"
clipped_HHR_branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHR.shp"

intersected_branches = clip_branches(
    in_branches_path=HHR_branches_path, buffered_branches=buffered_old_rm_branches
)
intersected_branches.to_file(clipped_HHR_branches_path)

# %% HHSK
print("HHSK")

HHSK_branches_path = r"D:\Work\Project\P1414\GIS\HHSK\Legger\Hoofdwatergang.shp"
clipped_HHSK_branches_path = r"D:\Work\Project\P1414\GIS\Uitgesneden watergangen\HHSK.shp"

intersected_branches = clip_branches(
    in_branches_path=HHSK_branches_path, buffered_branches=buffered_old_rm_branches
)
intersected_branches.to_file(clipped_HHSK_branches_path)
