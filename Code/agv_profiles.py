import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import LineString
from tqdm import tqdm

profile_points_path = (
    r"D:\Work\Project\P1414\GIS\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped.shp"
)
profiles_output_path = r"D:\Work\Project\P1414\GIS\WAGV\profielen.gpkg"

# Load shapefile with profile points & group them by metingprof attribute
profile_points = gpd.read_file(profile_points_path)
grouped_points = profile_points.groupby(by="metingprof")

# initialize empty lists
profile_groups = []
profile_lines = []
profile_points = []
roughnes_profiles = []

# Loop over the grouped profile points to create a line from the points
count = 0
for name, group in tqdm(grouped_points):
    profile_line_id = r"AGV_" + str(name)
    list_of_points = []

    _profile_group = dict(
        [("globalid", profile_line_id), ("brugid", None), ("stuwid", None), ("geometry", None)]
    )
    profile_groups.append(_profile_group)

    # sort points based on codevolg number
    sorted_group = group.sort_values("codevolgnu")

    # if both bottom and sludge measurements are presten, select the first
    if sorted_group["typebodemi"].isin([1]).any() and sorted_group["typebodemi"].isin([2]).any():
        sorted_group.drop(sorted_group[sorted_group["typebodemi"] == 2].index, inplace=True)

    # add points to line
    for ix, row in sorted_group.iterrows():
        # loop to turn multipoint into point. Probably redundant, but safer
        for point in row.geometry.geoms:
            # append point to list for line generation
            list_of_points.append(point)

            # create entry for point shape and append
            point_id = "AGV_" + str(row["code"])
            _point = dict(
                [
                    ("code", row["code"]),
                    ("globalid", point_id),
                    ("profiellijnid", profile_line_id),
                    ("codevolgnummer", row["codevolgnu"]),
                    ("geometry", point),
                ]
            )
            profile_points.append(_point)

            # create roughness table entry
            _roughness_profile = dict(
                [
                    ("profielpuntid", point_id),
                    ("typeruwheid", row["ruwheidsty"]),
                    ("ruwheidhoog", row["ruwheidswa"]),
                    ("ruwheidlaag", row["ruwheidsw0"]),
                    ("geometry", None),
                ]
            )
            roughnes_profiles.append(_roughness_profile)

    if len(list_of_points) > 1:
        # add line to list
        profile_line = LineString(list_of_points)
        _profile_line = dict(
            [("globalid", profile_line_id), ("profielgroepid", name), ("geometry", profile_line)]
        )
        profile_lines.append(_profile_line)

    count += 1
    # if count >= 10:
    #     break

# Create GeoDataFrame from list of dicts
profile_groups = gpd.GeoDataFrame(profile_groups, geometry="geometry", crs=28992)
profile_lines = gpd.GeoDataFrame(profile_lines, geometry="geometry", crs=28992)
profile_points = gpd.GeoDataFrame(profile_points, geometry="geometry", crs=28992)
roughnes_profiles = gpd.GeoDataFrame(roughnes_profiles, geometry="geometry", crs=28992)

# print(profile_groups)
# print(roughnes_profiles)

layers = dict(
    [
        ("profielgroep", profile_groups),
        ("profiellijn", profile_lines),
        ("profielpunt", profile_points),
        ("ruwheidsprofiel", roughnes_profiles),
    ]
)

for name, layer in layers.items():
    layer.to_file(filename=profiles_output_path, driver="GPKG", layer=name)

# profile_points.plot()
# plt.show()
