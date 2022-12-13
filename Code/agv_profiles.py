import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import LineString, MultiPoint, Point
from tqdm import tqdm

dist_tol = 0.25

roughness_mapping = {
    "Chezy": "Chezy",
    "Manning": "Manning",
    "StricklerKn": "StricklerNikuradse",
    "StricklerKs": "Strickler",
    "White Colebrook": "WhiteColebrook",
    "Bos en Bijkerk": "deBosBijkerk",
    "Onbekend": "Strickler",
    "Overig": "Strickler",
}
roughness_mapping = list(roughness_mapping)

profile_points_path = (
    r"D:\Work\Project\P1414\GIS\WAGV\metingprofielpunt_v13\metingprofielpunt_v13_clipped_rm.shp"
)
profiles_output_path = r"D:\Work\Project\P1414\GIS\WAGV\profielen_light.gpkg"

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
    # sort points based on codevolg number
    sorted_group = group.sort_values("codevolgnu")

    # if both bottom and sludge measurements are present, select the first
    if sorted_group["typebodemi"].isin([1]).any() and sorted_group["typebodemi"].isin([2]).any():
        sorted_group.drop(sorted_group[sorted_group["typebodemi"] == 2].index, inplace=True)

    # skip 'line' if it only has one point
    if sorted_group.shape[0] < 5:
        continue

    # Initialize name and list of points
    profile_line_id = r"AGV_" + str(name)
    list_of_points = []

    # Add info to profile_group table
    profile_group_id = profile_line_id + "groep"
    _profile_group = dict(
        [
            ("globalid", profile_group_id),
            ("brugid", ""),
            ("stuwid", ""),
            ("geometry", None),
        ]
    )
    # _profile_group = dict([("globalid", profile_line_id), ("geometry", None)])
    profile_groups.append(_profile_group)

    code_volg_nr = 0
    # add points to line
    for ix, row in sorted_group.iterrows():
        if type(row.geometry) == Point:
            l_points = [row.geometry]
        elif type(row.geometry) == MultiPoint:
            l_points = row.geometry.geoms

        for point in l_points:
            # Check if profile points are too close together, and skip if that is the case
            # if ix < sorted_group.shape[0]:
            if code_volg_nr > 0:
                p_0 = list_of_points[code_volg_nr - 1]
                p_1 = point
                p_dist = p_0.distance(p_1)

                if p_dist < dist_tol:
                    continue

            # Because sorted_group has been sorted on codevolgnr, we can assume sequentiallity
            code_volg_nr += 1

            # append point to list for line generation
            list_of_points.append(point)

            # create entry for point shape and append
            # point_id = "AGV_" + str(row["code"])
            point_id = "AGV_" + str(name) + r"_" + str(code_volg_nr)
            _point = dict(
                [
                    ("code", row["code"]),
                    ("globalid", point_id),
                    ("profiellijnid", profile_line_id),
                    # ("codevolgnummer", row["codevolgnu"]),
                    ("codevolgnummer", code_volg_nr),
                    ("geometry", point),
                ]
            )
            profile_points.append(_point)

            # create roughness table entry
            _roughness_profile = dict(
                [
                    ("code", row["code"]),
                    ("profielpuntid", point_id),
                    ("typeruwheid", roughness_mapping[int(row["ruwheidsty"]) - 1]),
                    ("ruwheidhoog", float(row["ruwheidswa"])),
                    ("ruwheidlaag", float(row["ruwheidsw0"])),
                    ("geometry", None),
                ]
            )
            roughnes_profiles.append(_roughness_profile)

    # # Check for line or point
    # if len(list_of_points) > 1:

    #     profile_line = LineString(list_of_points)
    # else:
    #     profile_line = point

    # Convert points to line
    profile_line = LineString(list_of_points)

    # add line to list
    _profile_line = dict(
        [
            ("globalid", profile_line_id),
            ("profielgroepid", profile_group_id),
            ("geometry", profile_line),
        ]
    )
    profile_lines.append(_profile_line)

    # break if required
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
