import io
import os.path
import zipfile
from pathlib import Path
from zipfile import BadZipFile

import geopandas as gpd
import rasterio
import requests
from rasterio.merge import merge
from tqdm import tqdm

# set parameters
ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4"
ahn_merged_path = ahn_path + r"\ahn4_merged.TIF"
ahn_tiles_path = r"D:\Work\Project\P1414\GIS\AHN\AHN_bladen.shp"
ahn_url = r"https://ns_hwh.fundaments.nl/hwh-ahn/ahn4/02b_DTM_5m/"

branches_path = r"D:\Work\Project\P1414\GIS\Randstadmodel_oud\rm_Branches_28992_edited.shp"

max_download_attempts = 2

# create path
Path(ahn_path).mkdir(exist_ok=True)

# determine extend of AHN to download
branches = gpd.read_file(branches_path)
convex_hull = branches.dissolve(by=None).envelope

# Load ahn tiles and clip with conex_hull of branches
ahn_tiles = gpd.read_file(ahn_tiles_path)
clipped_ahn_tiles = ahn_tiles.clip(convex_hull)

# open the final file
# dst = rasterio.open(ahn_merged_path, 'w')

# download all clipped ahn tiles
pbar = tqdm(clipped_ahn_tiles.iterrows(), total=clipped_ahn_tiles.shape[0])
path_list = []
for ix, tile in pbar:
    bladnr = tile["bladnr"]
    file_path = ahn_path + r"\M5_{}.TIF".format(str(bladnr).upper())

    # check if file already exists locally
    if os.path.isfile(file_path):
        pbar.set_description("found AHN tile {}".format(bladnr))

    # else download
    else:
        pbar.set_description("downloading AHN tile {}".format(bladnr))
        tile_url = ahn_url + r"M5_{}.zip".format(str(bladnr).upper())

        try:
            # try to download as many times as max_download_attempts
            for jx in range(max_download_attempts):
                r = requests.get(tile_url)
                if r.ok:
                    break

            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(ahn_path)

        # if file does not exis, skip
        except BadZipFile:
            pbar.write("{} not found".format(bladnr))
            continue

    path_list.append(file_path)

ahn_merged, out_transform = merge(path_list)
count, height, width = ahn_merged.shape
profile = {
    "count": count,
    "height": height,
    "width": width,
    "transform": out_transform,
    "dtype": "float32",
}

with rasterio.open(ahn_merged_path, "w", **profile) as dst:
    dst.write(ahn_merged)
