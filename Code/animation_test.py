import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from tqdm import tqdm

sys.path.append(r"D:\Work\git\GIS_tools\HydroLogic_Inundation_Toolbox\Readers")

from flowmeshreader import load_map_data, mesh_to_tiff
from plotting import raster_plot_with_context

# set paths
input_file_path = r"D:\Work\Project\P1414\Models\AGV\V5a\dflowfm\output\DFM_map.nc"
output_gif_path = r"D:\Work\Project\P1414\Models\AGV\V5a\dflowfm\output\waterdepth.webp"
output_fig_path = r"D:\Work\Project\P1414\Models\AGV\V5a\dflowfm\output\fig"
Path(output_fig_path).mkdir(exist_ok=True)

# raster options
resolution = 30  # m
dhydro_resolution = 100
distance_tol = np.ceil(np.sqrt(2 * dhydro_resolution**2))  # m
interpolation = r"nearest"
vmin = 1
vmax = 1.6


# load mesh coordinates and data from netCDF
variable = r"Mesh2d_waterdepth"
map_data = load_map_data(input_file_path, variable)
map_data[map_data < 0.01] = np.nan

frames = []
for ix in tqdm(range(25)):
    output_png_file_path = output_fig_path + r"\{}.png".format(ix)
    output_tiff_file_path = output_fig_path + r"\{}.tiff".format(ix)

    if Path(output_png_file_path).is_file():
        tqdm.write("{}.png found".format(ix))

    elif Path(output_tiff_file_path).is_file():
        tqdm.write("{}.tiff found".format(ix))
        fig, ax = raster_plot_with_context(
            raster_path=output_tiff_file_path,
            epsg=28992,
            clabel="water depth (m)",
            cmap="Reds",
            title="Water depth",
            vmin=vmin,
            vmax=vmax,
        )
        fig.savefig(output_png_file_path, dpi=300)
        plt.close(fig)

    else:
        tqdm.write("{}.tiff not found".format(ix))
        # convert to raster and save as tiff
        _, _, _ = mesh_to_tiff(
            map_data[ix, :],
            input_file_path,
            output_tiff_file_path,
            resolution,
            distance_tol,
            interpolation=interpolation,
        )

        fig, ax = raster_plot_with_context(
            raster_path=output_tiff_file_path,
            epsg=28992,
            clabel="water depth (m)",
            cmap="Reds",
            title="Water depth",
            vmin=vmin,
            vmax=vmax,
        )
        fig.savefig(output_png_file_path, dpi=300)
        plt.close(fig)

    new_frame = Image.open(output_png_file_path)
    frames.append(new_frame)

frames[0].save(
    output_gif_path,
    format="webp",
    append_images=frames[1:],
    disposal=2,
    duration=500,
    loop=1,
    optimize=True,
    save_all=True,
    transparency=1,
)
