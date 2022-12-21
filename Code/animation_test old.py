import os
import sys
from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from tqdm import tqdm

sys.path.append(r"D:\Work\git\GIS_tools\HydroLogic_Inundation_Toolbox\Readers")

from flowmeshreader import load_map_data, mesh_to_tiff
from plotting import raster_plot_with_context

# set paths
input_file_path = r"D:\Work\Project\P1414\Models\AGV\V5a\dflowfm\output\DFM_map.nc"
output_gif_path = r"D:\Work\Project\P1414\Models\AGV\V5a\dflowfm\output\waterdepth.gif"
output_tiff_path = r"D:\Work\Project\P1414\Models\AGV\V5a\dflowfm\output\tiff"
Path(output_tiff_path).mkdir(exist_ok=True)

# raster options
resolution = 30  # m
dhydro_resolution = 100
distance_tol = np.ceil(np.sqrt(2 * dhydro_resolution**2))  # m
interpolation = r"nearest"


# load mesh coordinates and data from netCDF
variable = r"Mesh2d_waterdepth"
map_data = load_map_data(input_file_path, variable)


def animate(ix, ax, blit=False):
    # ax.clear()

    output_file_path = output_tiff_path + r"\{}.tiff".format(ix)

    if Path(output_file_path).is_file():
        print("{}.tiff found".format(ix))
        fig, ax = raster_plot_with_context(
            raster_path=output_file_path,
            epsg=28992,
            clabel="water depth (m)",
            cmap="Reds",
            title="Water depth at last time step",
        )
    else:
        print("{}.tiff not found".format(ix))
        # convert to raster and save as tiff
        _, _, _ = mesh_to_tiff(
            map_data[ix, :],
            input_file_path,
            output_file_path,
            resolution,
            distance_tol,
            interpolation=interpolation,
        )

        fig, ax = raster_plot_with_context(
            raster_path=output_file_path,
            epsg=28992,
            clabel="water depth (m)",
            cmap="Reds",
            title="Water depth at last time step",
        )

    if blit:
        artists = ax.get_children()
        plt.close()

        return artists


def animate_0(ix=0):
    output_file_path = output_tiff_path + r"\{}.tiff".format(ix)

    if Path(output_file_path).is_file():
        print("{}.tiff found".format(ix))
        fig, ax = raster_plot_with_context(
            raster_path=output_file_path,
            epsg=28992,
            clabel="water depth (m)",
            cmap="Reds",
            title="Water depth at last time step",
        )

    else:
        print("{}.tiff not found".format(ix))
        # convert to raster and save as tiff
        _, _, _ = mesh_to_tiff(
            map_data[ix, :],
            input_file_path,
            output_file_path,
            resolution,
            distance_tol,
            interpolation=interpolation,
        )

        fig, ax = raster_plot_with_context(
            raster_path=output_file_path,
            epsg=28992,
            clabel="water depth (m)",
            cmap="Reds",
            title="Water depth at last time step",
        )

    return fig, ax


fig, ax = animate_0(0)

bblit = True
print(ax.get_children())
ani = FuncAnimation(fig, partial(animate, ax=ax, blit=bblit), frames=5, interval=200, blit=bblit)
ani.save(output_gif_path, dpi=300, writer=PillowWriter(fps=5))
