#!/usr/bin/env python3
"""Demonstrate PCAT's deterministic Gaussian-PSF image forward model."""

from tdpy.verbosity import print

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tdpy.cli import add_plot_arguments

from pcat.demo import configure_example_style, save_example_figure
from pcat.image import forward_model_image, make_gaussian_scene

EXAMPLE_PATH = Path(__file__).resolve().parent
IMAGE_SIZE_PIXELS = 51
PSF_SIGMA_PIXELS = 1.5
BACKGROUND_RATE = 5.0  # [electron s^-1 pixel^-1]


def create_simulated_scene():
    """Return an explicitly assumed two-source scene [electron s^-1 pixel^-1]."""
    return make_gaussian_scene(
        (IMAGE_SIZE_PIXELS, IMAGE_SIZE_PIXELS),
        ((24.0, 25.0, 1000.0), (28.0, 25.0, 650.0)),
        intrinsic_sigma_pixels=0.35,
    )


def run_example(typefileplot="png", typeplotback="white"):
    """Forward-model the simulated scene and save its component images."""
    foreground = configure_example_style(typeplotback)
    source_image = create_simulated_scene()
    products = forward_model_image(
        source_image,
        psf_sigma_pixels=PSF_SIGMA_PIXELS,
        background=BACKGROUND_RATE,
    )

    figure, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), facecolor="white")
    scene_slice = slice(17, 35)
    scene_extent = (16.5, 34.5, 16.5, 34.5)
    kernel_half_width = products["psf_kernel"].shape[0] / 2.0
    kernel_extent = (-kernel_half_width, kernel_half_width, -kernel_half_width, kernel_half_width)
    panels = (
        (products["source_image"][scene_slice, scene_slice], "Intrinsic source scene",
         "Flux rate [electron s$^{-1}$ pixel$^{-1}$]", scene_extent, "Detector"),
        (products["psf_kernel"], "Normalized Gaussian PSF", "Kernel weight", kernel_extent, "PSF"),
        (products["observed_image"][scene_slice, scene_slice], "Predicted detector image",
         "Flux rate [electron s$^{-1}$ pixel$^{-1}$]", scene_extent, "Detector"),
    )
    for axis, (image, title, colorbar_label, extent, coordinate_name) in zip(axes, panels):
        artist = axis.imshow(image, origin="lower", interpolation="nearest", cmap="magma", extent=extent)
        figure.colorbar(artist, ax=axis, label=colorbar_label)
        axis.set_xlabel(f"{coordinate_name} x [pixel]", color=foreground)
        axis.set_ylabel(f"{coordinate_name} y [pixel]", color=foreground)
        axis.set_title(title)
        axis.grid(False)
    figure.suptitle(f"PCAT image formation, PSF sigma {PSF_SIGMA_PIXELS:.1f} pixels")
    figure.tight_layout()
    path = save_example_figure(
        figure, EXAMPLE_PATH, "pcat_image_forward_model", typefileplot
    )
    return products, path


def run_bokeh_example(products, output_path):
    """Write the scene, PSF, and detector image as linked interactive panels."""
    from bokeh.layouts import row
    from bokeh.models import ColorBar, LinearColorMapper
    from bokeh.palettes import Magma256
    from bokeh.plotting import figure

    from tdpy.plotting import save_bokeh_html

    scene_slice = slice(17, 35)
    kernel_half_width = products["psf_kernel"].shape[0] / 2.0
    panels = (
        (products["source_image"][scene_slice, scene_slice], "Intrinsic source scene",
         (16.5, 34.5, 16.5, 34.5)),
        (products["psf_kernel"], "Normalized Gaussian PSF",
         (-kernel_half_width, kernel_half_width, -kernel_half_width, kernel_half_width)),
        (products["observed_image"][scene_slice, scene_slice], "Predicted detector image",
         (16.5, 34.5, 16.5, 34.5)),
    )
    plots = []
    for image, title, extent in panels:
        low = 0.0
        high = max(float(np.max(image)), np.finfo(float).eps)
        color_mapper = LinearColorMapper(palette=Magma256, low=low, high=high)
        plot = figure(
            title=title,
            width=360,
            height=340,
            x_range=(extent[0], extent[1]),
            y_range=(extent[2], extent[3]),
            match_aspect=True,
            tools="pan,wheel_zoom,box_zoom,reset,save",
        )
        plot.image(
            image=[image],
            x=extent[0],
            y=extent[2],
            dw=extent[1] - extent[0],
            dh=extent[3] - extent[2],
            color_mapper=color_mapper,
        )
        plot.add_layout(ColorBar(color_mapper=color_mapper, title="Image value"), "right")
        plot.xaxis.axis_label = "x [pixel]"
        plot.yaxis.axis_label = "y [pixel]"
        plot.grid.visible = False
        plots.append(plot)
    return save_bokeh_html(
        row(*plots, sizing_mode="stretch_width"),
        output_path,
        title="PCAT deterministic image formation",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_plot_arguments(parser)
    parser.add_argument("--typeplotback", choices=("white", "dark"), default="white")
    parser.add_argument("--bokeh", action="store_true", help="also write an interactive HTML version")
    arguments = parser.parse_args()
    products, _ = run_example(arguments.typefileplot, arguments.typeplotback)
    if arguments.bokeh:
        run_bokeh_example(products, EXAMPLE_PATH / "visuals" / "pcat_image_forward_model.html")


if __name__ == "__main__":
    main()
