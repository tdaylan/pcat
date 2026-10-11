#!/usr/bin/env python3
"""Demonstrate PCAT's deterministic Gaussian-PSF image forward model."""

from tdpy.verbosity import print

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tdpy.cli import add_plot_arguments

from pcat.demo import configure_example_style, save_example_figure
from pcat.image import forward_model_image

EXAMPLE_PATH = Path(__file__).resolve().parent
IMAGE_SIZE_PIXELS = 51
PSF_SIGMA_PIXELS = 1.5
BACKGROUND_RATE = 5.0  # [electron s^-1 pixel^-1]


def create_simulated_scene():
    """Return an explicitly assumed two-source scene [electron s^-1 pixel^-1]."""
    coordinate = np.arange(IMAGE_SIZE_PIXELS)
    x_grid, y_grid = np.meshgrid(coordinate, coordinate)
    source_image = np.zeros((IMAGE_SIZE_PIXELS, IMAGE_SIZE_PIXELS))
    source_parameters = ((24.0, 25.0, 1000.0), (28.0, 25.0, 650.0))
    intrinsic_sigma_pixels = 0.35
    for x_position, y_position, flux_rate in source_parameters:
        profile = np.exp(
            -0.5
            * ((x_grid - x_position) ** 2 + (y_grid - y_position) ** 2)
            / intrinsic_sigma_pixels**2
        )
        source_image += flux_rate * profile / profile.sum()
    return source_image


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_plot_arguments(parser)
    parser.add_argument("--typeplotback", choices=("white", "dark"), default="white")
    arguments = parser.parse_args()
    run_example(arguments.typefileplot, arguments.typeplotback)


if __name__ == "__main__":
    main()
