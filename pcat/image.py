"""Reusable forward models for astronomical image formation."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve_fft


def make_gaussian_scene(
    image_shape: tuple[int, int],
    source_parameters: Sequence[tuple[float, float, float]],
    intrinsic_sigma_pixels: float,
) -> np.ndarray:
    """Render compact Gaussian sources with the requested in-frame fluxes.

    Each source is ``(x_pixel, y_pixel, flux)``. Flux uses the caller's units.
    """

    if len(image_shape) != 2 or any(
        not isinstance(size, (int, np.integer)) or size <= 0 for size in image_shape
    ):
        raise ValueError("image_shape must contain two positive integer dimensions")
    if not np.isfinite(intrinsic_sigma_pixels) or intrinsic_sigma_pixels <= 0.0:
        raise ValueError("intrinsic_sigma_pixels must be finite and positive")

    source_parameters = np.asarray(source_parameters, dtype=float)
    if (
        source_parameters.ndim != 2
        or source_parameters.shape[1] != 3
        or source_parameters.shape[0] == 0
        or not np.isfinite(source_parameters).all()
    ):
        raise ValueError(
            "source_parameters must be a nonempty finite array of (x, y, flux) rows"
        )
    height_pixels, width_pixels = image_shape
    x_positions, y_positions, fluxes = source_parameters.T
    if (
        np.any(x_positions < 0.0)
        or np.any(x_positions > width_pixels - 1)
        or np.any(y_positions < 0.0)
        or np.any(y_positions > height_pixels - 1)
        or np.any(fluxes < 0.0)
    ):
        raise ValueError(
            "source positions must be in-frame and source fluxes nonnegative"
        )

    x_grid, y_grid = np.meshgrid(
        np.arange(width_pixels),
        np.arange(height_pixels),
    )
    source_image = np.zeros(image_shape, dtype=float)
    for x_position, y_position, flux_rate in source_parameters:
        profile = np.exp(
            -0.5
            * ((x_grid - x_position) ** 2 + (y_grid - y_position) ** 2)
            / intrinsic_sigma_pixels**2
        )
        source_image += flux_rate * profile / profile.sum()
    return source_image


def forward_model_image(
    source_image: np.ndarray,
    psf_sigma_pixels: float,
    background: float = 0.0,
) -> dict[str, np.ndarray]:
    """Convolve a two-dimensional source scene with a normalized Gaussian PSF."""

    source_image = np.asarray(source_image, dtype=float)
    if source_image.ndim != 2:
        raise ValueError("source_image must be two-dimensional.")
    if not np.isfinite(source_image).all():
        raise ValueError("source_image must contain only finite values.")
    if not np.isfinite(psf_sigma_pixels) or psf_sigma_pixels <= 0.0:
        raise ValueError("psf_sigma_pixels must be finite and positive.")
    if not np.isfinite(background):
        raise ValueError("background must be finite.")

    psf_kernel = Gaussian2DKernel(x_stddev=psf_sigma_pixels)
    observed_image = convolve_fft(
        source_image,
        psf_kernel,
        boundary="fill",
        fill_value=0.0,
        normalize_kernel=True,
    )
    observed_image += background
    return {
        "source_image": source_image,
        "psf_kernel": np.asarray(psf_kernel.array),
        "observed_image": observed_image,
    }
