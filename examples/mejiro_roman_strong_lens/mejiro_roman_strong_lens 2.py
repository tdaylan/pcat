#!/usr/bin/env python3
"""Infer a mejiro-simulated Roman strong lens with PCAT.

mejiro (Wedig et al., https://github.com/AstroMusers/mejiro) renders the SLSim
galaxy-galaxy lens ``SampleGG`` in the Roman Wide Field Instrument (WFI) F129 band
with its Roman point-spread function (PSF) width, zero point, and sky background.
PCAT samples the singular isothermal ellipsoid (SIE) lens, external shear, and
source position from that Poisson exposure, starting from a random prior draw.
Requires ``pip install -e <mejiro>`` and ``pip install roman-technical-information``.
"""

from __future__ import annotations

from tdpy.verbosity import print

import argparse
import shutil
from pathlib import Path

import numpy as np
from tdpy.cli import add_plot_arguments

from pcat.diagnostics import gelman_rubin
from pcat.main import retr_pathrun
from pcat.mejiro_lens import LABELS, render_mejiro_counts, run_mejiro_lens_inference, simulate_mejiro_exposure
from pcat.plotting import (
    animation_phase_label, animation_states, make_image_sequence_animation,
    plot_lens_image_fit, plot_lens_parameter_recovery,
)

EXAMPLE_PATH = Path(__file__).resolve().parent
VISUAL_PATH = EXAMPLE_PATH / "visuals"
RUN_NAME = "mejiro_roman_strong_lens"
NUMBER_CHAINS = 4


def run_example(smoke: bool = False, typefileplot: str = "png") -> dict:
    """Simulate the mejiro exposure, sample it with PCAT, and write figures."""
    specification = simulate_mejiro_exposure()
    run_root = Path(retr_pathrun(EXAMPLE_PATH, RUN_NAME))
    if run_root.exists():
        print(f"Removing cached PCAT run {run_root}...")
        shutil.rmtree(run_root)
    posterior = run_mejiro_lens_inference(
        specification, EXAMPLE_PATH, RUN_NAME, numbproc=NUMBER_CHAINS, numbframanim=24,
        numbswep=6000 if smoke else 60000, numbburn=4000 if smoke else 40000,
        numbsamp=500 if smoke else 2000,
    )
    # retained samples are ordered sample-major across chains
    draws = np.asarray(posterior.listpostparagenrscalbase)
    chains = draws.reshape(-1, NUMBER_CHAINS, draws.shape[-1])  # [sample, chain, parameter]
    median_counts = render_mejiro_counts(specification, np.median(draws, axis=0))  # [counts pixel^-1]
    figure_paths = [
        plot_lens_image_fit(VISUAL_PATH / f"{RUN_NAME}_image_fit", specification["observed_counts"],
                            median_counts, median_counts, specification["pixel_scale"], typefileplot),
        plot_lens_parameter_recovery(VISUAL_PATH / f"{RUN_NAME}_posterior", draws,
                                     specification["true_parameters"], typefileplot, labels=LABELS),
    ]
    figure_paths.append(write_replicated_image_animation(posterior, specification))
    summary = {
        "maximum_rhat": max(gelman_rubin(chains[:, :, index]) for index in range(chains.shape[-1])),
        "acceptance_fraction": float(np.mean(posterior.listpostboolpropaccp)),
        "true_parameters": specification["true_parameters"],
        "posterior_medians": np.median(draws, axis=0),
        "posterior_stdvs": draws.std(axis=0),
        "figure_paths": figure_paths,
    }
    return summary


def write_replicated_image_animation(posterior: object, specification: dict) -> Path:
    """Animate replicated mejiro exposures from the initial prior draw to posterior samples."""
    # one fixed noise realization, so frames differ only through the lens model
    noise_seed = 2027
    images, labels = [], []
    for snapshot in animation_states(posterior):
        counts = render_mejiro_counts(specification, snapshot["paragenrscalfull"][:len(LABELS)])
        # counts in units of the sky level, with an arcsinh stretch to show faint arcs beside the bright lens
        images.append(np.arcsinh(np.random.default_rng(noise_seed).poisson(counts) / specification["sky_counts"]))
        labels.append(animation_phase_label(snapshot))
    return make_image_sequence_animation(
        images, labels, VISUAL_PATH / f"{RUN_NAME}_replicated.gif", duration_ms=300,
        title="mejiro Roman F129 lens, PCAT replicated exposure",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="Run a short pipeline check.")
    add_plot_arguments(parser)
    arguments = parser.parse_args()
    summary = run_example(arguments.smoke, arguments.typefileplot)
    for label, truth, median, stdv in zip(LABELS, summary["true_parameters"], summary["posterior_medians"],
                                          summary["posterior_stdvs"]):
        print(f"{label}: injected {truth:+.4f}, posterior {median:+.4f} +/- {stdv:.4f}")
    print(f"Maximum Gelman-Rubin R-hat {summary['maximum_rhat']:.3f}, "
          f"acceptance fraction {summary['acceptance_fraction']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
