"""Shared likelihood calculations for PCAT inference workflows."""

import numpy as np
from scipy.special import xlogy


def poisson_log_likelihood(observed_counts, expected_counts):
    """Return per-bin Poisson log likelihoods, omitting data-only constants."""
    observed = np.asarray(observed_counts, dtype=float)
    expected = np.asarray(expected_counts, dtype=float)
    if observed.shape != expected.shape:
        raise ValueError("Poisson observations and expected counts must have matching shapes")
    if not np.isfinite(observed).all() or (observed < 0.0).any():
        raise ValueError("Poisson observations must be finite and nonnegative")
    if not np.isfinite(expected).all() or (expected < 0.0).any():
        raise ValueError("Poisson expected counts must be finite and nonnegative")
    return xlogy(observed, expected) - expected