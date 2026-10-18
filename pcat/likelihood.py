"""Shared likelihood calculations for PCAT inference workflows."""

import numpy as np
from scipy.special import gammaln, xlogy


def poisson_log_likelihood(observed_counts, expected_counts, *, include_constant=False):
    """Return per-bin Poisson log likelihoods, optionally including log(k!).

    Data-only constants are omitted by default. Set ``include_constant=True``
    for normalized count probabilities, including count priors.
    """
    observed = np.asarray(observed_counts, dtype=float)
    expected = np.asarray(expected_counts, dtype=float)
    if observed.shape != expected.shape:
        raise ValueError("Poisson observations and expected counts must have matching shapes")
    if not np.isfinite(observed).all() or (observed < 0.0).any():
        raise ValueError("Poisson observations must be finite and nonnegative")
    if not np.isfinite(expected).all() or (expected < 0.0).any():
        raise ValueError("Poisson expected counts must be finite and nonnegative")
    result = xlogy(observed, expected) - expected
    if include_constant:
        result = result - gammaln(observed + 1.0)
    return result