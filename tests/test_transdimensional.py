import numpy as np
import pytest
from scipy.special import gammaln

from pcat.transdimensional import sample_catalog


def test_flat_likelihood_recovers_normalized_truncated_poisson_count_prior():
    mean = 1.4
    result = sample_catalog(lambda base, catalog: 0.0, [[-2.0, 2.0]], [[0.0, 5.0], [-3.0, 3.0]],
                            maximum_components=3, count_prior_mean=mean,
                            number_steps=40000, burn_steps=2000, thinning=2, seed=7)
    counts = np.arange(4)
    expected = np.exp(counts * np.log(mean) - gammaln(counts + 1))
    expected /= expected.sum()
    measured = np.bincount(result['component_counts'], minlength=4) / len(result['component_counts'])

    np.testing.assert_allclose(measured, expected, atol=0.025, rtol=0.0)
    assert result['acceptance_counts']['birth'] > 0
    assert result['acceptance_counts']['death'] > 0
    assert all(catalog.shape == (count, 2) for count, catalog in zip(result['component_counts'], result['catalog_samples']))


def test_catalog_sampler_reproducibility_and_zero_component_limit():
    arguments = dict(maximum_components=0, number_steps=100, burn_steps=20, thinning=4, seed=2)
    first = sample_catalog(lambda base, catalog: -0.5 * base[0]**2, [[-1.0, 1.0]], [[0.0, 1.0]], **arguments)
    second = sample_catalog(lambda base, catalog: -0.5 * base[0]**2, [[-1.0, 1.0]], [[0.0, 1.0]], **arguments)

    np.testing.assert_array_equal(first['base_samples'], second['base_samples'])
    assert np.all(first['component_counts'] == 0)
    assert first['proposal_counts']['birth'] == first['proposal_counts']['death'] == 0


def test_catalog_sampler_rejects_invalid_parameter_bounds():
    with pytest.raises(ValueError, match='bounds'):
        sample_catalog(lambda base, catalog: 0.0, [[1.0, 0.0]], [[0.0, 1.0]])