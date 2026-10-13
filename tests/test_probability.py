from types import SimpleNamespace

import numpy as np
import pytest

from pcat.main import pdfn_gaus, retr_gradmaps, retr_llik_bind, retr_llik_eval, retr_lpostotl, retr_lprielem


def test_pdfn_gaus_is_finite_at_mean():
    np.testing.assert_allclose(pdfn_gaus(0., 0., 1.), 1. / np.sqrt(2. * np.pi))


def test_retr_lprielem_supports_log_uniform_prior():
    model = SimpleNamespace(
        minmpara=SimpleNamespace(amplitude=1.),
        maxmpara=SimpleNamespace(amplitude=100.),
    )
    gdat = SimpleNamespace(fitt=model)
    features = np.array([1., 10.])

    result = retr_lprielem(
        gdat,
        'fitt',
        0,
        0,
        'amplitude',
        'logt',
        np.array([]),
        [{'amplitude': features}],
        [features.size],
    )

    expected = np.sum(np.log(1. / np.log(100.) / features))
    np.testing.assert_allclose(result, expected)


def test_retr_llik_eval_invokes_custom_likelihood():
    calls = []

    def custom_likelihood(gdat, strgmodl, cntpmodl):
        calls.append((gdat, strgmodl, cntpmodl))
        return np.full_like(cntpmodl, -7.)

    gdat = SimpleNamespace(retr_llik=custom_likelihood)
    model_counts = np.array([1., 2.])

    result = retr_llik_eval(gdat, 'fitt', model_counts)

    np.testing.assert_array_equal(result, [-7., -7.])
    assert calls == [(gdat, 'fitt', model_counts)]


def test_poisson_likelihood_preserves_zero_model_rates():
    gdat = SimpleNamespace(liketype='pois', cntpdata=np.array([0.0, 1.0]))

    result = retr_llik_bind(gdat, 'fitt', np.array([0.0, 0.0]))

    assert result[0] == 0.0
    assert np.isneginf(result[1])


def test_poisson_likelihood_rejects_invalid_model_counts():
    gdat = SimpleNamespace(liketype='pois', cntpdata=np.array([0.0, 1.0]))

    with pytest.raises(ValueError, match='nonnegative'):
        retr_llik_bind(gdat, 'fitt', np.array([0.0, -1.0]))


def test_relevance_gradient_rejects_nonfinite_input_maps():
    gdat = SimpleNamespace(numbsidecart=3, numbpixlcart=9, sizepixl=1.0)
    maps = np.ones(9)
    maps[4] = np.nan

    with pytest.raises(ValueError, match='non-finite pixels'):
        retr_gradmaps(gdat, maps)


def test_retr_lpostotl_uses_requested_target():
    assert retr_lpostotl(-2., -5., 'prio') == -2.
    assert retr_lpostotl(-2., -5., 'post') == -7.
