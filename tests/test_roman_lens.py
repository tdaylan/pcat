import numpy as np
import pytest
import pickle
from types import SimpleNamespace

from pcat.roman_lens import (
    RomanLensConfig,
    gaussian_lens_log_likelihood,
    infer_catalog_probability,
    poisson_lens_log_likelihood,
    render_lens,
    render_lens_counts,
    run_lens_image_pipeline,
    simulate_population,
    summarize_population,
)


def test_lens_image_pipeline_writes_fit_and_parameter_visuals(tmp_path, monkeypatch):
    from pcat import main

    config = RomanLensConfig(number_side=24, pixel_scale=0.2, psf_fwhm=0.7)
    truth = np.array((1.2, 0.1, -0.2))  # [arcsec]
    observed = render_lens_counts(config, truth, 0.3, 0.7, 0.4)
    draws = np.tile(truth, (20, 1)) + np.linspace(-0.02, 0.02, 20)[:, None]
    posterior = SimpleNamespace(listpostparagenrscalbase=draws)
    def sample_with_visible_ring(**configuration):
        assert configuration["prior_minima"][0] < truth[0] < configuration["prior_maxima"][0]
        assert configuration["prior_maxima"][0] < config.number_side * config.pixel_scale / 2
        return posterior

    monkeypatch.setattr(main, "sample", sample_with_visible_ring)

    result = run_lens_image_pipeline(
        config=config,
        observed_counts=observed,
        source_size=0.3,
        source_axis_ratio=0.7,
        source_angle=0.4,
        true_parameters=truth,
        output_root=tmp_path,
        run_name="lens_test",
        visual_stem="lens_test",
    )

    assert result.posterior is posterior
    assert result.image_fit_path.is_file()
    assert result.parameter_path.is_file()
    assert result.model_counts.shape == observed.shape


def test_cluster_poisson_likelihood_prefers_the_injected_lens():
    config = RomanLensConfig(number_side=32, pixel_scale=0.2, psf_fwhm=0.7)
    true_parameters = np.array((1.2, 0.1, -0.2))  # [arcsec]
    source_size = 0.3  # [arcsec]
    observed = render_lens_counts(config, true_parameters, source_size, 0.7, 0.4)
    state = SimpleNamespace(
        lens_config=config,
        lens_observed_counts=observed,
        lens_source_size=source_size,
        lens_source_axis_ratio=0.7,
        lens_source_angle=0.4,
    )

    assert poisson_lens_log_likelihood(state, "fitt", true_parameters) > poisson_lens_log_likelihood(
        state, "fitt", np.array((1.6, 0.1, -0.2))
    )
    assert pickle.loads(pickle.dumps(poisson_lens_log_likelihood)) is poisson_lens_log_likelihood


def test_roman_poisson_likelihood_retains_omitted_constants_and_zero_pixels(monkeypatch):
    import pcat.roman_lens as roman_lens
    from scipy.special import gammaln
    from scipy.stats import poisson

    model = np.array([[0.0, 2.0], [3.0, 7.0]])
    observed = np.array([[0.0, 1.0], [4.0, 9.0]])
    monkeypatch.setattr(roman_lens, "render_lens_counts", lambda *args: model)
    state = SimpleNamespace(
        lens_config=None, lens_observed_counts=observed,
        lens_source_size=0.3, lens_source_axis_ratio=0.7, lens_source_angle=0.4,
    )
    expected = np.sum(poisson.logpmf(observed, model) + gammaln(observed + 1.0))
    assert poisson_lens_log_likelihood(state, "fitt", np.array([])) == pytest.approx(expected, rel=1e-14)
    observed[0, 0] = 1.0
    assert np.isneginf(poisson_lens_log_likelihood(state, "fitt", np.array([])))
    model[0, 0] = -1.0
    with pytest.raises(ValueError, match="nonnegative"):
        poisson_lens_log_likelihood(state, "fitt", np.array([]))


def test_simulated_roman_source_produces_bright_einstein_arc():
    config = RomanLensConfig()
    image = render_lens(config, 0.9, 0.06, -0.04, 0.09, 0.8, 0.35)
    row, column = np.indices(image.shape)
    radius = np.hypot(row - 15, column - 15) * config.pixel_scale  # [arcsec]
    assert image[(radius > 0.8) & (radius < 1.2)].mean() > 3 * image[radius < 0.2].mean()


def test_elliptical_psf_render_preserves_flux_and_uses_orientation():
    circular_config = RomanLensConfig(
        number_side=64, pixel_scale=0.08, psf_fwhm=0.16, source_counts=1.0
    )
    elliptical_config = RomanLensConfig(
        number_side=64,
        pixel_scale=0.08,
        psf_fwhm=0.16,
        psf_fwhm_minor=0.10,
        psf_angle=np.pi / 4,
        source_counts=1.0,
    )
    parameters = (1.2, 0.0, 0.0, 0.12, 0.65, 0.3)

    circular = render_lens(circular_config, *parameters)
    elliptical = render_lens(elliptical_config, *parameters)

    assert np.isclose(elliptical.sum(), circular.sum(), rtol=1e-3)
    assert not np.allclose(elliptical, circular)


def test_gaussian_lens_likelihood_uses_variance_and_ignores_invalid_pixels():
    config = RomanLensConfig(number_side=24, background=0.0)
    true_parameters = np.array((0.9, 0.04, -0.03))  # [arcsec]
    image = render_lens_counts(config, true_parameters, 0.1, 0.8, 0.2)
    variance = np.full_like(image, 4.0)  # [nJy^2 pixel^-2]
    variance[0, 0] = np.nan
    state = SimpleNamespace(
        lens_config=config,
        lens_observed_image=image,
        lens_variance=variance,
        lens_source_size=0.1,
        lens_source_axis_ratio=0.8,
        lens_source_angle=0.2,
    )

    assert gaussian_lens_log_likelihood(state, "fitt", true_parameters) > gaussian_lens_log_likelihood(
        state, "fitt", np.array((1.2, 0.04, -0.03))
    )


def test_injected_perturber_raises_catalog_probability():
    config = RomanLensConfig(number_side=21, source_counts=2.0e5, number_candidates=1)
    parameters = (config, 0.8, 0.05, -0.03, 0.09, 0.7, 0.4)
    macro = render_lens(*parameters)
    perturbed = render_lens(*parameters, subhalo=(0.8, 0.0, 0.03))
    variance = perturbed + config.background + config.read_noise**2
    template = (perturbed - macro)[None, :, :]
    null_probability = infer_catalog_probability(macro + config.background, macro + config.background, template, variance)[0]
    perturbed_probability = infer_catalog_probability(perturbed + config.background, macro + config.background, template, variance)[0]

    assert np.isfinite(perturbed_probability)
    assert perturbed_probability > null_probability


@pytest.mark.parametrize('signal', [-4.0, 0.0, 1.0, 10.0, 1e4])
def test_catalog_probability_preserves_penalties_amplitude_and_clipping(signal):
    data = np.full((2, 2), signal)
    templates = np.stack((np.ones((2, 2)), np.eye(2)))
    amplitude = max(signal, 0.0)
    log_odds = 0.5 * (4.0 * amplitude**2 - 3.0 * np.log(4.0) - 2.0 * np.log(2.0))
    expected = 1.0 / (1.0 + np.exp(-np.clip(log_odds, -40.0, 40.0)))

    probability, best_index, fitted_amplitude = infer_catalog_probability(
        data, np.zeros_like(data), templates, np.ones_like(data),
    )

    assert type(probability) is float
    assert probability == pytest.approx(expected, rel=1e-14)
    assert best_index == 0
    assert fitted_amplitude == amplitude


def test_zero_information_template_has_finite_catalog_probability():
    data = np.ones((5, 5))
    templates = np.zeros((1, 5, 5))

    probability, best_index, amplitude = infer_catalog_probability(
        data, data, templates, np.ones_like(data)
    )

    assert np.isfinite(probability)
    assert probability < 0.5
    assert best_index == 0
    assert amplitude == 0.0


def test_population_records_generated_scene_diagnostics():
    records, examples = simulate_population(number_lenses=4, seed=814)
    summary = summarize_population(records)

    assert records[0]["injected_signal_to_noise"] > 0.0
    assert records[-1]["injected_signal_to_noise"] == 0.0
    assert 0.6 <= records[0]["macro_einstein_radius_arcsec"] <= 1.2
    assert 0.06 <= records[0]["source_size_arcsec"] <= 0.14
    assert examples["lens_0_data"].shape == (31, 31)
    assert examples["lens_0_macro"].shape == (31, 31)
    assert np.any(examples["lens_0_residual"] != 0.0)
    assert summary["true_positive_rate_lower"] <= summary["true_positive_rate"]
    assert summary["true_positive_rate"] <= summary["true_positive_rate_upper"]
    assert summary["false_positive_rate_lower"] <= summary["false_positive_rate"]
    assert summary["false_positive_rate"] <= summary["false_positive_rate_upper"]


def test_population_threshold_controls_detection_rates():
    records, _ = simulate_population(number_lenses=20, seed=814)
    permissive = summarize_population(records, detection_threshold=0.1)
    conservative = summarize_population(records, detection_threshold=0.9)

    assert permissive["true_positive_rate"] >= conservative["true_positive_rate"]
    assert permissive["false_positive_rate"] >= conservative["false_positive_rate"]


def test_population_rejects_invalid_detection_threshold():
    records, _ = simulate_population(number_lenses=4, seed=814)

    with pytest.raises(ValueError, match="between zero and one"):
        summarize_population(records, detection_threshold=1.1)


@pytest.mark.parametrize('threshold', [
    np.complex128(0.5 + 0.1j), 0.5 + 0.0j, np.array([0.5]), np.array([[0.5]]),
    [0.5], '0.5', None, np.ma.masked, np.ma.array(0.5, mask=True),
    np.nan, np.inf, -np.inf, -0.1, 1.1,
])
@pytest.mark.parametrize('plot', [False, True])
def test_population_rejects_nonreal_or_nonscalar_thresholds_before_rendering(tmp_path, monkeypatch, threshold, plot):
    import warnings
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid threshold reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.5, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(ValueError, match='detection_threshold'):
        warnings.simplefilter('error')
        if plot:
            plot_detection_diagnostic(records, tmp_path / 'bad_scalar_threshold.png', detection_threshold=threshold)
        else:
            summarize_population(records, detection_threshold=threshold)
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('threshold', [0, 1, False, True, np.float32(0.5), np.float64(0.5),
                                     np.array(0.5), np.ma.array(0.5, mask=False)])
def test_population_preserves_real_scalar_thresholds(threshold):
    import warnings

    records = [dict(has_subhalo=True, posterior_one=0.5, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    expected = summarize_population(records, detection_threshold=float(threshold))
    originals = [record.copy() for record in records]
    with warnings.catch_warnings(), np.errstate(all='raise'):
        warnings.simplefilter('error')
        assert summarize_population(records, detection_threshold=threshold) == expected
    assert records == originals


def test_population_requires_injected_and_null_cohorts():
    with pytest.raises(ValueError, match="injected and null"):
        simulate_population(number_lenses=1)


@pytest.mark.parametrize('bad_probability', [
    np.nan, np.inf, -np.inf, -0.01, 1.01, np.ma.masked,
    np.ma.array(0.5, mask=True), np.array([0.5]),
    0.5 + 0.0j, np.complex64(0.5 + 0.1j), np.complex128(0.5 + 0.1j),
    np.complex128(0.5 + 0.0j), np.array(0.5 + 0.0j), np.ma.array(0.5 + 0.1j, mask=False),
])
@pytest.mark.parametrize('record_index', [0, 1])
def test_population_summary_rejects_invalid_probabilities(bad_probability, record_index):
    import warnings

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    records[record_index]['posterior_one'] = bad_probability
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(ValueError, match='Posterior probabilities'):
        warnings.simplefilter('error')
        summarize_population(records)


@pytest.mark.parametrize('threshold', [0.0, 0.5, 1.0])
@pytest.mark.parametrize('container', ['numeric', 'unmasked', 'string'])
def test_population_summary_preserves_probability_boundaries(threshold, container):
    from pcat.diagnostics import binomial_wilson_interval

    values = [0.0, 0.5, 1.0]
    if container == 'unmasked':
        probabilities = [np.ma.array(value, mask=False) for value in values]
    elif container == 'string':
        probabilities = [str(value) for value in values]
    else:
        probabilities = values
    records = [dict(has_subhalo=truth, posterior_one=value, injected_index=0 if truth else -1, recovered_index=0)
               for truth in (True, False) for value in probabilities]
    originals = [record.copy() for record in records]
    summary = summarize_population(records, detection_threshold=threshold)
    detections = sum(value >= threshold for value in values)
    lower, upper = binomial_wilson_interval(detections, 3)

    assert summary['number_lenses'] == 6
    assert summary['number_injected'] == 3
    assert summary['localization_rate'] == 1.0
    for cohort in ('true_positive', 'false_positive'):
        assert summary[f'{cohort}_rate'] == detections / 3
        assert summary[f'{cohort}_rate_lower'] == lower
        assert summary[f'{cohort}_rate_upper'] == upper
    assert summary['mean_posterior_injected'] == summary['mean_posterior_null'] == 0.5
    assert records == originals


@pytest.mark.parametrize('bad_probability', [np.ma.masked, np.complex128(0.9 + 0.4j), 0.9 + 0.0j])
def test_population_plot_validates_probabilities_before_conversion_or_rendering(tmp_path, monkeypatch, bad_probability):
    import warnings
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid probability reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=bad_probability, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    for record, signal_to_noise in zip(records, (5.0, 0.0)):
        record['injected_signal_to_noise'] = signal_to_noise
    with warnings.catch_warnings(), pytest.raises(ValueError, match='Posterior probabilities'):
        warnings.simplefilter('error')
        plot_detection_diagnostic(records, tmp_path / 'invalid_probability.png')
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('bad_truth', [
    'False', 'True', '0', '1', None, np.nan, np.inf, -np.inf,
    -1, 2, 0.5, 1.0 + 0.0j, np.ma.masked,
    np.ma.array(False, mask=True), np.array([True]),
])
@pytest.mark.parametrize('record_index', [0, 1])
def test_population_summary_rejects_invalid_truth_flags(bad_truth, record_index):
    import warnings

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    records[record_index]['has_subhalo'] = bad_truth
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(ValueError, match='has_subhalo'):
        warnings.simplefilter('error')
        summarize_population(records)


@pytest.mark.parametrize('truth_flags', [
    (True, False), (np.bool_(True), np.bool_(False)), (1, 0), (1.0, 0.0),
    (np.int64(1), np.uint8(0)), (np.float32(1.0), np.float64(0.0)),
    (np.array(1), np.array(0)),
    (np.ma.array(True, mask=False), np.ma.array(False, mask=False)),
])
def test_population_summary_preserves_binary_truth_flags(truth_flags):
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    expected = summarize_population(records)
    for record, flag in zip(records, truth_flags):
        record['has_subhalo'] = flag
    originals = [record.copy() for record in records]
    assert summarize_population(records) == expected
    assert records == originals


@pytest.mark.parametrize('bad_truth', ['False', np.nan, np.ma.masked])
def test_population_plot_validates_truth_flags_before_conversion_or_rendering(tmp_path, monkeypatch, bad_truth):
    import warnings
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid truth flag reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=bad_truth, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    with warnings.catch_warnings(), pytest.raises(ValueError, match='has_subhalo'):
        warnings.simplefilter('error')
        plot_detection_diagnostic(records, tmp_path / 'invalid_truth.png')
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('bad_index', [
    0.1, -0.1, 0.9, np.nan, np.inf, -np.inf, np.ma.masked,
    np.ma.array(0, mask=True), np.array([0]), 1.0 + 0.0j, None,
])
@pytest.mark.parametrize('field', ['injected_index', 'recovered_index'])
@pytest.mark.parametrize('record_index', [0, 1])
def test_population_summary_rejects_invalid_localization_indices(bad_index, field, record_index):
    import warnings

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    records[record_index][field] = bad_index
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(ValueError, match=field):
        warnings.simplefilter('error')
        summarize_population(records)


@pytest.mark.parametrize('convert', [int, float, str, np.asarray, np.int64,
                                    lambda value: np.ma.array(value, mask=False)])
def test_population_summary_preserves_exact_localization_indices(convert):
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0),
               dict(has_subhalo=True, posterior_one=0.8, injected_index=1, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=1)]
    expected = summarize_population(records)
    assert expected['localization_rate'] == 0.5
    for record in records:
        for field in ('injected_index', 'recovered_index'):
            record[field] = convert(record[field])
    originals = [record.copy() for record in records]
    assert summarize_population(records) == expected
    assert records == originals


@pytest.mark.parametrize('convert', [np.uint64, np.float64])
@pytest.mark.parametrize('field', ['injected_index', 'recovered_index'])
def test_population_summary_rejects_out_of_range_localization_indices(convert, field):
    import warnings

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0)]
    records[0][field] = convert(np.iinfo(np.int_).max + 1)
    originals = [record.copy() for record in records]
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(OverflowError):
        warnings.simplefilter('error')
        summarize_population(records)
    assert records == originals


@pytest.mark.parametrize('field', ['injected_index', 'recovered_index'])
def test_population_plot_validates_localization_indices_before_rendering(tmp_path, monkeypatch, field):
    import warnings
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid localization index reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    records[0][field] = 0.1
    with warnings.catch_warnings(), pytest.raises(ValueError, match=field):
        warnings.simplefilter('error')
        plot_detection_diagnostic(records, tmp_path / 'invalid_index.png')
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('threshold', [0.0, 0.1, 0.5, 0.9, 1.0])
@pytest.mark.parametrize('with_examples', [False, True])
@pytest.mark.parametrize('typefileplot', ['png', 'pdf'])
def test_population_plot_uses_shared_detection_threshold(tmp_path, monkeypatch, threshold, with_examples, typefileplot):
    from matplotlib import pyplot as plt
    from pcat import plotting

    records = [dict(has_subhalo=truth, posterior_one=probability, injected_index=0 if truth else -1,
                    recovered_index=0, injected_signal_to_noise=5.0 if truth else 0.0)
               for truth, probability in ((True, 0.2), (True, 0.8), (False, 0.1), (False, 0.6))]
    expected = summarize_population(records, detection_threshold=threshold)
    originals = [record.copy() for record in records]
    examples = None
    if with_examples:
        data = np.array([[2.0, 3.0], [4.0, 5.0]])  # [electron pixel^-1]
        examples = {'lens_0_data': data, 'lens_0_macro': data - 0.5,
                    'lens_0_residual': np.full((2, 2), 0.5)}
    captured = []
    save = plotting.save_figure
    before = set(plt.get_fignums())

    def capture(figure, path, **kwargs):
        axis = next(axis for axis in figure.axes if axis.get_title() == 'Catalog inference')
        assert len(axis.lines) == 1
        assert axis.lines[0].get_label() == 'Detection threshold'
        np.testing.assert_array_equal(axis.lines[0].get_ydata(), [threshold, threshold])
        text = axis.texts[0].get_text().splitlines()
        for line, label, key in zip(text[:2], ('TPR', 'FPR'), ('true_positive', 'false_positive')):
            assert line == (f"{label} {100 * expected[f'{key}_rate']:.0f}% "
                            f"[{100 * expected[f'{key}_rate_lower']:.0f}, {100 * expected[f'{key}_rate_upper']:.0f}]")
        assert text[2] == 'Wilson intervals ($z=1$)'
        captured.append(figure)
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_figure', capture)
    options = {} if threshold == 0.5 else {'detection_threshold': threshold}
    plotting.plot_detection_diagnostic(records, tmp_path / f'population.{typefileplot}', examples, **options)
    assert len(captured) == 1
    assert records == originals
    assert set(plt.get_fignums()) == before


@pytest.mark.parametrize('threshold', [-0.1, 1.1, np.nan, np.inf, -np.inf])
def test_population_plot_rejects_invalid_detection_threshold_before_rendering(tmp_path, monkeypatch, threshold):
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid detection threshold reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    with pytest.raises(ValueError, match='detection_threshold'):
        plot_detection_diagnostic(records, tmp_path / 'bad_threshold.png', detection_threshold=threshold)
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('field', ['data', 'macro', 'residual'])
@pytest.mark.parametrize('shape', [(4, 2), (2, 1), (3, 3), (0, 4), (2, 0), (8,), (2, 4, 3)])
def test_population_plot_rejects_incompatible_representative_grids(tmp_path, monkeypatch, field, shape):
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid representative grid reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    examples = {f'lens_0_{key}': np.ones((2, 4)) for key in ('data', 'macro', 'residual')}
    examples[f'lens_0_{field}'] = np.ones(shape)
    with pytest.raises(ValueError, match='nonempty two-dimensional shape'):
        plot_detection_diagnostic(records, tmp_path / 'invalid_representative_grid.png', examples)
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('shape', [(2, 4), (4, 2), (3, 3)])
@pytest.mark.parametrize('typefileplot', ['png', 'pdf'])
def test_population_plot_preserves_representative_pixel_grids(tmp_path, monkeypatch, shape, typefileplot):
    from matplotlib import pyplot as plt
    from pcat import plotting

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    data = np.arange(1.0, np.prod(shape) + 1.0).reshape(shape)  # [electron pixel^-1]
    examples = {'lens_0_data': data, 'lens_0_macro': data - 0.5,
                'lens_0_residual': np.full(shape, 0.5)}
    originals = {key: value.copy() for key, value in examples.items()}
    captured = []
    save = plotting.save_figure
    before = set(plt.get_fignums())

    def capture(figure, path, **kwargs):
        axes = [axis for axis in figure.axes if axis.images]
        assert len(axes) == 3
        for axis, key in zip(axes, ('lens_0_data', 'lens_0_macro', 'lens_0_residual')):
            np.testing.assert_array_equal(axis.images[0].get_array(), examples[key])
            assert axis.images[0].get_array().shape == shape
            assert axis.get_xlabel() == 'Detector x [pixel]'
            assert axis.get_ylabel() == 'Detector y [pixel]'
        captured.append(figure)
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_figure', capture)
    plotting.plot_detection_diagnostic(records, tmp_path / f'pixel_grid.{typefileplot}', examples)
    assert len(captured) == 1
    for key, original in originals.items():
        np.testing.assert_array_equal(examples[key], original)
    assert set(plt.get_fignums()) == before


@pytest.mark.parametrize('hidden_value', [1e12, np.nan, np.inf, -np.inf])
@pytest.mark.parametrize('typefileplot', ['png', 'pdf'])
def test_population_plot_scales_only_visible_representative_pixels(tmp_path, monkeypatch, hidden_value, typefileplot):
    import warnings
    from matplotlib import pyplot as plt
    from pcat import plotting

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    examples = {
        'lens_0_data': np.ma.array([[hidden_value, 2.0], [3.0, 4.0]], mask=[[True, False], [False, False]]),
        'lens_0_macro': np.ma.array([[0.5, 1.0], [hidden_value, 2.0]], mask=[[False, False], [True, False]]),
        'lens_0_residual': np.ma.array([[0.5, hidden_value], [-1.0, 2.0]], mask=[[False, True], [False, False]]),
    }
    originals = {key: value.copy() for key, value in examples.items()}
    minimum = min(np.percentile([2.0, 3.0, 4.0], 1.0), np.percentile([0.5, 1.0, 2.0], 1.0))  # [electron pixel^-1]
    maximum = max(np.percentile([2.0, 3.0, 4.0], 99.5), np.percentile([0.5, 1.0, 2.0], 99.5))  # [electron pixel^-1]
    captured = []
    save = plotting.save_figure
    before = set(plt.get_fignums())

    def capture(figure, path, **kwargs):
        axes = [axis for axis in figure.axes if axis.images]
        assert len(axes) == 3
        for axis, key, limits in zip(axes, examples, ((minimum, maximum), (minimum, maximum), (-2.0, 2.0))):
            image = axis.images[0]
            assert image.get_clim() == limits
            np.testing.assert_array_equal(np.ma.getmaskarray(image.get_array()), examples[key].mask)
            np.testing.assert_array_equal(image.get_array().compressed(), examples[key].compressed())
        captured.append(figure)
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_figure', capture)
    with warnings.catch_warnings(), np.errstate(all='raise'):
        warnings.simplefilter('error')
        plotting.plot_detection_diagnostic(records, tmp_path / f'masked_pixel_scales.{typefileplot}', examples)
    assert len(captured) == 1
    for key, original in originals.items():
        np.testing.assert_array_equal(examples[key].data, original.data)
        np.testing.assert_array_equal(examples[key].mask, original.mask)
    assert set(plt.get_fignums()) == before


@pytest.mark.parametrize('field', ['data', 'macro', 'residual'])
@pytest.mark.parametrize('bad_image', [np.full((2, 2), np.nan), np.full((2, 2), np.inf),
                                     np.full((2, 2), -np.inf), np.ma.masked_all((2, 2)),
                                     np.ones((2, 2), dtype=complex)])
def test_population_plot_rejects_invalid_visible_representative_pixels(tmp_path, monkeypatch, field, bad_image):
    import warnings
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid visible pixels reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    examples = {f'lens_0_{key}': np.ones((2, 2)) for key in ('data', 'macro', 'residual')}
    examples[f'lens_0_{field}'] = bad_image
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(ValueError, match='finite real unmasked pixels'):
        warnings.simplefilter('error')
        plot_detection_diagnostic(records, tmp_path / 'bad_visible_pixels.png', examples)
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('dtype', [np.int8, np.uint8, np.int64, np.uint64])
@pytest.mark.parametrize('masked_extreme', [False, True])
@pytest.mark.parametrize('typefileplot', ['png', 'pdf'])
def test_population_plot_scales_integer_residuals_without_overflow(tmp_path, monkeypatch, dtype, masked_extreme, typefileplot):
    import warnings
    from matplotlib import pyplot as plt
    from pcat import plotting

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    limits = np.iinfo(dtype)
    extreme = limits.min if np.dtype(dtype).kind == 'i' else limits.max  # [electron pixel^-1]
    values = np.array([[extreme, 1], [2, 3]], dtype=dtype)  # [electron pixel^-1]
    mask = np.array([[masked_extreme, False], [False, False]])
    residual = np.ma.array(values, mask=mask) if masked_extreme else values
    original = np.ma.asarray(residual).copy()
    magnitude = float(max(abs(int(value)) for value in values[~mask]))  # [electron pixel^-1]
    examples = {'lens_0_data': np.array([[2.0, 3.0], [4.0, 5.0]]),
                'lens_0_macro': np.ones((2, 2)), 'lens_0_residual': residual}
    captured = []
    save = plotting.save_figure
    before = set(plt.get_fignums())

    def capture(figure, path, **kwargs):
        axis = next(axis for axis in figure.axes if axis.get_title() == 'Data - macro model')
        image = axis.images[0]
        assert image.get_clim() == (-magnitude, magnitude)
        assert image.get_array().dtype == values.dtype
        np.testing.assert_array_equal(np.ma.getmaskarray(image.get_array()), mask)
        np.testing.assert_array_equal(image.get_array().compressed(), values[~mask])
        captured.append(figure)
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_figure', capture)
    with warnings.catch_warnings(), np.errstate(all='raise'):
        warnings.simplefilter('error')
        plotting.plot_detection_diagnostic(records, tmp_path / f'integer_residual_scale.{typefileplot}', examples)
    assert len(captured) == 1
    np.testing.assert_array_equal(np.ma.getdata(residual), original.data)
    np.testing.assert_array_equal(np.ma.getmaskarray(residual), np.ma.getmaskarray(original))
    assert set(plt.get_fignums()) == before


@pytest.mark.parametrize('example_prefix', [None, 'lens_0', 'lens_2'])
@pytest.mark.parametrize('typefileplot', ['png', 'pdf'])
def test_population_plot_selects_representative_example_without_mutation(tmp_path, monkeypatch, example_prefix, typefileplot):
    from matplotlib import pyplot as plt
    from pcat import plotting

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    examples = {}
    for prefix, offset, difference in (('lens_2', 10.0, 3.0), ('lens_0', 0.0, 0.5)):  # offset, difference [electron pixel^-1]
        data = np.arange(1.0, 9.0).reshape(2, 4) + offset  # [electron pixel^-1]
        examples[f'{prefix}_data'] = data
        examples[f'{prefix}_macro'] = data - difference
        examples[f'{prefix}_residual'] = np.full(data.shape, difference)
    originals = {key: value.copy() for key, value in examples.items()}
    expected_prefix = 'lens_0' if example_prefix is None else example_prefix
    captured = []
    save = plotting.save_figure
    before = set(plt.get_fignums())

    def capture(figure, path, **kwargs):
        axes = [axis for axis in figure.axes if axis.images]
        assert len(axes) == 3
        for axis, suffix in zip(axes, ('data', 'macro', 'residual')):
            np.testing.assert_array_equal(axis.images[0].get_array(), examples[f'{expected_prefix}_{suffix}'])
        captured.append(figure)
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_figure', capture)
    options = {} if example_prefix is None else {'example_prefix': example_prefix}
    plotting.plot_detection_diagnostic(records, tmp_path / f'selected_example.{typefileplot}', examples, **options)
    assert len(captured) == 1
    for key, original in originals.items():
        np.testing.assert_array_equal(examples[key], original)
    assert set(plt.get_fignums()) == before


@pytest.mark.parametrize('case', ['none', 'empty', 'no_data', 'unknown'])
def test_population_plot_rejects_unavailable_example_selection_before_rendering(tmp_path, monkeypatch, case):
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Unavailable example reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    examples = {f'lens_0_{key}': np.ones((2, 2)) for key in ('data', 'macro', 'residual')}
    prefix = 'lens_0'
    if case == 'none':
        examples = None
    elif case == 'empty':
        examples = {}
    elif case == 'no_data':
        examples.pop('lens_0_data')
    else:
        prefix = 'lens_missing'
    with pytest.raises(ValueError, match='example_prefix|Representative examples'):
        plot_detection_diagnostic(records, tmp_path / 'unavailable_example.png', examples, example_prefix=prefix)
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('bad_signal', [np.nan, np.inf, -np.inf, -0.1, np.ma.masked,
                                     np.ma.array(5.0, mask=True), np.array([5.0]),
                                     np.complex128(5.0 + 0.2j), 5.0 + 0.0j])
@pytest.mark.parametrize('record_index', [0, 1])
def test_population_plot_rejects_invalid_signal_to_noise_before_rendering(tmp_path, monkeypatch, bad_signal, record_index):
    import warnings
    from unittest.mock import Mock
    from matplotlib import pyplot as plt
    from pcat.plotting import plot_detection_diagnostic

    render = Mock(side_effect=AssertionError('Invalid signal-to-noise reached figure creation'))
    monkeypatch.setattr(plt, 'subplots', render)
    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=5.0),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=0.0)]
    records[record_index]['injected_signal_to_noise'] = bad_signal
    with warnings.catch_warnings(), np.errstate(all='raise'), pytest.raises(ValueError, match='signal-to-noise'):
        warnings.simplefilter('error')
        plot_detection_diagnostic(records, tmp_path / 'bad_signal_to_noise.png')
    render.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('signals', [
    (5, 0), (5.25, 0.0), ('5.25', '0.0'), (np.float32(5.25), np.float32(0.0)),
    (np.array(5.25), np.array(0.0)),
    (np.ma.array(5.25, mask=False), np.ma.array(0.0, mask=False)),
])
@pytest.mark.parametrize('typefileplot', ['png', 'pdf'])
def test_population_plot_preserves_all_finite_signal_to_noise_points(tmp_path, monkeypatch, signals, typefileplot):
    import warnings
    from matplotlib import pyplot as plt
    from pcat import plotting

    records = [dict(has_subhalo=True, posterior_one=0.9, injected_index=0, recovered_index=0,
                    injected_signal_to_noise=signals[0]),
               dict(has_subhalo=False, posterior_one=0.1, injected_index=-1, recovered_index=0,
                    injected_signal_to_noise=signals[1])]
    originals = [record.copy() for record in records]
    captured = []
    save = plotting.save_figure
    before = set(plt.get_fignums())

    def capture(figure, path, **kwargs):
        axis = figure.axes[0]
        assert len(axis.collections) == 2
        for collection, index in zip(axis.collections, (1, 0)):
            offsets = collection.get_offsets()
            assert not np.ma.getmaskarray(offsets).any()
            np.testing.assert_array_equal(offsets, [[float(signals[index]), records[index]['posterior_one']]])
        assert axis.texts[0].get_text().startswith('TPR 100%')
        captured.append(figure)
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plotting, 'save_figure', capture)
    with warnings.catch_warnings(), np.errstate(all='raise'):
        warnings.simplefilter('error')
        plotting.plot_detection_diagnostic(records, tmp_path / f'finite_signal_to_noise.{typefileplot}')
    assert len(captured) == 1
    assert records == originals
    assert set(plt.get_fignums()) == before
