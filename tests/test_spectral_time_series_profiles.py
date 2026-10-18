import numpy as np
import pytest
from types import SimpleNamespace

from pcat.main import retr_spec as pcat_retr_spec
from pcat.spectral import (
    apply_gaussian_resolving_power,
    apply_line_spread_function,
    evaluate_line_profile,
    gaussian_lsf_kernel,
    spectral_profile_parameters,
)
from tdpy.spectral import (
    evaluate_continuum_spectrum,
    evaluate_line_profile as evaluate_line_profile_shared,
)
from pcat.time_series import (
    evaluate_flare_profile,
    evaluate_rotating_spot_profile,
    flare_profile_parameters,
    rotating_spot_parameters,
)
from pcat.demo import write_time_series_fits


def test_pcat_reexports_shared_tdpy_line_profile():
    assert evaluate_line_profile is evaluate_line_profile_shared


@pytest.mark.parametrize(
    ("profile", "parameters"),
    (
        ("powr", {}),
        ("curv", {"curv": np.array([0.1, 0.2])}),
        ("expc", {"expc": np.array([0.8, 1.2])}),
    ),
)
def test_pcat_continuum_models_use_tdpy_spectrum_kernel(profile, parameters):
    energy = np.linspace(0.1, 1.0, 21)
    flux = np.array([2.0, 3.0])
    spectral_index = np.array([1.0, 2.0])
    state = SimpleNamespace(
        numbener=energy.size,
        bctrpara=SimpleNamespace(ener=energy, enerplot=energy),
        enerpivt=0.5,
        factlogtenerpivt=np.log(energy / 0.5),
    )
    actual = pcat_retr_spec(
        state, flux, spectype=profile, sind=spectral_index, **parameters
    )
    expected = evaluate_continuum_spectrum(
        energy,
        flux,
        profile,
        0.5,
        spectral_index,
        curvature=parameters.get("curv"),
        cutoff_energy=parameters.get("expc"),
        log_energy_ratio=state.factlogtenerpivt,
    )

    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)


@pytest.mark.parametrize(
    ("profile", "parameters", "tolerance"),
    [
        ("gaus", {"gaussian_width": [1.0]}, 1e-6),
        ("lore", {"lorentz_width": [1.0]}, 2e-2),
        ("voig", {"gaussian_width": [1.0], "lorentz_width": [0.5]}, 2e-2),
        ("pvoi", {"gaussian_width": [1.0], "lorentz_width": [1.0], "mixing_fraction": [0.3]}, 2e-2),
        ("sinc", {"gaussian_width": [1.0]}, 3e-3),
        ("skew", {"gaussian_width": [1.0], "skewness": [4.0]}, 1e-6),
        ("toph", {"gaussian_width": [1.0]}, 2e-3),
    ],
)
def test_line_profiles_preserve_integrated_flux(profile, parameters, tolerance):
    axis = np.linspace(-50.0, 50.0, 100_001)
    values = evaluate_line_profile(axis, profile, [2.0], [0.0], **parameters)[:, 0]
    assert np.trapezoid(values, axis) == pytest.approx(2.0, rel=tolerance)


def test_skew_gaussian_is_asymmetric():
    axis = np.linspace(-8.0, 8.0, 10_001)
    values = evaluate_line_profile(axis, "skew", [1.0], [0.0], gaussian_width=[1.0], skewness=[5.0])[:, 0]
    assert np.trapezoid(values[axis > 0.0], axis[axis > 0.0]) > 10.0 * np.trapezoid(values[axis < 0.0], axis[axis < 0.0])


def test_fixed_and_resolving_power_lsf_preserve_sampled_flux():
    axis = np.linspace(1.0, 2.0, 2001)
    spectrum = evaluate_line_profile(axis, "gaus", [2.0], [1.5], gaussian_width=[0.01])
    fixed = apply_line_spread_function(spectrum, gaussian_lsf_kernel(3.0))
    resolving = apply_gaussian_resolving_power(axis, spectrum, [1.5], 500.0)
    sigma = 1.5 / (2.0 * np.sqrt(2.0 * np.log(2.0)) * 500.0 * np.median(np.diff(axis)))
    expected = apply_line_spread_function(spectrum, gaussian_lsf_kernel(sigma))
    np.testing.assert_allclose(resolving, expected, rtol=1e-13, atol=1e-15)
    assert np.sum(fixed) == pytest.approx(np.sum(spectrum))
    assert np.sum(resolving) == pytest.approx(np.sum(spectrum))
    assert np.max(resolving) < np.max(spectrum)


def test_tabulated_lsf_longer_than_spectrum_keeps_shape_and_flux():
    spectrum = np.array([0.0, 1.0, 0.0])
    convolved = apply_line_spread_function(spectrum, np.ones(9))
    assert convolved.shape == spectrum.shape
    assert np.sum(convolved) == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("profile", "parameters"),
    [
        ("flargauss", {"fwhm": [0.4]}),
        ("flarexpd", {"decay_time": [0.8]}),
        ("flarfred", {"rise_time": [0.2], "decay_time": [0.8]}),
        ("flardav", {"fwhm": [0.4]}),
    ],
)
def test_flare_profiles_are_nonnegative_and_peak_normalized(profile, parameters):
    time = np.linspace(-2.0, 4.0, 6001)
    values = evaluate_flare_profile(time, profile, [3.0], [0.0], **parameters)[:, 0]
    assert np.all(values >= 0.0)
    assert values[np.argmin(np.abs(time))] == pytest.approx(3.0, rel=0.02)


@pytest.mark.parametrize("profile", ["flare", "spot"])
def test_gaussian_transients_preserve_original_widths(profile):
    time = np.array([-0.4, -0.2, 0.0, 0.2, 0.4])  # [day]
    widths = np.array([0.4, 0.8])  # [day]
    amplitudes = np.array([3.0, 5.0])  # [count]
    sigma = widths / (2.0 * np.sqrt(2.0 * np.log(2.0)))  # [day]
    expected = amplitudes * np.exp(-0.5 * (time[:, None] / sigma)**2)
    if profile == "flare":
        actual = evaluate_flare_profile(time, "flargauss", amplitudes, [0.0, 0.0], fwhm=widths)
    else:
        actual = evaluate_rotating_spot_profile(time, amplitudes, [0.0, 0.0], widths, period_days=5.0)
        expected = -expected
    np.testing.assert_allclose(actual, expected, rtol=1e-13)


def test_fred_has_independent_rise_and_decay_times():
    values = evaluate_flare_profile(
        [-1.0, 0.0, 1.0], "flarfred", [1.0], [0.0], rise_time=[0.2], decay_time=[2.0]
    )[:, 0]
    assert values[0] < values[2] < values[1]


def test_rotating_spot_profile_repeats_and_only_dims_the_star():
    time_days = np.array([1.25, 4.45, 7.65, 2.85])  # [day]
    values = evaluate_rotating_spot_profile(
        time_days,
        [120.0],
        [0.25],
        [0.16],
        period_days=3.2,  # [day]
        reference_time_days=1.0,  # [day]
    )[:, 0]

    np.testing.assert_allclose(values[:3], -120.0)
    assert values[3] > -1e-8
    assert np.all(values <= 0.0)
    assert rotating_spot_parameters("spotrot") == ("flux", "elin", "fwhm")


def test_rotating_spot_profile_rejects_width_longer_than_period():
    with pytest.raises(ValueError, match="shorter than the rotation period"):
        evaluate_rotating_spot_profile(
            [0.0], [10.0], [0.0], [2.0], period_days=1.0  # [day]
        )


def test_rotating_spot_profile_supports_zero_spot_state():
    values = evaluate_rotating_spot_profile(
        [0.0, 1.0], [], [], [], period_days=2.0  # [day]
    )

    assert values.shape == (2, 0)


def test_native_parameter_sets_are_profile_specific():
    assert spectral_profile_parameters("pvoi") == ("flux", "elin", "sigm", "gamm", "frac")
    assert spectral_profile_parameters("skew") == ("flux", "elin", "sigm", "skew")
    assert flare_profile_parameters("flarfred") == ("flux", "elin", "scalrise", "scalfall")


@pytest.mark.parametrize("baseline_counts", [5.0, [5.0, 7.0], np.array([[[5.0]], [[7.0]]])])
def test_time_series_fits_use_each_bin_width_and_return_baseline(tmp_path, baseline_counts):
    from astropy.io import fits

    edges = np.array([1.0, 1.1, 1.4])
    counts = np.array([12.0, 24.0])
    baseline = write_time_series_fits(tmp_path, edges, counts, baseline_counts)

    np.testing.assert_array_equal(fits.getdata(tmp_path / "sbrt.fits").ravel(), counts)
    np.testing.assert_allclose(fits.getdata(tmp_path / "expo.fits").ravel(), [10.0, 1.0 / 0.3])
    expected = [5.0, 5.0] if np.ndim(baseline_counts) == 0 else np.asarray(baseline_counts).ravel()
    np.testing.assert_array_equal(baseline[:, 0, 0], expected)
    assert baseline.shape == (2, 1, 1)
    if isinstance(baseline_counts, np.ndarray):
        assert not np.shares_memory(baseline, baseline_counts)


def test_time_series_fits_delegate_overwrite_and_narration(tmp_path, monkeypatch, capsys):
    from astropy.io import fits
    from unittest.mock import Mock
    import tdpy.serialization

    monkeypatch.setenv('TDPY_VERBOSITY', '1')
    writer = Mock(wraps=tdpy.serialization.write_fits)
    monkeypatch.setattr(tdpy.serialization, 'write_fits', writer)
    directory = tmp_path / 'nested'
    edges = np.array([1.0, 1.1, 1.4])
    for counts in (np.array([12.0, 24.0]), np.array([30.0, 40.0])):
        baseline = write_time_series_fits(directory, edges, counts, 5.0)

        assert capsys.readouterr().out == ''.join(
            f'Writing to {directory / name}...\n' for name in ('sbrt.fits', 'expo.fits')
        )
        print(f'Reading from {directory / "sbrt.fits"}...')
        brightness = fits.getdata(directory / 'sbrt.fits')
        print(f'Reading from {directory / "expo.fits"}...')
        exposure = fits.getdata(directory / 'expo.fits')
        assert brightness.shape == (2, 1, 1, 1)
        assert exposure.shape == (2, 1, 1)
        np.testing.assert_array_equal(brightness[:, 0, 0, 0], counts)
        np.testing.assert_allclose(exposure[:, 0, 0], 1.0 / np.diff(edges))
        np.testing.assert_array_equal(baseline[:, 0, 0], [5.0, 5.0])
        capsys.readouterr()

    assert writer.call_count == 4
    assert all(call.kwargs == {'overwrite': True} for call in writer.call_args_list)


def test_time_series_fits_reject_nonincreasing_bin_edges(tmp_path):
    with pytest.raises(ValueError, match="strictly increasing"):
        write_time_series_fits(tmp_path, [1.0, 1.2, 1.2], [3.0, 4.0], 2.0)


@pytest.mark.parametrize("baseline", [[1.0, 2.0, 3.0], np.nan, [1.0, np.inf]])
def test_time_series_fits_validate_baseline_before_writing(tmp_path, baseline):
    with pytest.raises(ValueError, match="baseline counts"):
        write_time_series_fits(tmp_path, [1.0, 1.2, 1.4], [3.0, 4.0], baseline)
    assert not any(tmp_path.iterdir())