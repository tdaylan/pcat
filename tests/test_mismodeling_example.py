import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (Path(__file__).resolve().parents[1] / "examples" / "mismodeling_flare_catalog"
          / "mismodeling_flare_catalog.py")
SPECIFICATION = importlib.util.spec_from_file_location("mismodeling_flare_catalog", SCRIPT)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


def test_each_light_curve_has_a_correct_and_a_misspecified_fit():
    light_curves = example.simulate_light_curves()
    assert light_curves["counts"]["A"].shape == light_curves["counts"]["B"].shape == light_curves["time"].shape
    # light curve B differs from A only through the rotational modulation of the quiescent level
    modulation = light_curves["baselines"]["modulated"] / light_curves["baselines"]["constant"] - 1.0
    np.testing.assert_allclose(np.abs(modulation).max(), example.MODULATION_AMPLITUDE, rtol=1e-3)
    for light_curve in ("A", "B"):
        fits = [row for row in example.FITS if row[1] == light_curve]
        assert len(fits) == 2
        assert sum("misspecified" in row[4] for row in fits) == 1
        assert {(row[2], row[3]) for row in fits} != {fits[0][2:4]}
    # run names differ from the folder name, so cache cleanup cannot delete the example folder
    assert all(row[0] != SCRIPT.parent.name for row in example.FITS)


@pytest.mark.parametrize("baseline_kind", ["constant", "modulated"])
def test_example_reuses_shared_fits_writer_and_preserves_baseline_template(tmp_path, monkeypatch, baseline_kind):
    from unittest.mock import Mock
    from astropy.io import fits
    from pcat import sampling
    from pcat.demo import write_time_series_fits

    assert example.write_time_series_fits is write_time_series_fits
    writer = Mock(wraps=write_time_series_fits)
    monkeypatch.setattr(example, "write_time_series_fits", writer)
    monkeypatch.setattr(example, "EXAMPLE_PATH", tmp_path)
    captured = {}
    def sample(**configuration):
        captured.update(configuration)
    monkeypatch.setattr(sampling, "sample", sample)
    posterior = object()
    monkeypatch.setattr(example, "readfile", lambda path: posterior)
    light_curves = example.simulate_light_curves()
    result = example.run_fit(light_curves, "writer_test", "B", "flarfred", baseline_kind, numbswep=300)
    assert result is posterior
    writer.assert_called_once()
    input_path = Path(example.retr_pathrun(tmp_path, "writer_test")) / "data" / "inpt"
    observed = fits.getdata(input_path / "sbrt.fits")
    exposure = fits.getdata(input_path / "expo.fits")
    np.testing.assert_array_equal(observed[:, 0, 0, 0], light_curves["counts"]["B"])
    np.testing.assert_allclose(exposure[:, 0, 0] * np.diff(light_curves["edges"]), 1.0)
    np.testing.assert_array_equal(
        captured["dictfitt"]["sbrtbacknorm"][0][:, 0, 0], light_curves["baselines"][baseline_kind],
    )


def test_misspecified_profile_fit_samples_gaussian_flares(tmp_path, monkeypatch):
    monkeypatch.setattr(example, "EXAMPLE_PATH", tmp_path)
    monkeypatch.setattr(example, "NUMBER_CHAINS", 1)
    light_curves = example.simulate_light_curves()
    posterior = example.run_fit(light_curves, "profile_wrong_test", "A", "flargauss", "constant", numbswep=300)
    catalog = posterior.listpostdictelem[-1][0]
    assert "fwhm" in catalog and "scalrise" not in catalog
    number = np.asarray(posterior.listpostnumbelem).ravel()
    assert np.all((number >= 0) & (number <= example.MAXIMUM_NUMBER))
    model = np.asarray(posterior.listpostcntpmodl)[:, :, 0, 0]
    assert model.shape[1] == light_curves["time"].size and np.isfinite(model).all()


def test_residual_and_predictive_chi_squared_shape_diagnostics(tmp_path, monkeypatch):
    rng = np.random.default_rng(41)
    model = np.full(500, 1000.0)
    counts = rng.poisson(model)
    posterior_models = np.tile(model, (16, 1))
    diagnostics = example.calculate_residual_diagnostics(
        counts, posterior_models, np.random.default_rng(52), number_replicates=300
    )

    assert abs(diagnostics["residual"].mean()) < 0.2
    assert abs(diagnostics["residual"].std(ddof=1) - 1.0) < 0.1
    assert diagnostics["chi2_predictive"].shape == (300,)
    assert abs(diagnostics["chi2_predictive"].mean() - counts.size) < 20.0
    assert 0.0 <= diagnostics["chi2_upper_tail_fraction"] <= 1.0

    monkeypatch.setattr(example, "VISUAL_PATH", tmp_path)
    summary = {row[0]: diagnostics for row in example.FITS}
    figure_path = example.plot_residual_chi2_shapes(summary, "png")
    assert figure_path.is_file()
    assert figure_path.stat().st_size > 10_000
