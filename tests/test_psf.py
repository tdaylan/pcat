from types import SimpleNamespace

import numpy as np
import pytest

from pcat.psf import CubicSubpixelPSF, psf_poly_fit


def test_cubic_subpixel_evaluation_matches_normalized_psf_and_deepcopy():
	from copy import deepcopy
	from pcat.main import PointSpreadFunction

	grid_x, grid_y = np.meshgrid(np.arange(-4., 5.), np.arange(-4., 5.), indexing='ij')
	state = SimpleNamespace(numbener=2, numbdqlt=1, typeexpr='test')
	profile = PointSpreadFunction(state, [1.5, 2.0], 'singgaus')
	model = CubicSubpixelPSF(profile, grid_x.ravel(), grid_y.ravel(), 1.)
	indices = np.arange(grid_x.size)
	for xpos, ypos in ((0., 0.), (0.37, -0.29), (-0.49, 0.49), (12., 13.), (1e20, -1e20)):
		expected = profile(np.hypot(grid_x.ravel() - xpos, grid_y.ravel() - ypos))
		actual = model.evaluate(xpos, ypos, indices)
		np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-4)
		assert np.all(actual >= 0)
		assert actual.shape == (2, grid_x.size, 1)
	assert deepcopy(model) is model
	np.testing.assert_array_equal(model(np.array([0., 1., 2.])), profile(np.array([0., 1., 2.])))
	np.testing.assert_allclose(model.evaluate(0.1, 0.2, indices[::3]),
	                           model.evaluate(0.1, 0.2, indices)[:, ::3, :])


def test_cubic_subpixel_rejects_irregular_coordinates():
	with pytest.raises(ValueError, match='lattice'):
		CubicSubpixelPSF(None, [0., 1.1], [0., 1.], 1.)


@pytest.mark.parametrize('width', [0.1, 0.5, 2.0])
def test_cubic_subpixel_preserves_flux_continuity_and_serialization(width):
	import pickle
	from pcat.main import PointSpreadFunction

	grid_x, grid_y = np.meshgrid(np.arange(-5., 6.), np.arange(-5., 6.), indexing='ij')
	profile = PointSpreadFunction(SimpleNamespace(numbener=1, numbdqlt=1, typeexpr='test'), [width], 'singgaus')
	model = CubicSubpixelPSF(profile, grid_x.ravel(), grid_y.ravel(), 1.)
	restored = pickle.loads(pickle.dumps(model))
	indices = np.arange(grid_x.size)
	peak = profile(np.array([0.])).max()
	for xpos in (0.499999, 0.500001, -0.499999, -0.500001):
		expected = profile(np.hypot(grid_x.ravel() - xpos, grid_y.ravel() - 0.17))
		actual = model.evaluate(xpos, 0.17, indices)
		np.testing.assert_allclose(actual, expected, rtol=0, atol=1.1e-3 * peak)
		assert abs(actual.sum() - expected.sum()) < 2e-3 * max(expected.sum(), peak)
		np.testing.assert_array_equal(restored.evaluate(xpos, 0.17, indices), actual)
	left = model.evaluate(0.499999, 0.17, indices)
	right = model.evaluate(0.500001, 0.17, indices)
	assert np.max(np.abs(left - right)) < 2e-3 * peak


@pytest.mark.parametrize('value', [0, np.nan, np.inf, True, [1.0], 1 + 1j])
def test_cubic_subpixel_rejects_invalid_pixel_scale(value):
	with pytest.raises(ValueError, match='positive finite scalars'):
		CubicSubpixelPSF(None, [0.], [0.], value)


def test_point_source_pipeline_selects_cubic_evaluation_without_distance_recomputation(monkeypatch):
	from pcat import main

	grid_x, grid_y = np.meshgrid(np.arange(-3., 4.), np.arange(-3., 4.), indexing='ij')
	state = SimpleNamespace(numbener=2, numbdqlt=1, typeexpr='test', typepixl='cart',
	                        xposgrid=grid_x.ravel(), yposgrid=grid_y.ravel(), sizepixl=1., kernevaltype='ulip')
	profile = main.PointSpreadFunction(state, [1.5, 2.0], 'singgaus')
	model = CubicSubpixelPSF(profile, state.xposgrid, state.yposgrid, 1.)
	indices = np.arange(grid_x.size)
	spectrum = np.array([2., 3.])
	expected = main.retr_sbrtpnts(state, 0.23, -0.17, spectrum, profile, indices)
	state.kernevaltype = 'bspx'
	monkeypatch.setattr(main, 'retr_angldistunit', lambda *args: pytest.fail('Cubic path recomputed direct angular distances'))
	actual = main.retr_sbrtpnts(state, 0.23, -0.17, spectrum, model, indices)
	np.testing.assert_allclose(actual, expected, rtol=0, atol=6e-4)
	state.typepixl = 'heal'
	with pytest.raises(ValueError, match='Cartesian'):
		main.retr_sbrtpnts(state, 0., 0., spectrum, model, indices)


def test_native_cartesian_sampler_runs_with_cubic_psf(tmp_path):
	from pathlib import Path
	import runpy
	from pcat.main import readfile
	from pcat.sampling import sample

	script = Path(__file__).parents[1] / 'examples' / 'daylan+2017_fermi_point_sources' / 'daylan+2017_fermi_point_sources.py'
	example = runpy.run_path(str(script))
	configuration = example['build_configuration'](smoke=True)
	template = example['build_ngpc_diffuse_template'](8)
	model = dict(configuration['dicttrue'],
	             sbrtbacknorm=[np.ones_like(template), template],
	             typemodlpsfn='singgaus', psfpexpr=np.full(6, 0.13))
	configuration.update(
		kernevaltype='bspx', numbsidecart=8, truenumbelempop0=1,
		truemaxmnumbelempop0=1, numbelempop0reg0=1,
		fittminmnumbelempop0=0, fittmaxmnumbelempop0=2,
		dicttrue=model, dictfitt=model, fermscalfact=np.ones((3, 2)),
		numbswep=12, numbburn=4, numbsamp=8, numbproc=1,
		boolcondcatl=False, boolmakeplot=False, boolmakeplotinit=False,
		boolmakeplotfram=False, boolmakeplotfinlpost=False, makeanim=False,
		pathbase=str(tmp_path), strgcnfg='cubic_psf_native_test',
	)
	run = sample(**configuration)
	posterior = readfile(run.pathoutpcnfg + 'gdatfinlpost')
	assert posterior.kernevaltype == 'bspx'
	assert np.isfinite(posterior.listpostlpostotl).all()


def test_psf_poly_fit_reconstructs_cubic_profile():
	factusam = 4
	indxusam = np.arange(12, dtype=np.float32)
	psfnusam = (12. - indxusam) * (1. + 0.1 * indxusam + 0.01 * indxusam**2)
	gdat = SimpleNamespace()

	coefspix = psf_poly_fit(gdat, psfnusam, factusam)
	offs = np.linspace(0., 1., factusam + 1, dtype=np.float32)
	design = np.column_stack((np.ones_like(offs), offs, offs**2, offs**3))
	reconstructed = (design @ coefspix).T.reshape(-1)[
		np.r_[0:5, 6:10, 11:15]
	]
	expected = np.pad(psfnusam, (0, 1))

	assert coefspix.shape == (4, 3)
	assert coefspix.dtype == np.float32
	assert np.allclose(reconstructed, expected, rtol=2e-5, atol=2e-5)
	assert np.array_equal(gdat.enerpsfnusam, np.arange(13) - 6.)


@pytest.mark.parametrize('factusam', [0, 2, 3.5])
def test_psf_poly_fit_rejects_invalid_oversampling(factusam):
	with pytest.raises(ValueError):
		psf_poly_fit(SimpleNamespace(), np.ones(12), factusam)


def test_psf_width_compatibility_wrapper_uses_shared_tdpy_interpolation():
	from pcat.main import retr_psfnwdth
	import tdpy

	angles = np.linspace(0.0, 4.0, 401)
	sigmas = np.array([[0.5, 1.0], [1.5, 2.0]])
	profiles = np.exp(-angles[None, :, None] ** 2 / (2.0 * sigmas[:, None, :] ** 2))
	state = SimpleNamespace(
		numbener=2,
		numbdqlt=2,
		indxener=np.arange(2),
		indxdqlt=np.arange(2),
		blimpara=SimpleNamespace(angl=angles),
	)

	np.testing.assert_array_equal(
		retr_psfnwdth(state, profiles, 0.5),
		tdpy.retr_psfnwdth(profiles, angles, 0.5),
	)


def test_bspx_point_spread_uses_cached_cubic_coefficients(monkeypatch):
	from pcat.main import PointSpreadFunction
	import tdpy.util

	parameters = np.array([0.5])
	exact_state = SimpleNamespace(numbener=1, numbdqlt=1, typeexpr='test', kernevaltype='ulip')
	interpolated_state = SimpleNamespace(
		numbener=1, numbdqlt=1, typeexpr='test', kernevaltype='bspx',
		typepixl='cart', sizepixl=0.2, maxmgangdata=2.0,
	)
	exact = PointSpreadFunction(exact_state, parameters, 'singgaus')
	interpolated = PointSpreadFunction(interpolated_state, parameters, 'singgaus')
	angles = np.linspace(0.0, 3.0, 121)
	expected = exact(angles)
	parameters[0] = 0.8
	assert interpolated.psfp[0] == 0.5

	def exact_evaluation_disabled(*args, **kwargs):
		raise RuntimeError('bspx evaluation called the analytic PSF')

	monkeypatch.setattr(tdpy.util, 'retr_singgaus', exact_evaluation_disabled)
	actual = interpolated(angles)

	assert actual.shape == expected.shape
	np.testing.assert_allclose(actual, expected, rtol=0.01, atol=1e-7)


def test_cubic_psf_benchmark_uses_actual_evaluators_and_writes_nonblank_figures(tmp_path):
	from pathlib import Path
	import runpy
	import matplotlib.image as mpimg

	script = Path(__file__).parents[1] / 'examples' / 'psf_cubic_interpolation' / 'psf_cubic_interpolation.py'
	example = runpy.run_path(str(script))
	example['run_benchmark'].__globals__['EXAMPLE_PATH'] = tmp_path
	result = example['run_benchmark'](evaluations=12, repeats=1)
	assert [record['energy_bins'] for record in result['records']] == [3, 16]
	for record in result['records']:
		assert record['direct_seconds'] > 0 and record['cubic_seconds'] > 0
		assert record['maximum_error_over_peak'] < 1e-3
	for path in (result['timing_figure'], result['error_figure']):
		assert path.parent == tmp_path / 'visuals'
		assert np.ptp(mpimg.imread(path)[..., :3]) > 0.2


def test_cubic_psf_benchmark_writes_interactive_bokeh_html(tmp_path):
	pytest.importorskip('bokeh')
	from pathlib import Path
	import runpy

	script = Path(__file__).parents[1] / 'examples' / 'psf_cubic_interpolation' / 'psf_cubic_interpolation.py'
	example = runpy.run_path(str(script))
	example['run_benchmark'].__globals__['EXAMPLE_PATH'] = tmp_path
	result = example['run_benchmark'](evaluations=12, repeats=1)
	output_path = example['run_bokeh_example'](result, tmp_path / 'cubic_psf_benchmark.html')

	assert output_path.is_file()
	html = output_path.read_text(encoding='utf-8')
	assert 'Synthetic double-King PSF evaluation cost' in html
	assert 'Held-out cubic PSF error' in html