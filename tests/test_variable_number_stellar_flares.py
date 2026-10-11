import astropy.io.fits
import h5py
import importlib.util
import numpy as np
from pathlib import Path
import sys
from types import ModuleType

from pcat import sampling
from pcat.time_series import FLARE_PROFILE_PARAMETERS


def _load_flare_example(monkeypatch):
    nicomedia = ModuleType("nicomedia")
    nicomedia.retr_lcurmodl_flarsing = lambda time, peak, amplitude, rise, fall: np.zeros_like(time)
    monkeypatch.setitem(sys.modules, "nicomedia", nicomedia)
    script = Path(__file__).resolve().parents[1] / "examples" / "variable_number_stellar_flares" / "variable_number_stellar_flares.py"
    specification = importlib.util.spec_from_file_location("pcat_flare_profile_example", script)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_variable_number_stellar_flares_uses_input_counts_bins_and_baseline_template(tmp_path):
    # a flat quiescent baseline of 1000 counts per bin with one injected burst, on time bins [day]
    edges = 1.0 + np.linspace(0.0, 1.0, 41)
    width = np.diff(edges)
    centers = 0.5 * (edges[1:] + edges[:-1])
    counts = np.round(1000. + 800. * np.exp(-0.5 * ((centers - 1.5) / 0.02)**2))
    exposure = 1. / width
    path_inpt = tmp_path / 'variable_number_stellar_flares' / 'data' / 'inpt'
    path_inpt.mkdir(parents=True)
    astropy.io.fits.writeto(path_inpt / 'sbrt.fits', counts[:, None, None, None].astype(float))
    astropy.io.fits.writeto(path_inpt / 'expo.fits', exposure[:, None, None])

    sampling.sample(
        typeexpr='fire', typedata='inpt', strgexprsbrt='sbrt.fits', typeexpo='file', strgexpo='expo.fits',
        binsenerfull=edges, spectype=['voig'], spatdisttype=['line'], typeelem=['lghtlinevoig'],
        dictfitt={'typeelem': ['lghtlinevoig'], 'spectype': ['voig'],
                  'sbrtbacknorm': [np.full((edges.size - 1, 1, 1), 1000.)], 'listnamediff': ['back0000']},
        limtparaelem={'flux': (1e1, 1e4), 'sigm': (5e-3, 5e-2), 'gamm': (1e-3, 2e-2)},
        maxmgangdata=1e-4, anlytype='spec', fittminmnumbelempop0=0, fittmaxmnumbelempop0=3,
        inittype='rand', numbswep=400, numbsamp=40, boolmakeplot=False, boolmakeplotinit=False,
        booldiag=False, typeverb=0, pathbase=str(tmp_path / 'variable_number_stellar_flares'),
        strgcnfg='variable_number_stellar_flares',
    )

    path_outp = tmp_path / 'variable_number_stellar_flares' / 'data' / 'outp' / 'variable_number_stellar_flares'
    with h5py.File(path_outp / 'gdatinit.h5', 'r') as file:
        np.testing.assert_allclose(file['cntpdata'][()].ravel(), counts)
    with h5py.File(path_outp / 'gdatfinlpost.h5', 'r') as file:
        background = file['listpostcntpback0000'][()]
        llik = file['listpostlliktotl'][()]
    # the baseline template enters the model with its input normalization
    np.testing.assert_allclose(background[-1].ravel(), 1000., rtol=1e-6)
    assert np.all(np.isfinite(llik))


def test_flare_example_builds_profile_specific_priors_proposals_and_templates(monkeypatch, tmp_path):
    example = _load_flare_example(monkeypatch)
    configurations = []
    monkeypatch.setattr(sampling, "sample", lambda **configuration: configurations.append(configuration))
    monkeypatch.setattr(example, "EXAMPLE_PATH", tmp_path)
    monkeypatch.setattr(example, "RUN_NAME", "flare_profile_options")
    edges = 1.0 + np.linspace(0.0, 0.2, 41)
    observed = np.full(edges.size - 1, example.BASELINE_COUNT_RATE)
    observed[15:20] += np.array([10.0, 80.0, 200.0, 80.0, 10.0])

    for profile, parameters in FLARE_PROFILE_PARAMETERS.items():
        limits, scales = example.flare_profile_configuration(profile)
        proposal = example.build_birth_proposal(edges, observed, profile)
        assert set(limits) == {"flux", *parameters[2:]}
        assert len(scales) == len(parameters)
        assert proposal["numbparapropelem"] == len(parameters)
        assert np.all(np.isfinite(proposal["pdfnposipropelem"]))

        example.run_pcat(edges, np.ones(edges.size - 1), 100, proposal, profile)
        configuration = configurations[-1]
        assert configuration["spectype"] == [profile]
        assert configuration["dictfitt"]["spectype"] == [profile]
        assert configuration["limtparaelem"] == limits
        assert configuration["stdvpropelemfire"] == scales
