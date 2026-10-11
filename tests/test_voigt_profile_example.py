import importlib.util
from pathlib import Path

import pytest


EXAMPLE_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "voigt_spectral_line_catalog"
    / "voigt_spectral_line_catalog.py"
)
SPECIFICATION = importlib.util.spec_from_file_location(
    "pcat_voigt_spectral_line_catalog", EXAMPLE_SCRIPT
)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


def test_voigt_profile_configuration_calls_pcat_pipeline(monkeypatch):
    calls = []

    def fake_sample(**configuration):
        calls.append(configuration)
        return {"configuration": configuration["strgcnfg"]}

    monkeypatch.setattr(example.sampling, "sample", fake_sample)

    result = example.run_voigt_profile_detection("nomi")

    configuration = calls[0]
    assert result == {"configuration": "voigt_nomi"}
    assert configuration["strgcnfg"] == "voigt_nomi"
    assert configuration["spectype"] == ["voig"]
    assert configuration["typeelem"] == ["lghtlinevoig"]
    assert configuration["strgexpo"] == pytest.approx(1.0e4)
    assert configuration["truenumbelempop0"] == 2
    assert configuration["typeseed"] == 0
    assert configuration["typeseedelem"] == 17
    assert configuration["inittype"] == "rand"
    assert configuration["numbframanim"] == 24
    assert configuration["dicttrue"]["typeelem"] == ["lghtlinevoig"]
    assert configuration["dictfitt"]["spectype"] == ["voig"]


def test_spectral_example_exposes_each_profile_specific_parameter_set():
    from pcat.spectral import spectral_profile_parameters

    for profile in example.SPECTRAL_EXAMPLE_PROFILES:
        common, _, _ = example.build_configurations(spectral_profile=profile)
        assert common["spectype"] == [profile]
        assert common["dictfitt"]["spectype"] == [profile]
        assert len(common["stdvpropelemfire"]) == len(spectral_profile_parameters(profile))
        assert common["limtparaelem"]["flux"] == (0.2, 10.0)
    assert example.build_configurations(spectral_profile="pvoi")[0]["limtparaelem"]["frac"] == (0.0, 1.0)
    assert example.build_configurations(spectral_profile="skew")[0]["limtparaelem"]["skew"] == (-10.0, 10.0)


def test_spectral_example_exposes_each_supported_instrument_lsf():
    gaussian = example.build_configurations(lsftype="gaus", lsfresolvingpower=4500.0)[0]
    tabulated = example.build_configurations(lsftype="tabu", lsfkernel_sigma=1.5)[0]

    assert gaussian["lsftype"] == "gaus"
    assert gaussian["lsfresolvingpower"] == pytest.approx(4500.0)
    assert tabulated["lsftype"] == "tabu"
    assert tabulated["lsfkernel"].sum() == pytest.approx(1.0)


def test_energy_dispersion_profile_requires_an_instrument_response():
    with pytest.raises(ValueError, match="requires an instrument response"):
        example.build_configurations(spectral_profile="edis")

    response = lambda centers: np.full_like(centers, 0.01)
    common, _, _ = example.build_configurations(spectral_profile="edis", edisintp=response)
    assert common["edisintp"] is response
