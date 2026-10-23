import importlib.util
from pathlib import Path

import numpy as np

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "kl_divergence_gaussian"
    / "kl_divergence_gaussian.py"
)
SPECIFICATION = importlib.util.spec_from_file_location("kl_divergence_example", SCRIPT)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


def test_gaussian_posterior_kl_matches_closed_form_and_writes_figure(tmp_path, monkeypatch):
    monkeypatch.setattr(example, "EXAMPLE_PATH", tmp_path)
    result = example.run_example(sample_count=500, burn_count=250, chain_count=2)

    assert np.isfinite(result["kl_divergence"])
    assert result['chains'].shape[0] == 2 and result['chains'].shape[2] == 2
    assert 0 < result['chains'].shape[1] <= 500
    assert result["standard_error"] > 0.0
    assert abs(result["kl_divergence"] - result["analytic_kl"]) < max(
        5.0 * result["standard_error"], 0.2
    )
    assert result["figure_path"].is_file()
    information = result['mutual_information']
    assert np.isfinite(information['mutual_information'])
    assert abs(information['mutual_information'] - result['analytic_mutual_information']) < max(
        5 * information['standard_error'], 0.03,
    )
    import matplotlib.image as mpimg
    for path in (result['figure_path'], result['mutual_information_figure_path']):
        assert path.parent == tmp_path / 'visuals'
        assert np.ptp(mpimg.imread(path)[..., :3]) > 0.2