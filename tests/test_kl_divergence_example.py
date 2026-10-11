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
    assert result["standard_error"] > 0.0
    assert abs(result["kl_divergence"] - result["analytic_kl"]) < max(
        5.0 * result["standard_error"], 0.2
    )
    assert result["figure_path"].is_file()