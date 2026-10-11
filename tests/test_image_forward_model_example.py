import importlib.util
from pathlib import Path

import numpy as np
from PIL import Image

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "examples"
    / "image_forward_model"
    / "image_forward_model.py"
)
SPECIFICATION = importlib.util.spec_from_file_location("image_forward_model_example", SCRIPT)
example = importlib.util.module_from_spec(SPECIFICATION)
SPECIFICATION.loader.exec_module(example)


def test_pcat_image_forward_model_example_preserves_scene_and_writes_figure(tmp_path, monkeypatch):
    monkeypatch.setattr(example, "EXAMPLE_PATH", tmp_path)
    source_image = example.create_simulated_scene()
    np.testing.assert_allclose(source_image.sum(), 1650.0)

    products, figure_path = example.run_example()

    assert products["source_image"].shape == (51, 51)
    np.testing.assert_allclose(products["psf_kernel"].sum(), 1.0)
    np.testing.assert_allclose(
        products["observed_image"].sum(),
        source_image.sum() + example.BACKGROUND_RATE * source_image.size,
    )
    assert figure_path.is_file()
    with Image.open(figure_path) as image:
        assert image.width > 1000 and image.height > 300
