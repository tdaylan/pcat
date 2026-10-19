import numpy as np
import pytest

from pcat import plot_population_grid


@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('inherit', [False, True])
def test_population_grid_writes_corner_histogram_and_pair_plots(tmp_path, monkeypatch, typeplotback, inherit):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    from PIL import Image
    from tdpy.plotting import plot_background_colors, plot_style

    rng = np.random.default_rng(1)
    listpara = [rng.standard_normal((50, 3)), rng.standard_normal((200, 3)) + 3.0]
    background, foreground = plot_background_colors(typeplotback)
    save = plt.Figure.savefig

    def capture_and_save(figure, path, **kwargs):
        assert figure.get_facecolor() == to_rgba(background)
        for axis in figure.axes:
            assert axis.get_facecolor() == to_rgba(background)
            assert axis.xaxis.label.get_color() == foreground
            assert axis.yaxis.label.get_color() == foreground
        return save(figure, path, **kwargs)

    monkeypatch.setattr(plt.Figure, 'savefig', capture_and_save)
    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original = dict(plt.rcParams)
        plot_population_grid(
            [['a', ''], ['b', ''], ['c', '']], listpara=listpara, listlablpopl=['one', 'two'],
            pathbase=f'{tmp_path}/', strgextn='test', boolplottria=True, boolplothistodim=True,
            boolplotpair=True, typeverb=0, **options,
        )
        assert dict(plt.rcParams) == original

    names = {path.name for path in tmp_path.glob('*.png')}
    extension = 'test_dark' if typeplotback == 'dark' else 'test'
    assert f'pmar_{extension}_scatscat.png' in names
    assert {f'hist_p00{k}_{extension}.png' for k in range(3)} <= names
    assert any(name.startswith(f'pmar_p001_p000_{extension}') for name in names)
    expected = tuple(int(channel * 255) for channel in to_rgba(background)[:3])
    for name in names:
        with Image.open(tmp_path / name) as image:
            pixels = np.asarray(image.convert('RGB'))
            assert tuple(pixels[0, 0]) == expected
            assert np.any(pixels != expected)
