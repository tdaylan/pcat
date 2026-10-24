"""Measure cubic subpixel PSF error and cost against direct double-King profiles."""

import argparse
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
from tdpy.cli import add_plot_arguments
from tdpy.verbosity import print

from pcat.demo import configure_example_style, save_example_figure
from pcat.main import PointSpreadFunction
from pcat.psf import CubicSubpixelPSF

EXAMPLE_PATH = Path(__file__).resolve().parent


def run_benchmark(evaluations=600, repeats=3, typefileplot='png'):
    """Benchmark fixed synthetic PSFs without a posterior or survey claim."""
    if evaluations < 2 or repeats < 1:
        raise ValueError('Require at least two evaluations and one timing repeat')
    configure_example_style('white')
    pixel_size = np.deg2rad(0.01) # [rad/pixel]
    grid_x, grid_y = np.meshgrid(np.arange(-10., 11.) * pixel_size,
                                np.arange(-10., 11.) * pixel_size, indexing='ij') # [rad]
    indices = np.flatnonzero((abs(grid_x.ravel()) <= 4 * pixel_size) & (abs(grid_y.ravel()) <= 4 * pixel_size))
    phases = np.random.default_rng(62).uniform(-0.49, 0.49, size=(evaluations, 2)) * pixel_size # [rad]
    records = []
    for energy_bins in (3, 16):
        state = SimpleNamespace(numbener=energy_bins, numbdqlt=2, typeexpr='test')
        parameters = np.tile([2 * pixel_size, 2.4, 4 * pixel_size, 2., 0.82], energy_bins * 2)
        profile = PointSpreadFunction(state, parameters, 'doubking')
        started = perf_counter() # [s]
        fitted = CubicSubpixelPSF(profile, grid_x.ravel(), grid_y.ravel(), pixel_size)
        setup_seconds = perf_counter() - started # [s]
        timings = []
        for repeat in range(repeats):
            started = perf_counter() # [s]
            direct = [profile(np.hypot(grid_x.ravel()[indices] - x, grid_y.ravel()[indices] - y)) for x, y in phases]
            direct_seconds = perf_counter() - started # [s]
            started = perf_counter() # [s]
            cubic = [fitted.evaluate(x, y, indices) for x, y in phases]
            cubic_seconds = perf_counter() - started # [s]
            timings.append((direct_seconds, cubic_seconds))
        peak = float(profile(np.array([0.])).max())
        errors = np.array([np.max(np.abs(first - second)) / peak for first, second in zip(direct, cubic)])
        direct_seconds, cubic_seconds = np.median(timings, axis=0) # [s]
        records.append({
            'energy_bins': energy_bins, 'setup_seconds': float(setup_seconds),
            'direct_seconds': float(direct_seconds), 'cubic_seconds': float(cubic_seconds),
            'evaluation_speed_ratio': float(direct_seconds / cubic_seconds),
            'maximum_error_over_peak': float(errors.max()), 'errors_over_peak': errors,
        })
        print(f'{energy_bins} energy bins, two quality classes: setup {setup_seconds:.3f} s, '
              f'evaluation ratio {direct_seconds / cubic_seconds:.2f}, '
              f'maximum absolute error / peak {errors.max():.3g}')

    figure, axis = plt.subplots(figsize=(7, 4.5))
    positions = np.arange(len(records))
    width = 0.24
    for offset, key, color, label in (
        (-width, 'setup_seconds', '#89939a', 'One-time cubic setup'),
        (0, 'direct_seconds', '#a51c30', f'{evaluations} direct evaluations'),
        (width, 'cubic_seconds', '#007360', f'{evaluations} cubic evaluations'),
    ):
        axis.bar(positions + offset, [record[key] for record in records], width,
                 color=color, label=label)
    axis.set_xticks(positions, [f'{record["energy_bins"]} energy bins' for record in records])
    axis.set(ylabel='Median evaluation or setup time [s]', xlabel='Two quality classes, 81-pixel cutout')
    axis.legend(loc='lower left', bbox_to_anchor=(0, 1.02), framealpha=1, fancybox=True)
    axis.grid(False)
    figure.tight_layout()
    timing_path = save_example_figure(figure, EXAMPLE_PATH, 'cubic_psf_cost', typefileplot)

    figure, axis = plt.subplots(figsize=(7, 4.0))
    for record, color in zip(records, ('#a51c30', '#007360')):
        axis.plot(np.arange(evaluations), record['errors_over_peak'], color=color,
                  label=f'{record["energy_bins"]} energy bins')
    axis.set(xlabel='Held-out random subpixel shift', ylabel='Maximum absolute density error / profile peak')
    axis.legend(loc='lower left', bbox_to_anchor=(0, 1.02), framealpha=1, fancybox=True)
    axis.grid(False)
    figure.tight_layout()
    error_path = save_example_figure(figure, EXAMPLE_PATH, 'cubic_psf_accuracy', typefileplot)
    return {'records': records, 'timing_figure': timing_path, 'error_figure': error_path}


def run_bokeh_example(result, output_path):
    """Write interactive timing and held-out error views from the same run."""
    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource, FixedTicker, HoverTool
    from bokeh.plotting import figure

    from tdpy.plotting import save_bokeh_html

    colors = {'setup_seconds': '#6C7373', 'direct_seconds': '#A51417', 'cubic_seconds': '#007360'}
    labels = {
        'setup_seconds': 'One-time cubic setup',
        'direct_seconds': 'Direct evaluation',
        'cubic_seconds': 'Cubic evaluation',
    }
    categories = [str(record['energy_bins']) for record in result['records']]
    timing = figure(
        x_range=(-0.6, len(categories) - 0.4), width=760, height=340,
        title='Synthetic double-King PSF evaluation cost',
        tools='pan,wheel_zoom,box_zoom,reset,save',
    )
    timing.add_tools(HoverTool(tooltips=[('Method', '@method'), ('Time [s]', '@seconds{0.0000}')]))
    for method, offset in zip(colors, (-0.22, 0.0, 0.22)):
        source = ColumnDataSource(data={
            'x': np.arange(len(categories), dtype=float) + offset,
            'seconds': [record[method] for record in result['records']],
            'method': [labels[method]] * len(categories),
        })
        timing.vbar(
            x='x', top='seconds', source=source,
            width=0.18, color=colors[method], legend_label=labels[method],
        )
    timing.xaxis.ticker = FixedTicker(ticks=list(range(len(categories))))
    timing.xaxis.major_label_overrides = {index: category for index, category in enumerate(categories)}
    timing.xaxis.axis_label = 'Energy bins, with two quality classes'
    timing.yaxis.axis_label = 'Median setup or evaluation time [s]'
    timing.legend.location = 'top_left'
    timing.legend.click_policy = 'hide'
    timing.grid.visible = False

    accuracy = figure(
        width=760, height=340,
        title='Held-out cubic PSF error',
        tools='pan,wheel_zoom,box_zoom,reset,save',
    )
    accuracy.add_tools(HoverTool(tooltips=[('Shift', '$x'), ('Error / peak', '$y{0.000000}')]))
    for record, color in zip(result['records'], ('#A51417', '#007360')):
        accuracy.line(
            np.arange(len(record['errors_over_peak'])), record['errors_over_peak'],
            line_width=2, color=color, legend_label=f"{record['energy_bins']} energy bins",
        )
    accuracy.xaxis.axis_label = 'Held-out random subpixel shift'
    accuracy.yaxis.axis_label = 'Maximum absolute density error / profile peak'
    accuracy.legend.location = 'top_left'
    accuracy.legend.click_policy = 'hide'
    accuracy.grid.visible = False
    return save_bokeh_html(
        column(timing, accuracy, sizing_mode='stretch_width'),
        output_path,
        title='PCAT cubic subpixel PSF benchmark',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_plot_arguments(parser)
    parser.add_argument('--evaluations', type=int, default=600)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--bokeh', action='store_true', help='also write an interactive HTML version')
    arguments = parser.parse_args()
    benchmark = run_benchmark(arguments.evaluations, arguments.repeats, arguments.typefileplot)
    if arguments.bokeh:
        run_bokeh_example(benchmark, EXAMPLE_PATH / 'visuals' / 'cubic_psf_benchmark.html')