from pathlib import Path
import numpy as np
import xarray as xr

import arviz as az
import arviz_plots as azp
from arviz_stats import summary
from arviz_base.labels import MapLabeller, NoVarLabeller

from lime import label_decomposition
from matplotlib import pyplot as plt, rc_context
from matplotlib.colors import to_hex
from specsy.io import SpecSyError, specsy_cfg
from lime.plotting.plots import save_close_fig_swicth
from itertools import chain
from specsy.plotting.bokeh_functions import update_bokeh_figure
from specsy.plotting.plots import theme
from lime import Spectrum

try:
    from bokeh.plotting import figure, output_file, save, show
    from bokeh.layouts import gridplot
    from bokeh.models import GlobalInlineStyleSheet, Span, BoxAnnotation
    from bokeh.io import curdoc
    bokeh_check = True
except ImportError:
    bokeh_check = False


# Sentinel object for non input figures
_NO_FIG = object()


def ref_values_dataset(true_values, var_names):
    comp = {var: xr.DataArray(true_values[var]) for var in var_names if var in true_values}
    return xr.Dataset(comp).expand_dims(column=['dist'])


def latex_labeller(notation_map, backend):
    # Matplotlib mathtext wants $...$, bokeh wants $$...$$
    delim = '$$' if backend == 'bokeh' else '$'
    var_map = {var: f'{delim}{latex.strip("$")}{delim}' for var, latex in notation_map.items()}
    return MapLabeller(var_name_map=var_map)


def plot_fitted_fluxes(trace_data, output_address=None, backend='matplotlib',
                       n_cols=5, n_rows=None, col_row_scale=(10, 4), in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None,
                       display_check=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is _NO_FIG else False

    # Load the inference data if necessary
    if isinstance(trace_data, (str, Path)):
        trace_data = az.from_netcdf(trace_data)

    # Unpack the plot data
    labels_arr = trace_data.observed_data['lines'].values
    obs_flux_arr = trace_data.constant_data['input_flux'].values
    obs_err_arr = trace_data.constant_data['input_err'].values

    # Generate the observed fluxes containers
    band = np.column_stack([obs_flux_arr - obs_err_arr, obs_flux_arr + obs_err_arr])
    band = xr.Dataset({'theo_flux': xr.DataArray(band, dims=('lines', 'band'), coords={'lines': labels_arr})})
    ref_value = xr.Dataset({'theo_flux': xr.DataArray(obs_flux_arr, dims='lines', coords={'lines': labels_arr})})

    # Generate colors array for each species
    ion_arr, latex_arr = label_decomposition(labels_arr, params_list=['particle', 'latex_label'])
    unique_ions = np.unique(ion_arr)
    ion_cmap = {ion: to_hex(plt.cm.viridis(x)) for ion, x in zip(unique_ions, np.linspace(0, 0.9, unique_ions.size))}
    color_arr = [ion_cmap[ion] for ion in ion_arr]

    # Guess grid size
    n_rows = int(np.ceil(float(labels_arr.size)/float(n_cols)))
    n_cells = n_rows * n_cols

    # Plot configuration
    visuals = {'point_estimate': False, 'point_estimate_text': False, 'credible_interval': False, 'face': {'alpha': 0.3}}

    match backend:
        case 'matplotlib':

            # Set the plot format where the user's overwrites the default
            size_conf = {'figure.figsize': (n_cols, n_rows)}
            size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}
            plot_cfg = theme.fig_defaults(size_conf, fig_type='flux_grid')

            with rc_context(plot_cfg):
                func = azp.plot_dist(trace_data, var_names=['theo_flux'], kind='hist', labeller=NoVarLabeller(), color=color_arr,
                                     aes={'color': ['lines']}, visuals=visuals)
                azp.add_lines(func, values=ref_value, color=['black'])
                azp.add_bands(func, values=band, ref_dim=['band'])

                in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

        case 'bokeh':
            if bokeh_check:

                plot_cfg = theme.fig_defaults(user_fig=fig_cfg, fig_type='flux_grid', plot_lib='bokeh')

                in_fig = azp.plot_dist(trace_data, var_names=['theo_flux'], kind='hist', labeller=NoVarLabeller(),
                                      color=color_arr, backend='bokeh', aes={'color': ['lines']}, visuals=visuals)
                azp.add_lines(in_fig, values=ref_value, color=[theme.colors['fg']])
                azp.add_bands(in_fig, values=band, ref_dim=['band'], visuals={'ref_band': {'alpha': 0.2,
                                                                                           'color': theme.colors['fg']}})
                in_fig.viz["figure"].item()

                # Assign the figure format
                in_fig = in_fig.viz["figure"].item()
                update_bokeh_figure(in_fig, plot_cfg)

                if display_check:
                    show(in_fig)

            else:
                SpecSyError(f'Bokeh is not installed')

        case _:
            raise SpecSyError(f'Backend {backend} not supported. Please choose matplotlib or bokeh.')

    return in_fig




def plot_fitted_params(trace_data, true_values=None, output_address=None, backend='matplotlib', n_cols=3,
                       in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None, display_check=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is _NO_FIG else False

    # Load the inference data if necessary
    if isinstance(trace_data, (str, Path)):
        trace_data = az.from_netcdf(trace_data)

    # Parameter labeller
    my_labels = latex_labeller(specsy_cfg['latex_param_notation'], backend)

    # Everything except theo_flux
    summary_params = summary(trace_data, var_names=['~theo_flux'], filter_vars='like')
    var_names = summary_params.index.to_numpy()

    # True values reference dataset
    ref_values = ref_values_dataset(true_values, var_names) if true_values is not None else None

    visuals = {'point_estimate': {'color': theme.colors['fg']},
               'credible_interval': {'color': theme.colors['fg'], 'width': 1},
               'point_estimate_text': {'color': theme.colors['fg']}}

    match backend:
        case 'matplotlib':

            # Set the plot format where the user's overwrites the default
            plot_cfg = theme.fig_defaults(fig_cfg, fig_type=None)

            with rc_context(plot_cfg):
                func = azp.plot_dist(trace_data, var_names=var_names, labeller=my_labels, kind='hist', col_wrap=n_cols,
                                     visuals=visuals, stats={'point_estimate': {'round_to': 4}})
                if ref_values is not None:
                    azp.add_lines(func, values=ref_values, color=['black'])

                in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

        case 'bokeh':
            if bokeh_check:

                plot_cfg = theme.fig_defaults(user_fig=fig_cfg, fig_type=None, plot_lib='bokeh')
                in_fig = azp.plot_dist(trace_data, var_names=var_names, labeller=my_labels, kind='hist',
                                       backend='bokeh', col_wrap=n_cols, visuals=visuals, stats={'point_estimate':
                                                                                                     {'round_to': 'none'}})

                if ref_values is not None:
                    azp.add_lines(in_fig, values=ref_values, color=['black'])

                # Assign the figure format
                in_fig = in_fig.viz["figure"].item()
                update_bokeh_figure(in_fig, plot_cfg)

                if display_check:
                    show(in_fig)

            else:
                SpecSyError(f'Bokeh is not installed')

        case _:
            raise SpecSyError(f'Backend {backend} not supported. Please choose matplotlib or bokeh.')

    return in_fig


def plot_fitted_pairs(trace_data, var_names, true_values=None, output_address=None, backend='matplotlib',
                      in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None, display_check=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is _NO_FIG else False

    # Load the inference data if necessary
    if isinstance(trace_data, (str, Path)):
        trace_data = az.from_netcdf(trace_data)

    # Parameter labeller
    my_labels = latex_labeller(specsy_cfg['latex_param_notation'], backend)

    # True values reference dataset
    ref_values = ref_values_dataset(true_values, var_names) if true_values is not None else None

    # Plot configuration
    visuals = {'scatter': True, 'divergence': True}#, 'contourf': True}

    match backend:
        case 'matplotlib':

            # Set the plot format where the user's overwrites the default
            plot_cfg = theme.fig_defaults(fig_cfg, fig_type=None)

            with rc_context(plot_cfg):
                func = azp.plot_pair(trace_data, var_names=var_names, filter_vars='like', labeller=my_labels,
                                     visuals=visuals)
                if ref_values is not None:
                    azp.add_lines(func, values=ref_values, color=['black'])

                in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

        case 'bokeh':
            if bokeh_check:

                plot_cfg = theme.fig_defaults(user_fig=fig_cfg, fig_type=None, plot_lib='bokeh')

                in_fig = azp.plot_pair(trace_data, var_names=var_names, filter_vars='like', labeller=my_labels,
                                       visuals=visuals, backend='bokeh')

                if ref_values is not None:
                    azp.add_lines(in_fig, values=ref_values, color=['black'])

                # Assign the figure format
                in_fig = in_fig.viz["figure"].item()
                update_bokeh_figure(in_fig, plot_cfg)

                if display_check:
                    show(in_fig)

            else:
                SpecSyError(f'Bokeh is not installed')

        case _:
            raise SpecSyError(f'Backend {backend} not supported. Please choose matplotlib or bokeh.')

    return in_fig



def plot_prior_posterior(trace_data, var_names=None, true_values=None, output_address=None, backend='matplotlib',
                         in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None, display_check=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is _NO_FIG else False

    # Load the inference data if necessary
    if isinstance(trace_data, (str, Path)):
        trace_data = az.from_netcdf(trace_data)

    # Check the trace contains a prior group
    if not hasattr(trace_data, 'prior'):
        raise SpecSyError(f'The input trace does not contain a "prior" group for the prior-posterior comparison')

    # Parameter labeller
    my_labels = latex_labeller(specsy_cfg['latex_param_notation'], backend)

    # Default to all parameters except the line fluxes
    if var_names is None:
        var_names = [var for var in trace_data.posterior.data_vars if not var.startswith('theo_flux')]

    # True values reference dataset
    ref_values = ref_values_dataset(true_values, var_names) if true_values is not None else None

    match backend:
        case 'matplotlib':

            # Set the plot format where the user's overwrites the default
            plot_cfg = theme.fig_defaults(fig_cfg, fig_type=None)

            with rc_context(plot_cfg):
                func = azp.plot_prior_posterior(trace_data, var_names=var_names, labeller=my_labels)
                if ref_values is not None:
                    azp.add_lines(func, values=ref_values, color=['black'])

                in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

        case 'bokeh':
            if bokeh_check:

                plot_cfg = theme.fig_defaults(user_fig=fig_cfg, fig_type=None, plot_lib='bokeh')

                in_fig = azp.plot_prior_posterior(trace_data, var_names=var_names, labeller=my_labels,
                                                  backend='bokeh')
                if ref_values is not None:
                    azp.add_lines(in_fig, values=ref_values, color=['black'])

                # Assign the figure format
                in_fig = in_fig.viz["figure"].item()
                update_bokeh_figure(in_fig, plot_cfg)

                if display_check:
                    show(in_fig)

            else:
                SpecSyError(f'Bokeh is not installed')

        case _:
            raise SpecSyError(f'Backend {backend} not supported. Please choose matplotlib or bokeh.')

    return in_fig


def plot_traces(trace_data, var_names=None, true_values=None, output_address=None, backend='matplotlib',
                in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None, display_check=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is _NO_FIG else False

    # Load the inference data if necessary
    if isinstance(trace_data, (str, Path)):
        trace_data = az.from_netcdf(trace_data)

    # Default to all parameters except the line fluxes
    if var_names is None:
        var_names = [var for var in trace_data.posterior.data_vars if not var.startswith('theo_flux')]

    # Parameter labeller
    my_labels = latex_labeller(specsy_cfg['latex_param_notation'], backend)

    # True values reference dataset
    ref_values = ref_values_dataset(true_values, var_names) if true_values is not None else None

    # One viridis color per chain
    n_chains = trace_data.posterior.sizes['chain']
    viridis_colors = [to_hex(plt.cm.viridis(x)) for x in np.linspace(0, 1, len(var_names))]


    # Plot configuration
    visuals = ({'label': {'rotation': 0}} if backend == 'matplotlib'
               else {'label': {'orientation': 'parallel'}})

    match backend:
        case 'matplotlib':

            # Set the plot format where the user's overwrites the default
            plot_cfg = theme.fig_defaults(fig_cfg, fig_type=None)

            with rc_context(plot_cfg):
                func = azp.plot_trace_dist(trace_data, var_names=var_names, labeller=my_labels, visuals=visuals,
                                           color=viridis_colors)
                if ref_values is not None:
                    azp.add_lines(func, values=ref_values, color=['black'])

                in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

        case 'bokeh':
            if bokeh_check:

                plot_cfg = theme.fig_defaults(user_fig=fig_cfg, fig_type=None, plot_lib='bokeh')

                in_fig = azp.plot_trace_dist(trace_data, var_names=var_names, labeller=my_labels, visuals=visuals,
                                             color=viridis_colors, backend='bokeh')
                if ref_values is not None:
                    azp.add_lines(in_fig, values=ref_values, color=['black'])

                # Assign the figure format
                in_fig = in_fig.viz["figure"].item()
                update_bokeh_figure(in_fig, plot_cfg)

                if display_check:
                    show(in_fig)

            else:
                SpecSyError(f'Bokeh is not installed')

        case _:
            raise SpecSyError(f'Backend {backend} not supported. Please choose matplotlib or bokeh.')

    return in_fig



import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.colors import to_rgba, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import gaussian_filter

# Spectrum, theme and save_close_fig_swicth come from your existing imports

CORNER_COLORS = ['#6a4bb3', '#b22b2b', '#2b8cb2', '#3f9b4a']


def _edges_from_centers(centers):

    # Bin edges halfway between consecutive grid points (for discrete parameters)
    c = np.sort(np.asarray(centers, dtype=float))
    if c.size == 1:
        return np.array([c[0] - 0.5, c[0] + 0.5])
    mid = 0.5 * (c[1:] + c[:-1])
    return np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def _contour_2d(ax, x, y, x_edges, y_edges, color, smooth, sigmas):

    hist, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    if smooth:
        hist = gaussian_filter(hist, smooth)
    if hist.max() <= 0:
        return

    # Density thresholds enclosing the 2D-Gaussian-equivalent mass of each sigma
    mass = 1 - np.exp(-0.5 * np.square(sigmas))
    h_flat = np.sort(hist.ravel())[::-1]
    cdf = np.cumsum(h_flat) / h_flat.sum()
    idx = np.clip(np.searchsorted(cdf, mass, side='left'), 0, h_flat.size - 1)
    levels = np.unique(h_flat[idx])
    levels = levels[levels > 0]
    if levels.size == 0:
        return

    xc = 0.5 * (x_edges[1:] + x_edges[:-1])
    yc = 0.5 * (y_edges[1:] + y_edges[:-1])
    rgb = to_rgba(color)[:3]
    cmap = LinearSegmentedColormap.from_list('corner_fill', [(*rgb, 0.1), (*rgb, 0.85)])

    ax.contourf(xc, yc, hist.T, levels=np.append(levels, hist.max() * 1.01), cmap=cmap)
    ax.contour(xc, yc, hist.T, levels=levels, colors=[color], linewidths=1)

    return


def plot_corner(fig, subplot_spec, samples, labels=None, colors=None, set_labels=None, truths=None,
                grid_points=None, bins=20, smooth=1.0, sigmas=(0.5, 1.0, 1.5, 2.0), space=0.05):
    """
    Lower-triangle scatter plot matrix drawn inside a SubplotSpec of an existing figure.

    samples:     dict {param: 1D array} or list of such dicts (one per model set, same keys)
    labels:      dict {param: axis label}
    set_labels:  list of legend labels, one per sample set
    truths:      dict {param: value} -> black solid lines
    grid_points: dict {param: array} -> gray dotted lines on the diagonal. If no bins are given
                 for that parameter, the bin edges are placed halfway between grid points.
    bins:        int or dict {param: bin edges}
    """

    sets = [samples] if isinstance(samples, dict) else list(samples)
    params = list(sets[0].keys())
    n_par = len(params)
    colors = colors or CORNER_COLORS
    labels, truths, grid_points = labels or {}, truths or {}, grid_points or {}

    # Common bin edges per parameter for all the sample sets
    edges = {}
    for p in params:
        if isinstance(bins, dict) and p in bins:
            edges[p] = np.asarray(bins[p], dtype=float)
        elif p in grid_points:
            edges[p] = _edges_from_centers(grid_points[p])
        else:
            lo = min(np.nanmin(s[p]) for s in sets)
            hi = max(np.nanmax(s[p]) for s in sets)
            pad = 0.05 * (hi - lo) if hi > lo else 0.5
            n_bins = bins if isinstance(bins, int) else 20
            edges[p] = np.linspace(lo - pad, hi + pad, n_bins + 1)

    inner = subplot_spec.subgridspec(n_par, n_par, hspace=space, wspace=space)
    axes = np.full((n_par, n_par), None, dtype=object)

    for i, p_y in enumerate(params):
        for j, p_x in enumerate(params[:i + 1]):
            ax = fig.add_subplot(inner[i, j])
            axes[i, j] = ax
            diag = i == j

            # Distributions
            for k, s in enumerate(sets):
                c = colors[k % len(colors)]
                if diag:
                    ax.hist(s[p_x], bins=edges[p_x], histtype='step', density=True, color=c, linewidth=1.2)
                    ax.axvline(np.nanmedian(s[p_x]), color=c, linestyle='--', linewidth=1)
                else:
                    _contour_2d(ax, s[p_x], s[p_y], edges[p_x], edges[p_y], c, smooth, sigmas)

            # Reference lines
            if diag:
                for g in grid_points.get(p_x, []):
                    ax.axvline(g, color='0.6', linestyle=':', linewidth=0.7, zorder=0)
            if p_x in truths:
                ax.axvline(truths[p_x], color='k', linewidth=1)
            if not diag and p_y in truths:
                ax.axhline(truths[p_y], color='k', linewidth=1)

            # X axis: labels only on the bottom row
            ax.set_xlim(edges[p_x][0], edges[p_x][-1])
            ax.xaxis.set_major_locator(MaxNLocator(4, prune='both'))
            if i < n_par - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(labels.get(p_x, p_x))
                ax.tick_params(axis='x', labelrotation=45)

            # Y axis: no ticks on the diagonal, labels only on the first column
            if diag:
                ax.set_ylim(bottom=0)
                ax.set_yticks([])
            else:
                ax.set_ylim(edges[p_y][0], edges[p_y][-1])
                ax.yaxis.set_major_locator(MaxNLocator(4, prune='both'))
                if j == 0:
                    ax.set_ylabel(labels.get(p_y, p_y))
                    ax.tick_params(axis='y', labelrotation=45)
                else:
                    ax.tick_params(labelleft=False)

    # Legend in the empty upper triangle
    if set_labels is not None and n_par > 1:
        handles = [Line2D([], [], color=colors[k % len(colors)], label=lbl) for k, lbl in enumerate(set_labels)]
        axes[0, 0].legend(handles=handles, loc='upper right', frameon=False,
                          bbox_to_anchor=(n_par + (n_par - 1) * space, 1.0))

    return axes


def plot_ssp_diag(spectrum, samples=None, in_fig=None, fig_cfg=None, ax_cfg=None, corner_cfg=None, label=None,
                  fname=None):

    # Display check for the user figures
    display_check = True if in_fig is None else False

    # Load the spectrum file
    spec = Spectrum.from_file(spectrum, instrument='text')

    # Adjust the default theme
    PLT_CONF = theme.fig_defaults({"figure.dpi": 100, "figure.figsize": [16, 9], "font.size": 12, "axes.labelsize": 13,
                                   "axes.titlesize": 14, "legend.fontsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10})
    AXES_CONF = theme.ax_defaults(ax_cfg, observation=spec)

    # Create and fill the figure
    with (rc_context(PLT_CONF)):

        # Generate the figure object and figures
        if in_fig is None:
            in_fig = plt.figure()

        # Left column: spectrum panels, right column: scatter plot matrix
        outer_grid = in_fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1.3, 1], wspace=0.2)
        grid_ax = outer_grid[0].subgridspec(nrows=2, ncols=1, height_ratios=[2, 1])
        ax = in_fig.add_subplot(grid_ax[0])

        wave_plot = spec.wave.data
        flux_plot = spec.flux.data
        mask_plot = ~spec.flux.mask

        # Plot spectrum
        ax.step(wave_plot, flux_plot, label=label, where='mid', color=theme.colors['fg'],
                linewidth=theme.plt['spectrum_width'])

        ax.update(AXES_CONF)

        # Plot the scatter plot matrix
        if samples is not None:
            plot_corner(in_fig, outer_grid[1], samples, **(corner_cfg or {}))

        # By default, plot on screen unless an output address is provided
        in_fig = save_close_fig_swicth(fname, 'tight', in_fig, False, display_check)

    return