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