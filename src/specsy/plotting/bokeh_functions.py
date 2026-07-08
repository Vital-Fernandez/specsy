import logging
import arviz as az
import numpy as np
import lime
from pathlib import Path
from specsy.plotting.plots import theme
from specsy.io import SpecSyError, specsy_cfg
from lime.plotting.bokeh_plots import update_bokeh_figure
from itertools import chain
from matplotlib import pyplot as plt, rc_context, cm, colors

_logger = logging.getLogger('SpecSy')


try:
    from bokeh.plotting import figure, output_file, save, show
    from bokeh.layouts import gridplot
    from bokeh.models import GlobalInlineStyleSheet, Span, BoxAnnotation
    from bokeh.io import curdoc
except ImportError:
    bokeh_check = False


# Sentinel object for non input figures
_NO_FIG = object()



def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))


def bokeh_trace(sampling_result, fname=None, var_list=None, exclude=None, var_latex=None, fig_cfg=None, in_fig=_NO_FIG):

    if isinstance(sampling_result, (str, Path)):
        trace = az.from_netcdf(sampling_result)
    else:
        trace = sampling_result

    # Check the lines to plot
    if var_list is None:
        var_list = list(trace.posterior.data_vars)

    # Remove the fluxes
    exclude = ['theo_flux'] if exclude is None else exclude + ['theo_flux']
    if exclude is not None:
        var_list = list(set(var_list) - set(exclude))

    # Sort alphabetically and generate the latex notation
    var_list.sort()
    var_latex = var_latex if var_latex is not None else specsy_cfg['latex_param_notation']
    latex_list = None if var_latex is None else ['$' + var_latex.get(var, var) + '$' for var in var_list]
    mapped_variables = az.labels.MapLabeller(var_name_map=dict(zip(var_list, latex_list)))

    # Display check for input figures and legend
    display_check = True if in_fig is _NO_FIG else False
    legend_check = False

    # Adjust the default theme
    PLT_CONF_bokeh = theme.fig_defaults(user_fig=fig_cfg, fig_type='trace_plot', plot_lib='bokeh')
    tools = PLT_CONF_bokeh.get('tools', "pan,wheel_zoom,box_zoom,reset,save")

    # Generate the figure (arviz uses matplotlib parameteres...)
    with rc_context({'axes.prop_cycle': plt.cycler(color=[theme.colors['fg']])}):

        # Establish figure
        if (in_fig is None) or (in_fig is _NO_FIG):
            fig_arr = [[figure(tools=tools), figure(tools=tools)] for _ in var_list]
        else:
            fig_arr = in_fig

        # Trace plotting function
        az.plot_trace(trace, var_names=var_list, combined=True, labeller=mapped_variables, figsize=(3,2),
                      backend="bokeh", axes=fig_arr, show=False)

        # Assign the figure format
        flat_arr = list(chain.from_iterable(fig_arr))
        update_bokeh_figure(flat_arr, PLT_CONF_bokeh)
        fig_arr = gridplot(fig_arr)

        # Save, display or keep
        if display_check:
            show(fig_arr)

    return fig_arr



def bokeh_scatter_matrix(sampling_result, fname=None, var_list=None, exclude=None, var_latex=None, fig_cfg=None,
                         in_fig=_NO_FIG):

    if isinstance(sampling_result, (str, Path)):
        trace = az.from_netcdf(sampling_result)
    else:
        trace = sampling_result

    # Check the lines to plot
    if var_list is None:
        var_list = list(trace.posterior.data_vars)

    # Remove the fluxes
    exclude = ['theo_flux'] if exclude is None else exclude + ['theo_flux']
    if exclude is not None:
        var_list = list(set(var_list) - set(exclude))

    # Sort alphabetically and generate the latex notation
    var_list.sort()
    var_latex = var_latex if var_latex is not None else specsy_cfg['latex_param_notation']
    latex_list = None if var_latex is None else ['$' + var_latex.get(var, var) + '$' for var in var_list]
    mapped_variables = az.labels.MapLabeller(var_name_map=dict(zip(var_list, latex_list)))

    # Display check for input figures and legend
    display_check = True if in_fig is _NO_FIG else False
    legend_check = False

    # Adjust the default theme
    PLT_CONF_bokeh = theme.fig_defaults(user_fig=fig_cfg, fig_type='scatter_matrix', plot_lib='bokeh')

    # Generate the figure (arviz uses matplotlib parameteres...)
    with rc_context({}):

        fig_arr = az.plot_pair(trace, var_names=var_list, kind='kde', labeller=mapped_variables, backend='bokeh',
                               show=False)

        flat_arr = list(chain.from_iterable(fig_arr))
        update_bokeh_figure(flat_arr, PLT_CONF_bokeh)
        fig_arr = gridplot(fig_arr.tolist())

        # Save, display or keep
        if display_check:
            show(fig_arr)

    return fig_arr


def bokeh_flux_grid(sampling_result, fname=None, n_cols=5, fig_cfg=None, in_fig=_NO_FIG):

    if isinstance(sampling_result, (str, Path)):
        trace = az.from_netcdf(sampling_result)
    else:
        trace = sampling_result

    if 'theo_flux' not in list(trace.posterior.data_vars):
        _logger.warning(f'The input trace does not include the "theo_flux" grid from the sampling')
        return

    # Unpack the flux posterior data
    labels_arr = trace.observed_data['lines'].values
    obs_flux_arr = trace.constant_data['input_flux'].values
    obs_err_arr = trace.constant_data['input_err'].values
    theo_grid = trace.posterior['theo_flux'].values

    # Get the ionic species
    line_list = lime.Line.from_list(labels_arr)
    ion_list = [line.particle.label for line in line_list]
    ion_list = np.unique(ion_list)

    # Declare plot grid size
    n_lines = theo_grid.shape[2]
    n_rows = int(np.ceil(float(n_lines)/float(n_cols)))
    n_cells = n_rows * n_cols

    # Generate color dict
    colorNorm = colors.Normalize(0, ion_list.size)
    cmap = cm.get_cmap(name=theme.colors['mask_map'])
    color_dict = dict(zip(ion_list, np.arange(ion_list.size)))

    # Display check for input figures and legend
    display_check = True if in_fig is _NO_FIG else False
    legend_check = False

    # Set the plot format where the user's overwrites the default
    PLT_CONF = theme.fig_defaults(user_fig=fig_cfg, fig_type='flux_grid', plot_lib='bokeh')

    # Loop through the lines to generate the figures
    fig_grid, row = [], []
    for i in range(n_cells):
        if i < n_lines:

            # Current line
            label = line_list[i].label
            ion = line_list[i].particle.label
            ion_color = rgba_to_hex(cmap(colorNorm(color_dict[ion])))

            # Create figure
            fig = figure(title=label, height=150, width=150, background_fill_color=theme.colors['bg'],
                         border_fill_color=theme.colors['bg'])
            fig.yaxis.visible = False

            # Plot histogram as quad glyphs
            trace = theo_grid[:, :, i].reshape(-1)
            hist, edges = np.histogram(trace, bins=35, density=False)
            fig.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], color=ion_color, fill_alpha=1)

            # Plot observed flux
            if obs_flux_arr is not None:
                inFlux, inErr = obs_flux_arr[i], obs_err_arr[i]
                fig.add_layout(Span(location=inFlux, dimension='height', line_color=theme.colors['fg'], line_width=1.5))
                fig.add_layout(BoxAnnotation(left=inFlux - inErr, right=inFlux + inErr,
                                      fill_color=theme.colors['fg'], fill_alpha=0.2,
                                      line_color=theme.colors['fg'], line_dash='dotted',
                                      line_width=1.5))

        else:
            fig = None

        row.append(fig)
        if len(row) == n_cols:
            fig_grid.append(row)
            row = []

    # Append last incomplete row if any
    if row:
        fig_grid.append(row)

    # Update the figure format
    flat_arr = list(chain.from_iterable(fig_grid))
    update_bokeh_figure(flat_arr, PLT_CONF)
    fig_arr = gridplot(fig_grid)

    # Save, display or keep
    if display_check:
        show(fig_arr)

    return fig_arr