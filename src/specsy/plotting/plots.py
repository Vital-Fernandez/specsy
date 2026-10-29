from logging import getLogger

_logger = getLogger('SpecSy')

import numpy as np
import arviz as az

from pathlib import Path
from functools import partial, update_wrapper
from types import SimpleNamespace
from pandas import Series, DataFrame

from matplotlib.font_manager import FontProperties
from matplotlib import pyplot as plt, gridspec, rc_context, cm, colors, rcParams
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea, VPacker
from matplotlib.patches import Rectangle

import lime
from innate import load_dataset
from lime.transitions import label_decomposition
from lime.plotting.plots import save_close_fig_swicth
from lime.plotting.utils import PlotAccessor, _NO_FIG
from lime.plotting.format import Themer
from specsy.variables import Z_SUN


from specsy.plotting.format import theme, plot_vars
from specsy.io import SpecSyError, specsy_cfg, load_trace
from specsy.tools import linear_regression
from specsy.variables import VMS_ZERO_AGE_MYR, Q_LIMITS_LAMBDA, HC_EV_ANGSTROM, MIRI_ION_LINES

# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
import specsy.models.ssp as ssp_lib

try:
  from bokeh.plotting import figure, show, save
  from bokeh.models import ColumnDataSource, Whisker, HoverTool, TeeHead, CDSView, GroupFilter   # (existing)
  from bokeh.models import BoxAnnotation, Span, CustomJS, Div, Range1d, Scatter, Spacer
  from bokeh.layouts import row, column
  from bokeh.resources import CDN
  from lime.plotting.bokeh_plots import update_bokeh_figure
  from lime.plotting.format import theme as lime_theme
  bokeh_check = True
except ImportError:
  bokeh_check = False


try:
    from bokeh.plotting import figure, show, save
    from bokeh.models import ColumnDataSource, Whisker, HoverTool, TeeHead, CDSView, GroupFilter
    bokeh_check = True
except ImportError:
    bokeh_check = False

try:
    import corner
except ImportError:
    corner_check = False




# class PlotAccessor:
#
#     def __init__(self, obj):
#         self._obj = obj
#         self.fig = self.ax = None
#         return
#
#     def show(self):
#         plt.show(self.fig)
#         return



def extinction_gradient(cHbeta, cHbeta_err, frame, rel_Hbeta=True, fname=None, fig_cfg=None, ax_cfg=None, return_fig=False):

    # Extract the parameters from the log
    norm_label = frame.loc[frame.extinction_idcs == 0].index[0]
    norm_latex = frame.loc[frame.extinction_idcs == 0].latex_label.values[0]
    coeff_label = r'Hβ' if rel_Hbeta else norm_latex
    norm_latex = norm_latex.replace('$', '')

    idcs_lines = frame.loc[~np.isnan(frame.extinction_idcs)].index.to_numpy()
    idcs_valid = (frame.loc[idcs_lines].extinction_idcs < 2).to_numpy()

    x_arr = frame.loc[idcs_lines, 'f_lambda'].to_numpy() - frame.loc[norm_label, 'f_lambda']
    y_arr = np.log10(frame.loc[idcs_lines, 'theo_ratio'].to_numpy()) - np.log10(frame.loc[idcs_lines, 'obs_ratio'].to_numpy())
    y_err = (1 / np.log(10)) * (frame.loc[idcs_lines, 'obs_ratio_err'].to_numpy()/frame.loc[idcs_lines, 'obs_ratio'].to_numpy())

    # Linear fitting
    c_HI, c_HI_err, intercept, intercept_err = linear_regression(x_arr[idcs_valid], y_arr[idcs_valid], y_err[idcs_valid])
    linear_fit = c_HI * x_arr + intercept
    linear_label = f'c({coeff_label})={cHbeta:.3f}±{cHbeta_err:.3f}'

    # Adjust the axis labels to include the reference line
    x_label = r'$f_{\lambda} - ' + f'f_{{{norm_latex}}}$'
    y_label = (fr"$\log\left(\frac{{I_{{\lambda}}}}{{I_{{{norm_latex}}}}}\right)_{{\mathrm{{theo}}}}"
               fr"-\log\left(\frac{{F_{{\lambda}}}}{{F_{{{norm_latex}}}}}\right)_{{\mathrm{{obs}}}}$")
    title = fr'c({coeff_label}) extinction calculation'

    if theme.default_lib == 'matplotlib':

        # User configuration overrites user
        PLOT_CONF = theme.fig_defaults()
        AXES_CONF = {'xlabel': x_label, 'ylabel': y_label, 'title': title}

        # Draw the figure
        with rc_context(PLOT_CONF):

            fig, ax = plt.subplots()
            AXES_CONF = AXES_CONF
            ax.set(**AXES_CONF)

            # Plot valid entries
            ax.errorbar(x_arr, y_arr, y_err, fmt='o')

            # Plot excluded entries
            if not np.all(idcs_valid):
                ax.errorbar(x_arr[~idcs_valid], y_arr[~idcs_valid], y_err[~idcs_valid], fmt='o', color='tab:red',
                            label='excluded lines')

            ax.plot(x_arr, linear_fit, linestyle='--', label=linear_label)

            # Labels for the lines
            for i, line_label in enumerate(frame.loc[idcs_lines].index):
                ax.annotate(frame.loc[line_label].latex_label,
                            xy=(x_arr[i], y_arr[i]), xytext=(x_arr[i], y_arr[i] + 1.25 * y_err[i]),
                            horizontalalignment="center", rotation=90,
                            xycoords='data', textcoords=("data", "data"))

            # Legend
            ax.legend(loc=8, ncol=2)

            # Display/save the figure
            save_close_fig_swicth(fname, 'tight', fig_obj=fig)

    elif theme.default_lib == 'bokeh':

        if bokeh_check:

            PLT_CONF = theme.fig_defaults(fig_cfg, plot_lib='bokeh')
            fig = figure(tools=PLT_CONF.get('tools', "pan,wheel_zoom,box_zoom,reset,save"))

            # Valid data set
            source_1 = ColumnDataSource(data=dict(x=x_arr[idcs_valid], y=y_arr[idcs_valid], label=idcs_lines[idcs_valid],
                                        y_upper=y_arr[idcs_valid]+y_err[idcs_valid], y_lower=y_arr[idcs_valid]-y_err[idcs_valid]))
            r1 = fig.scatter('x', 'y', size=8, source=source_1)
            fig.add_layout(Whisker(base="x", upper="y_upper", lower="y_lower", source=source_1, line_color=theme.colors['fg'],
                                   upper_head=TeeHead(size=10, line_color=theme.colors['fg']),
                                   lower_head=TeeHead(size=10, line_color=theme.colors['fg'])))

            if not np.all(idcs_valid):
                source_2 = ColumnDataSource(data=dict(x=x_arr[~idcs_valid], y=y_arr[~idcs_valid], label=idcs_lines[~idcs_valid],
                                            y_upper=y_arr[~idcs_valid]+y_err[~idcs_valid], y_lower=y_arr[~idcs_valid]-y_err[~idcs_valid]))
                r2 = fig.scatter('x', 'y', size=8, source=source_2,  legend_label='Excluded lines', color='red')
                fig.add_layout(Whisker(base="x", upper="y_upper", lower="y_lower", source=source_2, line_color=theme.colors['fg'],
                                       upper_head=TeeHead(size=10, line_color=theme.colors['fg']),
                                       lower_head=TeeHead(size=10, line_color=theme.colors['fg'])))
                renderers = [r1, r2]
            else:
                renderers = [r1]

            # Hover tool for all points
            fig.add_tools(HoverTool(renderers=renderers, tooltips=[("Label", "@label")]))

            # Linear fit
            fig.line(x_arr, linear_fit, line_dash='dashed', line_width=2, legend_label=linear_label.replace('$', '$$'),
                     color=theme.colors['regression'])

            fig.legend.location = "bottom_center"
            fig.legend.orientation = "horizontal"

            # Plot labels
            fig.xaxis.axis_label = x_label.replace("$", "$$")
            fig.yaxis.axis_label = y_label.replace("$", "$$")

            # # Adjust the format of the plot TODO fix this one
            # update_bokeh_figure(fig, PLT_CONF)

            # Save or display the plot
            if return_fig:
                return fig

            elif fname is not None:
                save(fname, filename=fname)

            else:
                # output_notebook()
                show(fig)

        else:
            raise SpecSyError(f'Please install Bokeh or use "plot" (matplotlib) as the default library"')

    else:
        raise SpecSyError(f'Default plotting library: {theme.default_lib} is not recognized please use "plot" or "bokeh"')

    return


def parameter_notation(param, mean, std):

    # Label for the plot
    if mean > 10:
        label = r'{} = ${:.0f}$$\pm${:.0f}'.format(_setup_cfg['latex'][param], mean, std)
    else:
        label = r'{} = ${:.3f}$$\pm${:.3f}'.format(_setup_cfg['latex'][param], mean, std)

    return label


def numberStringFormat(value, cifras = 4):
    if value > 0.001:
        newFormat = f'{value:.{cifras}f}'
    else:
        newFormat = f'{value:.{cifras}e}'

    return newFormat


# ----------------- Chemical sampler plots -----------------

def plot_traces(fname, output_address=None, params_list=None, true_values=None, n_cols=1, n_rows=None, col_row_scale=(10, 4),
                in_fig=None, fig_cfg=None, ax_cfg=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is None else False

    # Load the inference data
    # infer_db = load_inference_data(fname)
    infer_db = load_dataset(fname)

    # Check for true values
    if true_values is None:
        if 'true_values' in infer_db:
             true_values = dict(zip(infer_db.true_values.parameters.values, infer_db.true_values.magnitude.values))

    # Set the number of parameters to plot the
    chain_params = list(infer_db.posterior.data_vars)
    if params_list is None:
        input_params = [param for param in chain_params if '_Op' not in param]
    else:
        input_params = []
        for param in params_list:
            if param in chain_params:
                input_params.append(param)
    n_traces = len(input_params)

    # Compute the number of rows configuration
    if n_traces > n_cols:
        if n_rows is None:
            n_rows = int(np.ceil(n_traces / n_cols))
    else:
        n_cols, n_rows = n_traces, 1
    n_grid = n_cols * n_rows

    # Set the plot format where the user's overwrites the default
    size_conf = {'figure.figsize': (8, n_traces)}
    size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}

    plot_cfg = theme.fig_defaults(size_conf, fig_type='traces')
    # ax_cfg = theme.ax_defaults(fig_type='traces')

    # Initialize the figure
    with (rc_context(plot_cfg)):

        # Plot format
        # Generate the figure if not provided
        if in_fig is None:
            in_fig = plt.figure()
        gs = gridspec.GridSpec(n_traces * 2, 4)
        gs.update(wspace=0.2, hspace=1.8)

        # Colors
        colorNorm = colors.Normalize(0, n_traces)
        cmap = cm.get_cmap(name=theme.colors['mask_map'])

        for i in range(n_grid):

            if i < n_traces:

                param = input_params[i]
                trace_array = infer_db.posterior[param].values
                trace_array = trace_array.reshape(-1)

                mean_value = np.mean(trace_array)
                std_dev = np.std(trace_array)

                axTrace = in_fig.add_subplot(gs[2 * i:2 * (1 + i), :3])
                axPoterior = in_fig.add_subplot(gs[2 * i:2 * (1 + i), 3])

                param_latex = _setup_cfg['latex'][param]
                label_measurement = parameter_notation(param, mean_value, std_dev)

                # Plot the traces
                axTrace.plot(trace_array, label=label_measurement, color=cmap(colorNorm(i)))
                axTrace.axhline(y=mean_value, color=cmap(colorNorm(i)), linestyle='--')
                axTrace.set_ylabel(param_latex)

                # Plot the histograms
                axPoterior.hist(trace_array, bins=50, histtype='step', color=cmap(colorNorm(i)), align='left')

                # Plot the axis as percentile
                median, percentile16th, percentile84th = np.percentile(trace_array, (50, 16, 84))

                # Add true value if available
                if true_values is not None:
                    if param in true_values:
                        value_param = true_values[param]

                        # Nominal value and uncertainty
                        if isinstance(value_param, (list, tuple, np.ndarray)):
                            nominal_value, std_value = value_param[0], 0.0 if len(value_param) == 1 else value_param[1]
                            axPoterior.axvline(x=nominal_value, color=theme.colors['fg'], linestyle='solid')
                            axPoterior.axvspan(nominal_value - std_value, nominal_value + std_value, alpha=0.5,
                                               color=cmap(colorNorm(i)))

                        # Nominal value only
                        else:
                            nominal_value = value_param
                            axPoterior.axvline(x=nominal_value, color=theme.colors['fg'], linestyle='solid')

                # Add legend
                axTrace.legend(loc=7)

                # Remove ticks and labels
                if i < n_traces - 1:
                    axTrace.get_xaxis().set_visible(False)
                    axTrace.set_xticks([])

                axPoterior.yaxis.set_major_formatter(plt.NullFormatter())
                axPoterior.set_yticks([])

                axPoterior.set_xticks([percentile16th, median, percentile84th])
                round_n = 0 if median > 10 else 3
                axPoterior.set_xticklabels(['', numberStringFormat(median, round_n), ''])

                axTrace.set_yticks((percentile16th, median, percentile84th))
                round_n = 0 if median > 10 else 3
                axTrace.set_yticklabels((numberStringFormat(percentile16th, round_n), '',
                                         numberStringFormat(percentile84th, round_n)))

        # Show or save the image
        in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

    return in_fig


def plot_flux_grid(fname, output_address=None, line_list=None, obs_values=None, n_cols=8,  n_rows=None, combined_dict=None,
                   in_fig=None, fig_cfg=None, ax_cfg=None, maximize=False):

    # Display check for the user figures
    display_check = True if in_fig is None else False

    # Load the inference data
    if isinstance(fname, (str, Path)):
        infer_db = az.from_netcdf(fname)
    else:
        infer_db = fname

    # Recover the fluxes
    input_traces = infer_db.posterior.calcFluxes_Op.values

    # Check for input fluxes and errors
    input_lines = infer_db.inputs.labels.values
    input_fluxes = infer_db.inputs.fluxes.values
    input_errs = infer_db.inputs.errs.values

    # Crop the line list if requested
    if line_list is not None:
        mask = np.isin(input_lines, line_list)
        input_lines, input_fluxes, input_errs = input_lines[mask], input_fluxes[mask], input_errs[mask]
        input_traces = input_traces[:, :, mask]

    # Recover the true values
    if obs_values is None:
        if 'inputs' in infer_db:
            labels_true = infer_db.inputs.labels.values
            fluxes_true, err_true = infer_db.inputs.fluxes.values, infer_db.inputs.errs.values
            obs_values = {key: (val1, val2) for key, val1, val2 in zip(labels_true, fluxes_true, err_true)}

    # Get ions to group the colors
    ion_array, latexLabel_array = label_decomposition(input_lines, params_list=('particle', 'latex_label'))

    # Declare plot grid size
    n_lines = len(input_lines)
    n_rows = int(np.ceil(float(n_lines)/float(n_cols)))
    n_cells = n_rows * n_cols

    # Set the plot format where the user's overwrites the default
    size_conf = {'figure.figsize': (n_cols, n_rows)}
    size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}

    plot_cfg = theme.fig_defaults(size_conf, fig_type='flux_grid')
    # ax_cfg = theme.ax_defaults(fig_type='traces')

    # Initialize the figure
    with (rc_context(plot_cfg)):

        # Generate the color dict
        input_ions = np.unique(ion_array)
        colorNorm = colors.Normalize(0, input_ions.size)
        cmap = cm.get_cmap(name=theme.colors['mask_map'])
        color_dict = dict(zip(input_ions, np.arange(input_ions.size)))

        # self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=n_columns, n_rows=n_rows)
        if in_fig is None:
            in_fig = plt.figure()

        axes = in_fig.subplots(n_rows, n_cols)
        # axes = plt.subplots(n_rows, n_cols)
        axes = axes.ravel()

        # Plot individual traces
        for i in range(n_cells):

            if i < n_lines:

                # Current line
                label, ion = input_lines[i], ion_array[i]
                ion_color = cmap(colorNorm(color_dict[ion]))

                # Plot histogram
                trace = input_traces[:, :, i].reshape(-1)
                median_flux = np.median(trace)

                axes[i].hist(trace, histtype='stepfilled', bins=35, alpha=.7, color=ion_color, density=False)

                # Plot observed flux
                if obs_values is not None:
                    if label in obs_values:

                        inFlux, inErr = obs_values[label]
                        label_mean = 'Mean value: {}'.format(np.around(median_flux, 4))
                        label_true = 'True value: {}'.format(np.around(inFlux, 3))
                        axes[i].axvline(x=inFlux, label=label_true, color=theme.colors['fg'], linestyle='solid')
                        axes[i].axvspan(inFlux - inErr, inFlux + inErr, alpha=0.5, color='grey')

                # Plot formating
                axes[i].get_yaxis().set_visible(False)
                axes[i].set_yticks([])
                axes[i].set_title(latexLabel_array[i])

            else:
                in_fig.delaxes(axes[i])

        # Show or save the image
        in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)

    return in_fig


def plot_corner_matrix(fname, output_address=None, params_list=None, true_values=None, in_fig=None, fig_cfg=None,
                       ax_cfg=None, maximize=False):

    if corner_check:

        # Display check for the user figures
        display_check = True if in_fig is None else False

        # Load the inference data
        if isinstance(fname, (str, Path)):
            infer_db = az.from_netcdf(fname)
        else:
            infer_db = fname

        # Set the parameters to plot
        chain_params = np.array(list(infer_db.posterior.data_vars))
        if params_list is not None:
            idcs_plot = np.isin(chain_params, params_list)
        else:
            idcs_plot = np.char.find(chain_params, '_Op') == -1
        input_params = chain_params[idcs_plot]

        # Check for true values
        if true_values is None:
            if 'true_values' in infer_db:
                true_values = dict(zip(infer_db.true_values.parameters.values, infer_db.true_values.magnitude.values))

        # Prepare corner arrays
        labels_list, traces_list = [], []
        true_array = None if true_values is None else []

        for i, param in enumerate(input_params):
            labels_list.append(_setup_cfg['latex'].get(param, param))
            traces_list.append(infer_db.posterior[param].values.reshape(-1))

            if true_array is not None:
                true_array.append(true_values[param])

        # Change to numpy and transpose
        labels_list, traces_list = np.array(labels_list), np.array(traces_list).T
        true_array = None if true_array is None else np.array(true_array)

        # Set the plot format where the user's overwrites the default
        plot_cfg = theme.fig_defaults(fig_cfg)
        # ax_cfg = theme.ax_defaults()

        # Initialize the figure
        with (rc_context(plot_cfg)):

            # Generate the plot
            corner.corner(traces_list, fontsize=30, labels=labels_list, quantiles=[0.16, 0.5, 0.84],
                                show_titles=True, title_args={"fontsize": 200}, truths=true_array,
                                truth_color=theme.colors['fg'], title_fmt='0.3f', fig=in_fig)

            # Show or save the image
            in_fig = save_close_fig_swicth(output_address, 'tight', in_fig, maximize, display_check)


        # Dark models
        # # Declare figure format
        # background = np.array((43, 43, 43)) / 255.0
        # foreground = np.array((179, 199, 216)) / 255.0
        #
        # figConf = {'text.color': foreground,
        #            'figure.figsize': (16, 10),
        #            'figure.facecolor': background,
        #            'axes.facecolor': background,
        #            'axes.edgecolor': foreground,
        #            'axes.labelcolor': foreground,
        #            'axes.labelsize': 30,
        #            'xtick.labelsize': 12,
        #            'ytick.labelsize': 12,
        #            'xtick.color': foreground,
        #            'ytick.color': foreground,
        #            'legend.edgecolor': 'inherit',
        #            'legend.facecolor': 'inherit',
        #            'legend.fontsize': 16,
        #            'legend.loc': "center right"}
        # rcParams.update(figConf)
        # # Generate the plot
        # mykwargs = {'no_fill_contours':True, 'fill_contours':True}
        # self.Fig = corner.corner(traces_array[:, :], fontsize=30, labels=labels_list, quantiles=[0.16, 0.5, 0.84],
        #                          show_titles=True, title_args={"fontsize": 200},
        #                          truth_color='#ae3135', title_fmt='0.3f', color=foreground, **mykwargs)#, hist2d_kwargs = {'cmap':'RdGy',
        #                                                                                    #'fill_contours':False,
        #                                                                                    #'plot_contours':False,
        #                                                                                    #'plot_datapoints':False})



        # plt.savefig(plot_address, dpi=100, bbox_inches='tight')
        # plt.close(fig)

    else:
        SpecSyError(f'Please install corner to generate the scatter plot matrix')

    return in_fig


def az_trace(sampling_result, fname=None, in_fig=None, var_list=None, exclude=['theo_flux'],
             var_latex=plot_vars['latex'], fig_cfg=None, maximize=False, display_check=True):

    if isinstance(sampling_result, (str, Path)):
        trace = az.from_netcdf(sampling_result)
    else:
        trace = sampling_result

    # Check the lines to plot
    if var_list is None:
        var_list = list(trace.posterior.data_vars)

    # Remove the fluxes
    if exclude is not None:
        var_list = list(set(var_list) - set(exclude))

    # Sort alphabetically and generate the latex notation
    var_list.sort()
    latex_list = None if var_latex is None else [var_latex.get(var, var) for var in var_list]

    mapped_variables = az.labels.MapLabeller(var_name_map=dict(zip(var_list, latex_list)))

    # Set the plot format where the user's overwrites the default
    size_conf = {'figure.figsize': (6, len(var_list) * 2),
                 'axes.prop_cycle': plt.cycler(color=[theme.colors['fg']]),
                 'axes.titlesize': 16, 'axes.titlepad': 10}
    size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}

    plot_cfg = theme.fig_defaults(size_conf, fig_type='traces')
    print(plot_cfg)

    with rc_context(plot_cfg):

        # Arviz function
        az.plot_trace(trace, var_names=var_list, combined=True, labeller=mapped_variables)
        plt.subplots_adjust(hspace=0.8)

        # Display | saving logic
        in_fig = save_close_fig_swicth(fname, 'tight', in_fig, maximize, display_check)




    return


def az_scatter_matrix(sampling_result, fname=None, in_fig=None, var_list=None, exclude=['theo_flux'],
                      var_latex=plot_vars['latex'], fig_cfg=None, maximize=False, display_check=True, n_cols=5):

    if isinstance(sampling_result, (str, Path)):
        trace = az.from_netcdf(sampling_result)
    else:
        trace = sampling_result

    # Check the lines to plot
    if var_list is None:
        var_list = list(trace.posterior.data_vars)

    # Remove the fluxes
    if exclude is not None:
        var_list = list(set(var_list) - set(exclude))

    # Sort alphabetically and generate the latex notation
    var_list.sort(reverse=True)
    latex_list = None if var_latex is None else [var_latex.get(var, var) for var in var_list]

    mapped_variables = az.labels.MapLabeller(var_name_map=dict(zip(var_list, latex_list)))


    size_conf = {}
    size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}

    plot_cfg = theme.fig_defaults(size_conf, fig_type='traces')


    with rc_context(plot_cfg):

        # Arviz function
        az.plot_pair(trace, var_names=var_list, kind='kde', labeller=mapped_variables)

        # Display | saving logic
        in_fig = save_close_fig_swicth(fname, 'tight', in_fig, maximize, display_check)

    return


def az_flux_grid(sampling_result, fname=None, n_cols=5, fig_cfg=None, in_fig=None,  maximize=False,
                 display_check=True):

    if isinstance(sampling_result, (str, Path)):
        trace = az.from_netcdf(sampling_result)
    else:
        trace = sampling_result

    if 'theo_flux' not in list(trace.posterior.data_vars):
        raise SpecSyError(f'The input trace does not include the "theo_flux" grid')

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

    # Set the plot format where the user's overwrites the default
    size_conf = {'figure.figsize': (n_cols*2, n_rows*2)}
    size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}

    plot_cfg = theme.fig_defaults(size_conf, fig_type='flux_grid')

    # Initialize the figure
    with (rc_context(plot_cfg)):

        # Generate the color dict
        colorNorm = colors.Normalize(0, ion_list.size)
        cmap = cm.get_cmap(name=theme.colors['mask_map'])
        color_dict = dict(zip(ion_list, np.arange(ion_list.size)))

        # self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=n_columns, n_rows=n_rows)
        if in_fig is None:
            in_fig = plt.figure()

        axes = in_fig.subplots(n_rows, n_cols)
        axes = axes.ravel()

        # Plot individual traces
        for i in range(n_cells):

            if i < n_lines:

                # Current line
                label, ion = line_list[i].label, line_list[i].particle.label
                ion_color = cmap(colorNorm(color_dict[ion]))

                # Plot histogram
                trace = theo_grid[:, :, i].reshape(-1)
                axes[i].hist(trace, histtype='stepfilled', bins=35, alpha=.7, color=ion_color, density=False)

                # Plot observed flux
                if obs_flux_arr is not None:
                    inFlux, inErr = obs_flux_arr[i], obs_err_arr[i]
                    axes[i].axvline(x=inFlux, color=theme.colors['fg'], linestyle='solid')
                    axes[i].axvspan(inFlux - inErr, inFlux + inErr, alpha=0.2, edgecolor=theme.colors['fg'],
                                    linewidth=1.5, linestyle=':', color=theme.colors['fg'])

                # Plot formating
                axes[i].get_yaxis().set_visible(False)
                axes[i].set_yticks([])
                axes[i].set_title(label)

            else:
                in_fig.delaxes(axes[i])

        # Bigger head spaces
        plt.subplots_adjust(hspace=0.5)

        # Show or save the image
        in_fig = save_close_fig_swicth(fname, 'tight', in_fig, maximize, display_check)

    return


# ----------------- Stellar binaries plots -----------------
def plot_ionizing_hardness(self, ages_myr, libraries=None, energies=None, ion_lines=MIRI_ION_LINES, z_sun=None,
                           ylim=None, lib_labels=None, fig_cfg=None, fname=None):

    # Loaded spectra of the requested libraries and ages
    selection = self.frame.loc[self._library_mask(libraries, 'plot') & self._loaded_mask()]
    grid_ages = np.unique(selection['age'])

    # Nearest grid age to each request
    ages = np.log10(np.where(np.atleast_1d(ages_myr) == 0, VMS_ZERO_AGE_MYR, np.atleast_1d(ages_myr)) * 1e6)
    nearest = grid_ages[np.abs(grid_ages[:, None] - ages[None, :]).argmin(axis=0)]
    missing = ages[np.abs(nearest - ages) > 0.02]
    if missing.size > 0:
        raise SpecSyError(f'No grid age within 0.02 dex of log(age)={missing.tolist()}. Available ages (Myr): '
                          f'{[float(f"{10 ** a / 1e6:.3g}") for a in grid_ages]}')

    # Unpack the values
    ages = np.unique(nearest)
    mets = np.unique(selection['metallicity'])
    libs = sorted(selection['library'].unique(), key=str)

    # Check if the library has the nebular continuum
    neb_libs = [lib for lib in selection['library'].unique() if str(lib).endswith('_neb')]
    if neb_libs:
        _logger.warning(f'The selection includes the nebular libraries {neb_libs}. Their flux below '
                        f'{Q_LIMITS_LAMBDA["HI"]:.1f} Å may not be the stellar ionizing output.')

    # Compute the ionization steps
    energies = (13.6, 60.01) if energies is None else energies
    energies = np.arange(energies[0], energies[1], 0.2)
    edges = HC_EV_ANGSTROM / energies

    # Figure format variables
    n_cols = min(len(ages), 2)
    n_rows = int(np.ceil(len(ages) / n_cols))

    colors = {met: f'C{i}' for i, met in enumerate(mets)}
    styles = {lib: ('-', '--', ':', '-.')[i % 4] for i, lib in enumerate(libs)}
    lib_labels = {} if lib_labels is None else lib_labels
    z_label = (lambda met: rf'$Z = {met / z_sun:.2g}\,Z_\odot$') if z_sun is not None else (lambda met: f'Z = {met}')

    # Generate the figure
    with rc_context(fig_cfg):

        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, squeeze=False)

        # Check if the spectrum can produce ionization photons
        without_Q_counter = 0
        for ax, age in zip(axes.ravel(), ages):
            group = selection.loc[np.isclose(selection['age'], age)]
            for idx, row in group.iterrows():

                wave_arr, flux_arr = self._spectrum(idx)
                log_q0 = log_ionizing_photons(wave_arr, flux_arr, units_wave=self.units_wave, units_flux=self.units_flux)
                if not np.isfinite(log_q0):
                    without_Q_counter += 1
                    continue

                # Energies without photons are NaN, so the curve stops there
                log_q = np.array([log_ionizing_photons(wave_arr, flux_arr, self.units_wave, self.units_flux,
                                                       lambda_limit=edge) for edge in edges])
                ratio = np.where(np.isfinite(log_q), 10 ** (log_q - log_q0), np.nan)
                ax.plot(energies, ratio, color=colors[row['metallicity']], linestyle=styles[row['library']])

            # Ionization energies of the diagnostic lines
            for label, energy in (ion_lines or {}).items():
                if energies[0] <= energy <= energies[-1]:
                    ax.axvline(energy, color='black', linewidth=0.7, linestyle=':')
                    ax.text(energy, 0.03, label, rotation=90, ha='right', va='bottom', fontsize=8,
                            transform=ax.get_xaxis_transform())

            ax.set_yscale('log')
            ax.set_xlim(energies[0], energies[-1])
            if ylim is not None:
                ax.set_ylim(*ylim)
            ax.set_title(f'{10 ** age / 1e6:.3g} Myr')

        if without_Q_counter > 0:
            _logger.warning(f'{without_Q_counter} spectra have no flux below {Q_LIMITS_LAMBDA["HI"]:.1f} Å and were not plotted.')

        # Unused panels
        for ax in axes.ravel()[len(ages):]: ax.set_visible(False)

        # Label the axes
        for ax in axes[-1, :]: ax.set_xlabel('Energy, E [eV]')
        for ax in axes[:, 0]: ax.set_ylabel(r'$Q(>E)\,/\,Q(>13.6\,\mathrm{eV})$')

        # Colors for the metallicity and line styles for the library (if several)
        handles = [Line2D([], [], color=colors[met], label=z_label(met)) for met in mets]
        if len(libs) > 1:
            handles += [Line2D([], [], color='0.3', linestyle=styles[lib], label=lib_labels.get(lib, lib))
                        for lib in libs]

        # Legend in first plot
        axes[0, 0].legend(handles=handles, loc='upper right', fontsize='small')

        plt.tight_layout()
        if fname is None:
            plt.show()
        else:
            plt.savefig(fname, bbox_inches='tight')
            plt.close(fig)

    return


def plot_spectra(library_mask, loaded_mask, frame_df, source, spectrum, fname=None, metallicity=None, age=None, metallicity_range=None,
                 age_range=None, libraries=None, log_scale=False, title=None, fig_cfg=None, in_fig=_NO_FIG, in_ax=None, maximize=False):

    # Display check for input figures
    display_check = True if in_fig is _NO_FIG else False

    # Select loaded rows matching the filters (exact values and inclusive ranges)
    mask = library_mask(libraries, 'plot') & loaded_mask()
    if metallicity is not None:
        mask &= frame_df['metallicity'].isin(np.atleast_1d(metallicity)).to_numpy()
    if metallicity_range is not None:
        mask &= frame_df['metallicity'].between(*metallicity_range).to_numpy()
    if age is not None:
        mask &= frame_df['age'].isin(np.atleast_1d(age)).to_numpy()
    if age_range is not None:
        mask &= frame_df['age'].between(*age_range).to_numpy()

    if not mask.any():
        raise SpecSyError('No loaded binaries match the requested filters.')

    selection = frame_df.loc[mask].sort_values(['metallicity', 'age'])

    # Figure format
    plt_cfg = theme.fig_defaults(fig_cfg, fig_type='metallicity_age_map') # TODO add units
    ax_cfg = {'xlabel': 'Wavelength (Å)', 'ylabel': 'Flux',
              'title': f'Stellar binary spectra — source: {source}' if title is None else title}

    with rc_context(plt_cfg):

        # Infer figure structure
        if in_ax is not None:
            fig, ax = in_ax.figure, in_ax
        else:
            fig = plt.figure() if (in_fig is None) or (in_fig is _NO_FIG) else in_fig
            ax = fig.add_subplot()
        ax.set(**ax_cfg)

        # Color cycle across unique metallicities for visual grouping
        unique_mets = np.unique(selection['metallicity'])
        met_color = {m: f'C{i}' for i, m in enumerate(unique_mets)}

        legend_entries = {}
        for idx, row in selection.iterrows():
            wave_arr, flux_arr = spectrum(idx)
            label_key = f'Z={row["metallicity"]:.5f} | {row["library"]}'
            ax.step(wave_arr, flux_arr, where='mid', alpha=0.7,
                         linewidth=0.8, color=met_color[row['metallicity']],
                         label=label_key if label_key not in legend_entries else None)
            legend_entries[label_key] = True

        if log_scale:
            ax.set_yscale('log')

        if len(legend_entries) <= 20:
            ax.legend(title='Metallicity | Library')

        save_close_fig_swicth(fname, 'tight', fig, maximize, display_check)

    return fig, ax



def plot_metallicity_age_map(frame_df, library_mask, libraries=None, in_fig=None, in_ax=None, fig_cfg=None, title=None,
                             fname=None, maximize=False):

    # Prepare the data
    selection = frame_df.loc[library_mask(libraries, 'plot')]
    library_groups = selection.groupby('library', sort=False, dropna=False)
    use_colors = library_groups.ngroups > 1
    sources = sorted(frame_df['source'].unique(), key=str)
    title = f"Age-Metallicity coverage — source: {', '.join(map(str, sources))}" if title is None else title

    # Figure format
    plt_cfg = theme.fig_defaults(fig_cfg, fig_type='metallicity_age_map')
    ax_cfg = {'title': title, 'xlabel': "log(Age)", 'ylabel': "Metallicity (Z)"}

    with rc_context(plt_cfg):

        # Infer figure structure
        if in_ax is not None:
            fig, ax = in_ax.figure, in_ax
        else:
            fig = plt.figure() if (in_fig is None) or (in_fig is _NO_FIG) else in_fig
            ax = fig.add_subplot()
        ax.set(**ax_cfg)

        for i, (lib, group) in enumerate(library_groups):
            ax.scatter(group['age'], group['metallicity'], alpha=0.6, s=15, color=f'C{i}' if use_colors else 'C0',
                       label=lib if use_colors else None)

        # Figure formating
        ax.update(ax_cfg)

        if use_colors:
            ax.legend(title='Library', markerscale=2)

        # By default, plot on screen unless an output address is provided
        if in_fig is not None:
            save_close_fig_swicth(fname, 'tight', fig, maximize)

    return fig, ax


class StellarBinariesPlots(PlotAccessor):

    """Plots of a StellarBinaries object: spectra, ionizing_hardness, metallicity_age_map."""

    def metallicity_age_map(self, fname=None, libraries=None, title=None, in_fig=_NO_FIG, in_ax=None, fig_cfg=None,
                            maximize=False):

        """
        Age-metallicity coverage of the stellar binaries, one point per spectrum.

        Parameters
        ----------
        fname : str, optional
            Output file. If None the figure is shown on screen.
        libraries : str or list, optional
            Libraries to plot (all if None). Different libraries are drawn with different colors.
        title : str, optional
            Figure title. The default lists the source(s).
        in_fig, in_ax : optional
            Existing figure or axes to draw on.
        fig_cfg : dict, optional
            Matplotlib rcParams overrides.
        maximize : bool, optional
            Maximize the figure window.

        Returns
        -------
        fig, ax
        """
        with rc_context():

            self.fig, self.ax = plot_metallicity_age_map(self._obj.frame, self._obj._library_mask, libraries,
                                                         in_fig=in_fig, in_ax=in_ax, fig_cfg=fig_cfg, title=title,
                                                         fname=fname, maximize=maximize)

            self.save_close_fig_swicth(fname, 'tight', self.fig, maximize, in_fig=in_fig, in_ax=in_ax)

        return

    def spectra(self, fname=None, metallicity=None, age=None, metallicity_range=None, age_range=None, libraries=None,
                log_scale=False, title=None, fig_cfg=None, in_fig=_NO_FIG, in_ax=None, maximize=False):


        with rc_context():

            self.fig, self.ax = plot_spectra(self._obj._library_mask, self._obj._loaded_mask, self._obj.frame,
                                             self._obj.source, self._obj._spectrum, fname=fname,
                                             metallicity=metallicity, age=age, metallicity_range=metallicity_range, age_range=age_range,
                                             libraries=libraries, log_scale=log_scale, title=title, fig_cfg=fig_cfg,
                                             in_fig=in_fig, in_ax=in_ax, maximize=maximize)

            self.save_close_fig_swicth(fname, 'tight', self.fig, maximize, in_fig=in_fig, in_ax=in_ax)


        return


# ----------------- SSP sampler -----------------


def _labels(z_sun):
    labels = dict(plot_vars['latex'])
    if z_sun is not None:
        labels['log_Z'] = r'Z/Z_\odot'
    return labels


def _fontsize(key, scale=1.0):
    """rcParams font size in points times `scale` (works with named sizes such as 'large' or 'medium')."""
    return FontProperties(size=rcParams[key]).get_size_in_points() * scale


def _relative_residuals(y, model):
    """(y - model) / y, with NaN where it is not finite (e.g. zero flux in chip gaps)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        res = (y - model) / y
    return np.where(np.isfinite(res), res, np.nan)


def _mask_spans(x, mask):
    """(start, end) wavelengths of each run of masked pixels at the pixel edges (mask True = fitted)."""
    edges = np.concatenate([[1.5 * x[0] - 0.5 * x[1]], 0.5 * (x[1:] + x[:-1]), [1.5 * x[-1] - 0.5 * x[-2]]])
    runs = np.flatnonzero(np.diff(np.r_[False, ~mask, False])).reshape(-1, 2)
    return [(edges[lo], edges[hi]) for lo, hi in runs]


def _value_label(name, p16, p50, p84, solar=False):
    """Two-line LaTeX label: the parameter name on the first line, 'p50^{+up}_{-low}' on the second.
    Metallicity is stored as log10(Z) and shown in linear scale; with `solar` the values are already Z/Z_sun."""

    if name == 'log_Z':
        if not solar:
            p16, p50, p84 = 10.0 ** p16, 10.0 ** p50, 10.0 ** p84
        label, max_d = (r'Z/Z_\odot' if solar else 'Z'), 5
    else:
        label, max_d = plot_vars['latex'].get(name, name), 4

    up, lo = p84 - p50, p50 - p16
    errs = [e for e in (up, lo) if e > 0]
    d = int(np.clip(1 - np.floor(np.log10(min(errs))), 1, max_d)) if errs else 2
    return f"${label}$\n${p50:.{d}f}^{{+{up:.{d}f}}}_{{-{lo:.{d}f}}}$"


def _conditions_title(fig, axes, fitted, mouse=None, mouse_color='magenta', z_sun=None):

    """
    Title centered over `axes` (the line zooms or the residuals panel) with the physical conditions: the fitted value
    of each parameter in the color of the fitted spectrum and, if there is a mouse selection, its value in parentheses
    in the color of the selection. It replaces the previous title. `fitted` and `mouse` are (log_age, log_Z, E(B-V),
    log_A); Z is shown in linear scale.
    """

    for old in [artist for artist in fig.get_children() if artist.get_gid() == 'conditions_title']:
        old.remove()

    plain, fit_color = rcParams['text.color'], theme.colors['fitted_SSP']
    labels, formats = plot_vars['bokeh']['CONDITION_LABELS'], plot_vars['bokeh']['CONDITION_FORMATS']
    if z_sun is not None:
        labels[1], formats[1] = r'Z/Z$_\odot$', '.3f'
    z_norm = 1.0 if z_sun is None else z_sun
    to_linear = lambda values: (values[0], 10.0 ** values[1] / z_norm, values[2], values[3])

    # Colored pieces of text, parameter by parameter
    fitted_values = to_linear(fitted)
    mouse_values = None if mouse is None else to_linear(mouse)
    pieces = []
    for k, (label, fmt) in enumerate(zip(labels, formats)):
        pieces += [(f'{label} = ', plain), (format(fitted_values[k], fmt), fit_color)]
        if mouse_values is not None:
            pieces += [(' (', plain), (format(mouse_values[k], fmt), mouse_color), (')', plain)]
        if k < len(plot_vars['bokeh']['CONDITION_LABELS']) - 1:
            pieces.append((',  ', plain))

    # Matplotlib text has a single color: the pieces are packed side by side on a common baseline
    size = _fontsize('axes.titlesize', 0.8)
    packer = HPacker(children=[TextArea(text, textprops=dict(color=color, fontsize=size)) for text, color in pieces],
                     align='center', pad=0, sep=0)

    x0, x1 = axes[0].get_position().x0, axes[-1].get_position().x1
    y1 = max(ax.get_position().y1 for ax in axes)
    title = AnchoredOffsetbox(loc='lower center', child=packer, pad=0, borderpad=0, frameon=False,
                              bbox_to_anchor=((x0 + x1) / 2, y1 + 0.01), bbox_transform=fig.transFigure)
    title.set_gid('conditions_title')
    fig.add_artist(title)

    return title


def _reduced_chi2(grid, model, n_par, window=None):
    """Reduced chi2 of a model spectrum on the fitted pixels, optionally only on those inside `window` (boolean pixel
    mask). n_dof = fitted pixels - n_par (NaN if there are none left)."""
    good = grid.mask & np.isfinite(grid.lum) & (grid.lum_err > 0)
    if window is not None:
        good = good & window
    chi2 = np.sum(((grid.lum[good] - model[good]) / grid.lum_err[good]) ** 2)
    n_dof = int(good.sum()) - n_par
    return float(chi2 / n_dof) if n_dof > 0 else np.nan


def _interval_chi2(grid, model, w3, w4):
    """Reduced chi2 on the fitted pixels of [w3, w4]: chi2 / pixels, since the free parameters are constrained by the
    whole spectrum and not by a single line."""
    return _reduced_chi2(grid, model, 0, window=(grid.wave >= w3) & (grid.wave <= w4))


def _line_label(ax, name, chi2_fit, chi2_mouse=None, mouse_color='magenta', zoom_factor=0.8):

    """
    Box at the lower right of a line axis (it looks like a legend, but legends cannot color parts of their text): the
    line name followed by the reduced chi2 of the fit in the color of the fitted spectrum and, if there is a mouse
    selection, that of the selection in parentheses in its color. It replaces the previous box of the axis.
    """

    for old in [artist for artist in ax.artists if artist.get_gid() == 'line_label']:
        old.remove()

    plain, fit_color = rcParams['text.color'], theme.colors['fitted_SSP']
    pieces = [(f'{name},  ', plain), (r'$\chi^2_{red}$ = ', plain), (f'{chi2_fit:.2f}', fit_color)]
    if chi2_mouse is not None:
        pieces += [(' (', plain), (f'{chi2_mouse:.2f}', mouse_color), (')', plain)]

    size = _fontsize('legend.fontsize', zoom_factor)
    box = HPacker(children=[TextArea(text, textprops=dict(color=color, fontsize=size)) for text, color in pieces],
                  align='center', pad=0, sep=0)

    label = AnchoredOffsetbox(loc='lower right', child=box, pad=0.3, borderpad=0.4, frameon=True)
    edge = rcParams['legend.edgecolor']
    label.patch.set(facecolor=ax.get_facecolor(), edgecolor=rcParams['axes.edgecolor'] if edge == 'inherit' else edge,
                    alpha=1)
    label.patch.set_boxstyle('round', pad=0, rounding_size=0.2)
    label.set_gid('line_label')
    label.set_zorder(1000)  # artists are drawn by zorder: the lines go up to zorder=600
    ax.add_artist(label)

    return label


def _connect_draw_picker(fig, ax0, ax1, mat, theta, grid, best_theta, color='magenta', pick_radius=10, rc=None,
                         line_axes=(), lines=None, best_model=None, zoom_factor=1, z_sun=None, plot_coords=None):

    """
    Left click on an off-diagonal scatter-matrix panel overplots in `ax0` the model of the nearest posterior draw
    (all four parameters of that draw, with the same nearest-node lookup as the likelihood) labelled with its reduced
    chi2 and marks the draw on every panel. Its residuals go in `ax1` or, if that panel was replaced by the `line_axes`,
    its model is overplotted in every line axis. The title over `ax1` (or over the line axes) adds the physical
    conditions of the draw, in parentheses, to those of the best fit `best_theta`. Right click clears the selection.
    `mat` is the (n, n) axes array of scatter_matrix and `theta` the (chain, draw, param) posterior array, both in the
    order of `best_theta` (the sampled values, i.e. log10(Z)). `plot_coords` maps sampled values to the plotted
    coordinates (Z/Z_sun when the metallicity is shown in solar units) and `z_sun` is the solar value for the title.
    Clicks farther than `pick_radius` pixels from any draw are ignored.
    """

    # Lower triangle of the scatter matrix: row i = y parameter, column j = x parameter
    n = mat.shape[0]
    panels = {(i, j): mat[i, j] for i in range(n) for j in range(i + 1)}
    flat = theta.reshape(-1, theta.shape[-1])
    plot_coords = (lambda t: np.asarray(t, dtype=float)) if plot_coords is None else plot_coords
    flat_plot = plot_coords(flat)

    spec_axes = [ax for ax in (ax0, ax1, *line_axes) if ax is not None]
    title_axes = list(line_axes) if len(line_axes) > 0 else [ax1]

    # Line zooms: name, interval and reduced chi2 of the fit on it
    intervals = [] if len(line_axes) == 0 else list(lines[['w3', 'w4']].itertuples(name=None))
    fit_chi2 = [_interval_chi2(grid, best_model, w3, w4) for _, w3, w4 in intervals]

    selected = []

    def clear():
        while selected:
            selected.pop().remove()
        _conditions_title(fig, title_axes, best_theta, z_sun=z_sun)
        for ax, (name, _, _), chi2_fit in zip(line_axes, intervals, fit_chi2):
            _line_label(ax, name, chi2_fit, zoom_factor=zoom_factor)

    def on_click(event):

        # Clicks while zooming or panning with the toolbar are not selections
        toolbar = getattr(fig.canvas, 'toolbar', None)
        if (toolbar is not None) and (toolbar.mode != ''):
            return

        with rc_context(rc):

            # Right click: clear the selection
            if event.button == 3:
                clear()
                ax0.legend(loc=1)
                fig.canvas.draw_idle()
                return

            key = next((k for k, ax in panels.items() if ax is event.inaxes), None)
            if (event.button != 1) or (key is None) or (key[0] == key[1]):
                return

            # Nearest draw in display pixels, so the parameter scales do not matter
            i, j = key
            pix = event.inaxes.transData.transform(flat_plot[:, [j, i]])
            dist2 = np.sum((pix - (event.x, event.y)) ** 2, axis=1)
            k = int(np.argmin(dist2))
            if dist2[k] > pick_radius ** 2:
                return

            theta_k = flat[k]
            plot_k = plot_coords(theta_k)
            row = int(grid.rows(theta_k[0], theta_k[1])[0])
            if row < 0:
                chain, draw = np.unravel_index(k, theta.shape[:2])
                _logger.warning(f'Draw (chain {chain}, draw {draw}) at log(age)={theta_k[0]:.3f}, '
                                f'log10(Z)={theta_k[1]:.3f} falls on an (age, Z) combination missing from the grid.')
                return
            model = grid.model(theta_k)

            # Freeze the current limits so the selection does not rescale the panels
            clear()
            for ax in (*spec_axes, *panels.values()):
                ax.set_autoscale_on(False)

            # Draw location on every panel
            for (pi, pj), ax in panels.items():
                if pi == pj:
                    selected.append(ax.axvline(plot_k[pi], color=color, lw=1.2, zorder=20, linestyle='--'))
                else:
                    selected.extend(ax.plot(plot_k[pj], plot_k[pi], 'o', mfc='none', mec=color, ms=8, mew=1.5,
                                            zorder=20))

            # Model of the draw with its reduced chi2 in the spectrum legend
            label = rf'Mouse selection ($\chi^2_{{red}}$={_reduced_chi2(grid, model, theta.shape[-1]):.2f})'
            selected.extend(ax0.step(grid.wave, model, where='mid', color=color, lw=1, label=label, zorder=600, linestyle='--'))
            ax0.legend(loc=1)

            # Residuals of the draw
            if ax1 is not None:
                selected.extend(ax1.step(grid.wave, _relative_residuals(grid.lum, model), where='mid', color=color,
                                         lw=1, zorder=600))

            # Model of the draw on every line interval, with its reduced chi2 on it next to that of the fit
            for ax, (name, w3, w4), chi2_fit in zip(line_axes, intervals, fit_chi2):
                selected.extend(ax.step(grid.wave, model, where='mid', color=color, lw=1, zorder=600, linestyle='--'))
                _line_label(ax, name, chi2_fit, _interval_chi2(grid, model, w3, w4), mouse_color=color,
                            zoom_factor=zoom_factor)

            # Physical conditions of the draw (age and Z of the node it maps to, E(B-V) and log(A) of the draw itself)
            mouse = (grid.age_node[row], grid.z_node[row], theta_k[2], theta_k[3])
            _conditions_title(fig, title_axes, best_theta, mouse, mouse_color=color, z_sun=z_sun)

            fig.canvas.draw_idle()

    return fig.canvas.mpl_connect('button_press_event', on_click)


def scatter_matrix(fig, slot, samples, names, pct, p0=None, p0_arr=None, ref=None, labels=None, solar=False,
                   nodes=None):

    """Corner-style scatter matrix nested in a gridspec slot: histograms on the diagonal, scatter in the lower triangle.
    Each diagonal panel is titled with the parameter median and its 16-84% uncertainties. Returns the (n, n) axes
    array (None above the diagonal).

    fig     : parent figure
    slot    : SubplotSpec to fill (e.g. gs[0, 1])
    samples : {name: 1D posterior samples}
    pct     : {'p16'|'p50'|'p84': {name: value}} marginal percentiles
    p0      : {name: value} initial sampler values (optional)
    p0_arr  : {name: array} bootstrap initial values (optional)
    ref     : {name: value} true values (optional)
    labels  : {name: LaTeX label} axis labels (default: latex_notation)
    solar   : True if the log_Z entries hold Z/Z_sun (linear) instead of log10(Z)
    nodes   : {name: 1D values} grid node positions drawn as thin lines on the panels of that parameter (optional)
    """

    labels = plot_vars['latex'] if labels is None else labels
    n = len(names)
    sub = slot.subgridspec(n, n, wspace=0.25, hspace=0.08)
    mat = np.empty((n, n), dtype=object)
    title_size, tick_size = _fontsize('axes.titlesize', 0.8), _fontsize('xtick.labelsize', 0.5)
    c_sample, c_band, c_fit = theme.colors['sample_value'], theme.colors['predictive_band'], theme.colors['fitted_SSP']
    c_p0, c_ref = theme.colors['initial_value'], theme.colors['true_value']

    for i in range(n):
        for j in range(i + 1):

            # Each column shares x with its histogram, each row of scatter panels shares y
            ax = fig.add_subplot(sub[i, j], sharex=mat[j, j] if i > j else None,
                                 sharey=mat[i, 0] if 0 < j < i else None)
            mat[i, j] = ax
            x_name, y_name = names[j], names[i]

            if i == j:
                ax.hist(samples[x_name], bins=30, color=c_sample, alpha=0.7)
                ax.axvspan(pct['p16'][x_name], pct['p84'][x_name], color=c_band, alpha=0.3)
                ax.axvline(pct['p50'][x_name], color=c_fit, lw=1.5)
                if p0 is not None: ax.axvline(p0[x_name], color=c_p0, ls='--', lw=1.5)
                if ref is not None: ax.axvline(ref[x_name], color=c_ref, ls='-', lw=1.5)
                ax.set_yticks([])
                ax.set_title(_value_label(x_name, pct['p16'][x_name], pct['p50'][x_name], pct['p84'][x_name], solar=solar),
                             fontsize=title_size)
            else:
                ax.scatter(samples[x_name], samples[y_name], s=3, alpha=0.1, color=c_sample, rasterized=True)
                if p0_arr is not None:
                    ax.scatter(p0_arr[x_name], p0_arr[y_name], s=6, alpha=0.3, color=c_p0, lw=0, zorder=1.5,
                               rasterized=True)
                ax.axvspan(pct['p16'][x_name], pct['p84'][x_name], color=c_band, alpha=0.2)
                ax.axhspan(pct['p16'][y_name], pct['p84'][y_name], color=c_band, alpha=0.2)
                ax.axvline(pct['p50'][x_name], color=c_fit, lw=1.5)
                ax.axhline(pct['p50'][y_name], color=c_fit, lw=1.5)
                if p0 is not None:
                    ax.axvline(p0[x_name], color=c_p0, ls='--', lw=0.8)
                    ax.axhline(p0[y_name], color=c_p0, ls='--', lw=0.8)
                if ref is not None:
                    ax.axvline(ref[x_name], color=c_ref, ls='-', lw=0.8)
                    ax.axhline(ref[y_name], color=c_ref, ls='-', lw=0.8)
                    ax.scatter(ref[x_name], ref[y_name], marker='*', s=120, color=c_ref, edgecolor='white', zorder=11)

            # Axis labels and tick labels only on the outer edges
            if i == n - 1:
                ax.set_xlabel(f"${labels.get(x_name, x_name)}$")
            else:
                ax.tick_params(labelbottom=False)
            if j == 0 and i > 0:
                ax.set_ylabel(f"${labels.get(y_name, y_name)}$")
            elif j > 0:
                ax.tick_params(labelleft=False)

            ax.tick_params(labelsize=tick_size)

    # Grid nodes inside the current limits (the limits are frozen so the lines never rescale a panel)
    if nodes:
        style = dict(color=theme.colors['fg'], lw=0.5, alpha=0.45, zorder=0.5)
        for i in range(n):
            for j in range(i + 1):
                ax = mat[i, j]
                xlim, ylim = ax.get_xlim(), ax.get_ylim()
                for v in nodes.get(names[j], ()):
                    if xlim[0] <= v <= xlim[1]: ax.axvline(v, **style)
                if i != j:
                    for v in nodes.get(names[i], ()):
                        if ylim[0] <= v <= ylim[1]: ax.axhline(v, **style)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

    return mat


def residuals_plot(ax, x, y, best_model, band=None):
    """Fractional residuals (y - model) / y on an existing axes, optionally with the predictive band as residuals.

    ax         : axes to draw into
    x, y       : wavelength and observed flux
    best_model : best-fit model flux
    band       : (p16, p50, p84) predictive arrays, or None
    """

    if band is not None:
        ax.fill_between(x, _relative_residuals(y, band[0]), _relative_residuals(y, band[2]), step='mid', alpha=0.25,
                        color=theme.colors['predictive_band'], ec=None, zorder=0)

    ax.step(x, _relative_residuals(y, best_model), where='mid', color=theme.colors['fg'], lw=1)
    ax.axhline(0, color=theme.colors['fg'], linestyle=':')

    ax.set_ylim(-0.7, 0.7)
    ax.set_ylabel('Residuals')
    ax.set_xlabel(r"Wavelength ($\AA$)")

    return ax


def trace_plots(fig, slot, chains, names, pct, p0=None, ref=None, labels=None):

    """Row of trace plots (one per parameter, all chains overlaid) nested in a gridspec slot. Returns the axes list.

    fig    : parent figure
    slot   : SubplotSpec to fill (e.g. gs[1, 1])
    chains : {name: array (chain, draw)} posterior samples
    pct    : {'p16'|'p50'|'p84': {name: value}} marginal percentiles
    p0     : {name: value} initial sampler values (optional)
    ref    : {name: value} true values (optional)
    """

    labels = plot_vars['latex'] if labels is None else labels
    tick_size = _fontsize('ytick.labelsize', 0.5)

    # The empty top row of the subgrid leaves room for the x labels of the plot above
    sub = slot.subgridspec(2, len(names), height_ratios=[0.5, 1], wspace=0.25, hspace=0)
    axes = []
    for k, name in enumerate(names):
        ax = fig.add_subplot(sub[1, k])
        ax.plot(chains[name].T, lw=0.3, alpha=0.25, color=theme.colors['sample_value'], rasterized=True)
        ax.axhspan(pct['p16'][name], pct['p84'][name], color=theme.colors['predictive_band'], alpha=0.3)
        ax.axhline(pct['p50'][name], color=theme.colors['fitted_SSP'], lw=1.5)
        if p0 is not None: ax.axhline(p0[name], color=theme.colors['initial_value'], ls='--', lw=1.5)
        if ref is not None: ax.axhline(ref[name], color=theme.colors['true_value'], ls='-', lw=1.5)
        ax.set_xlabel(f"${labels.get(name, name)}$")
        ax.tick_params(labelsize=tick_size)
        ax.set_xticks([])
        axes.append(ax)

    return axes


def line_zoom_axes(fig, slot, x, y, best_model, lines, grid, wspace=0.1, zoom_factor=1):
    """Row of axes nested in a gridspec slot, one per line, zoomed on the [w3, w4] interval of each line with the input
    and the fitted spectrum. The y axes have no ticks (the flux range comes from the input and the fit inside the
    interval) and the line name goes in a box at the lower right, with the reduced chi2 of the fit on that interval.

    fig        : parent figure
    slot       : SubplotSpec to fill (e.g. gs_l[1])
    x, y       : wavelength and observed flux
    best_model : best-fit model flux
    lines      : DataFrame with the line names as index and the w3, w4 columns (one axis per row)
    grid       : SSPGrid with the observation (flux, uncertainty and fitted-pixel mask) for the chi2
    wspace     : horizontal space between the axes
    """

    sub = slot.subgridspec(1, len(lines), wspace=wspace)
    tick_size = _fontsize('axes.labelsize', 0.8)

    axes = []
    for k, (name, w3, w4) in enumerate(lines[['w3', 'w4']].itertuples(name=None)):
        ax = fig.add_subplot(sub[0, k])
        ax.step(x, y, where='mid', color=theme.colors['fg'], lw=1)
        ax.step(x, best_model, where='mid', color=theme.colors['fitted_SSP'], lw=1, zorder=500)
        ax.set_xlim(w3, w4)

        # Flux range of the input and the fit inside the interval (the mouse selection does not rescale it)
        window = (x >= w3) & (x <= w4)
        flux = np.concatenate([y[window], best_model[window]])
        flux = flux[np.isfinite(flux)]
        if flux.size > 0:
            span = np.ptp(flux) or abs(flux.max()) or 1.0
            free_bottom, free_top = 0.30, 0.20
            full = span / (1 - free_bottom - free_top)
            ax.set_ylim(flux.min() - free_bottom * full, flux.max() + free_top * full)

        ax.set_yticks([])
        ax.locator_params(axis='x', nbins=3)
        ax.ticklabel_format(axis='x', useOffset=False, style='plain')
        ax.tick_params(labelsize=tick_size)

        # Line name and reduced chi2 of the fit on the interval, above the data
        _line_label(ax, name, _interval_chi2(grid, best_model, w3, w4), zoom_factor=zoom_factor)

        axes.append(ax)

    return axes


def ssp_fitting_infograph(trace, fname=None, truth_dict=None, title=None, show_stellar=True, fig_cfg=None,
                          pick_draws=True, residual_df=None, solar_Z=False, z_sun=Z_SUN, show_grid=False, **kwargs):

    """
    Infograph of an SSP fit: spectrum and fit, residuals (or line zooms), scatter matrix and traces.

    solar_Z   : show the metallicity as a solar fraction, Z/Z_sun, on the scatter-matrix and trace axes and in the
                titles, instead of log10(Z) and Z
    z_sun     : solar metallicity used by `solar_Z`
    show_grid : draw the age and metallicity grid nodes as thin lines on the scatter-matrix panels of those parameters
    """

    # Unpack the trace: everything comes from the stored inputs/outputs datasets
    inputs, out = trace.inputs, trace.outputs
    names = list(out.attrs['param_names'])

    x, y, mask = inputs['wave'].values, inputs['lum'].values, inputs['mask'].values.astype(bool)
    best_model = out['best_model'].values
    band = tuple(out[f'band_{q}'].values for q in ('p16', 'p50', 'p84')) if 'band_p16' in out else None

    # Visited grid nodes with the sampler's nearest-node lookup (best-fit components and draw picker)
    grid = ssp_lib.SSPGrid.from_idata(trace)

    # Stored per-parameter vectors as {name: value} dicts, keyed by the 'param' coordinate
    as_dict = lambda key: out[key].to_series().to_dict() if key in out else None
    best_theta = out['theta'].sel(param=names).values
    pct = {q: as_dict(f'marginal_{q}') for q in ('p16', 'p50', 'p84')}
    p0 = as_dict('p0')
    p0_arr = dict(zip(names, out['p0_arr'].values.T)) if 'p0_arr' in out else None

    # True values: the argument takes precedence over those stored in the trace
    if truth_dict is not None:
        ref = {'log_age': float(truth_dict['log_age']), 'log_Z': float(np.log10(truth_dict['metallicity'])),
               'ebv': float(truth_dict['ebv']), 'log_A': float(np.log10(truth_dict['mass'] / 1e6))}
    else:
        ref = as_dict('truth')

    raw_chains = {k: trace.posterior[k].values for k in names}  # (chain, draw), as sampled

    # Metallicity in solar units: every plotted quantity uses Z/Z_sun instead of the sampled log10(Z), while the model
    # lookups (draw picker, best fit) keep using the sampled values
    z_sun = z_sun if solar_Z else None
    z_col = names.index('log_Z')

    def plot_coords(theta):
        """Sampled values (..., param) -> plotted coordinates."""
        theta = np.array(theta, dtype=float)
        if z_sun is not None:
            theta[..., z_col] = 10.0 ** theta[..., z_col] / z_sun
        return theta

    def plotted(d):
        """{name: value or array} of sampled values -> plotted values (None stays None)."""
        if (d is None) or (z_sun is None):
            return d
        return {k: (10.0 ** np.asarray(v) / z_sun if k == 'log_Z' else v) for k, v in d.items()}

    pct = {q: plotted(v) for q, v in pct.items()}
    p0, p0_arr, ref = plotted(p0), plotted(p0_arr), plotted(ref)
    chains = plotted(raw_chains)
    samples = {k: v.ravel() for k, v in chains.items()}
    labels = _labels(z_sun)
    nodes = {'log_age': grid.age_log, 'log_Z': plotted({'log_Z': grid.z_log})['log_Z'] if z_sun else grid.z_log} \
        if show_grid else None

    # Observed spectrum as a LiMe object (stored mask: True = fitted, LiMe pixel_mask: True = masked)
    spec = lime.Spectrum(x, y, inputs['lum_err'].values, redshift=0, pixel_mask=~mask,
                         units_wave=inputs.attrs['units_wave'], units_flux=inputs.attrs['units_flux'])

    # Line list: up to 4 lines inside the spectrum replace the residuals panel by a zoom on their [w3, w4] intervals
    if residual_df is not None:
        residual_df = spec.retrieve.lines_frame(ref_bands=residual_df, adjust_central_band=False, exclude_bands_masked=False)
        if (residual_df is None) or (len(residual_df) == 0):
            _logger.warning('None of the residual_df lines is inside the spectrum range: the residuals panel is '
                            'shown instead.')
            residual_df = None
        else:
            residual_df = residual_df.iloc[:4]
    zoom_factor = 1 if residual_df is None else min(1.0, 2.75 / len(residual_df))

    # Adjust the default theme
    plt_cfg = theme.fig_defaults(fig_cfg, fig_type='ssp_infograph')
    ax_cfg = None if title is None else {'title': title}

    # Generate the figure
    with rc_context(plt_cfg):

        fig = plt.figure()
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], width_ratios=[2, 1.3])
        ax0 = fig.add_subplot(gs[0, 0])

        # Lower-left slot: residuals panel or one zoom axis per line (top third left free for the conditions title)
        slot_low = gs[1, 0].subgridspec(2, 1, height_ratios=[0.5, 1], hspace=0)[1]
        if residual_df is None:
            ax1, line_axes = fig.add_subplot(slot_low, sharex=ax0), []
        else:
            ax1 = None
            line_axes = line_zoom_axes(fig, slot_low, x, y, best_model, residual_df, grid, zoom_factor=zoom_factor)
        spec_axes = [ax for ax in (ax0, ax1, *line_axes) if ax is not None]
        title_axes = line_axes if len(line_axes) > 0 else [ax1]

        # Input spectrum (the wavelength axis is only hidden if the residuals panel is below)
        spec.plot.spectrum(label='Input spectrum', show_masks=False, in_ax=ax0, fig_cfg=plt_cfg, title=title, **kwargs)

        if ax1 is not None:
            ax0.tick_params(labelbottom=False)
            ax0.set_xlabel('')

        # Shade masked-out regions at the pixel edges
        for w_lo, w_hi in _mask_spans(x, mask):
            for ax in spec_axes:
                ax.axvspan(w_lo, w_hi, alpha=0.3, color=theme.colors['masked_region'], lw=0)

        chi2_red = _reduced_chi2(grid, best_model, len(names))
        ax0.step(x, best_model, where='mid', color=theme.colors['fitted_SSP'], lw=1, zorder=500,
                 label=rf'Fitted spectrum ($\chi^2_{{red}}$={chi2_red:.2f})')

        # Stellar / Nebular components of the best fit
        star, neb = grid.components(best_theta)
        if neb is not None:
            if show_stellar:
                ax0.step(x, star, where='mid', color='darkorange', lw=0.8, label='Stellar', linestyle=':', zorder=400)
            ax0.step(x, neb, where='mid', color='crimson', lw=0.8, label='Nebular continuum', linestyle=':', zorder=400)

        # Predictive band
        if band is not None:
            ax0.fill_between(x, band[0], band[2], step='mid', alpha=0.6, color=theme.colors['predictive_band'], ec=None,
                             label='Predictive 16-84%', zorder=0)
        ax0.legend(loc=1)

        # Residuals panel
        if ax1 is not None:
            residuals_plot(ax1, x, y, best_model, band=band)

        # Scatter matrix and traces
        mat = scatter_matrix(fig, gs[0, 1], samples, names, pct, p0=p0, p0_arr=p0_arr, ref=ref, labels=labels,
                             solar=z_sun is not None, nodes=nodes)
        trace_axes = trace_plots(fig, gs[1, 1], chains, names, pct, p0=p0, ref=ref, labels=labels)

        fig.subplots_adjust(left=0.06, right=0.98, top=0.92 if title is not None else 0.95,
                            bottom=0.09, wspace=0.15, hspace=0 if ax1 is not None else 0.20)

        # Titles placed after subplots_adjust, since they use the final axes positions
        box = [ax.get_position() for ax in trace_axes]
        fig.text((box[0].x0 + box[-1].x1) / 2, max(b.y1 for b in box) + 0.01, 'Traces', ha='center', va='bottom',
                 fontsize=_fontsize('axes.titlesize', 0.8))
        _conditions_title(fig, title_axes, best_theta, z_sun=z_sun)

        # Click a draw in the scatter matrix to overplot its model (needs an interactive backend)
        if pick_draws and (fname is None):
            theta = np.stack([raw_chains[k] for k in names], axis=-1)
            fig._draw_picker_cid = _connect_draw_picker(fig, ax0, ax1, mat, theta, grid, best_theta, rc=plt_cfg,
                                                        line_axes=line_axes, lines=residual_df,
                                                        best_model=best_model, zoom_factor=zoom_factor,
                                                        z_sun=z_sun, plot_coords=plot_coords)

        if fname is not None:
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

    return fig



# ---------------------------------------------------------------------------------------------------------------------
# plots.py: two import changes
#   from pandas import Series, DataFrame
#   from specsy.io import SpecSyError, specsy_cfg, load_ssp_trace_results, load_trace
# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# plots.py: two import changes
#   from pandas import Series, DataFrame
#   from specsy.io import SpecSyError, specsy_cfg, load_ssp_trace_results, load_trace
# ---------------------------------------------------------------------------------------------------------------------


def ssp_line_grid(traces, line, bands, fname=None, n_cols=4, title=None, fig_cfg=None):

    """
    Grid with one panel per object zoomed on the [w3, w4] interval of the same line: input spectrum and best SSP fit,
    with the masked regions shaded and the reduced chi2 of the fit on the interval (as in line_zoom_axes). Each panel
    has its own interval and flux range, so the y axes have no ticks.

    traces  : {object name: trace or trace file}, in the panel order. With files, each trace is loaded only while its
              panel is drawn
    line    : line label in the bands frames (e.g. 'C4_1548A'), also the default figure title
    bands   : {object name: bands DataFrame} with `line` in the index and the w3, w4 columns (e.g. the residual_df of
              each object from spec.retrieve.lines_frame), or a single DataFrame or (w3, w4) interval for every object
    fname   : output file. If None the figure is shown on screen
    n_cols  : panels per row (the rows follow from the number of objects)
    title   : figure title (default: `line`)
    fig_cfg : matplotlib rcParams overrides
    """

    n_objs = len(traces)
    n_rows = int(np.ceil(n_objs / n_cols))

    # Set the plot format where the user's overwrites the default
    size_conf = {'figure.figsize': (3 * n_cols, 2.2 * n_rows)}
    size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}
    plt_cfg = theme.fig_defaults(size_conf)

    c_fg, c_fit, c_mask = theme.colors['fg'], theme.colors['fitted_SSP'], theme.colors['masked_region']
    no_band, no_coverage = [], []

    with rc_context(plt_cfg):

        fig, axes = plt.subplots(n_rows, n_cols, squeeze=False, layout='constrained')
        axes = axes.ravel()

        for ax, (obj, trace) in zip(axes, traces.items()):

            ax.set_title(obj)
            ax.set_yticks([])

            # Line interval of the object
            band = bands.get(obj) if isinstance(bands, dict) else bands
            if isinstance(band, DataFrame):
                band = band.loc[line, ['w3', 'w4']].to_numpy() if line in band.index else None
            if band is None:
                no_band.append(obj)
                ax.set_xticks([])
                ax.text(0.5, 0.5, f'No {line} band', transform=ax.transAxes, ha='center', va='center')
                continue
            w3, w4 = band

            # Observation and best fit from the stored inputs/outputs groups. The namespace has the attributes
            # _interval_chi2 reads from an SSPGrid (wave, lum, lum_err, mask), so the grid is not rebuilt
            trace = load_trace(trace) if isinstance(trace, (str, Path)) else trace
            inputs = trace.inputs
            obs = SimpleNamespace(wave=inputs['wave'].values, lum=inputs['lum'].values,
                                  lum_err=inputs['lum_err'].values, mask=inputs['mask'].values.astype(bool))
            best_model = trace.outputs['best_model'].values

            ax.step(obs.wave, obs.lum, where='mid', color=c_fg, lw=1)
            ax.step(obs.wave, best_model, where='mid', color=c_fit, lw=1, zorder=500)
            for w_lo, w_hi in _mask_spans(obs.wave, obs.mask):
                ax.axvspan(w_lo, w_hi, alpha=0.3, color=c_mask, lw=0)

            ax.set_xlim(w3, w4)
            ax.locator_params(axis='x', nbins=3)
            ax.ticklabel_format(axis='x', useOffset=False, style='plain')

            # Flux range of the input and the fit inside the interval, with room at the bottom for the chi2
            window = (obs.wave >= w3) & (obs.wave <= w4)
            ylim = _zoom_ylim(obs.lum[window], best_model[window])
            if ylim is None:
                no_coverage.append(obj)
                continue
            ax.set_ylim(*ylim)

            chi2 = _interval_chi2(obs, best_model, w3, w4)
            ax.text(0.97, 0.05, rf'$\chi^2_{{red}}$ = {chi2:.2f}', transform=ax.transAxes, ha='right', va='bottom',
                    color=c_fit, fontsize=_fontsize('legend.fontsize', 0.8))

        if no_band:
            _logger.warning(f'{len(no_band)} of {n_objs} object(s) without a {line} band in their bands frame: '
                            f'{no_band}')
        if no_coverage:
            _logger.warning(f'{len(no_coverage)} of {n_objs} object(s) without finite flux in their {line} '
                            f'interval: {no_coverage}')

        # Remove the empty panels
        for k in range(n_objs, n_rows * n_cols):
            axes[k].remove()

        fig.supxlabel(r'Wavelength ($\AA$)')
        fig.suptitle(line if title is None else title)
        fig.legend(handles=[Line2D([], [], color=c_fg, label='Input spectrum'),
                            Line2D([], [], color=c_fit, label='Fitted spectrum')],
                   loc='outside upper right', ncol=2, frameon=False)

        if fname is not None:
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

    return fig


# ---------------------------------------------------------------------------------------------------------------------
# 4_SPECSY_BPASS.py (import ssp_line_grid from specsy.plotting.plots)
# ---------------------------------------------------------------------------------------------------------------------

# residual_dict = {}
# for i, obj in enumerate(target_list):
#     ...
#     residual_df = spec_treated.retrieve.lines_frame(**cfg_sample['SSP_review'][obj]['residuals_df'])
#     residual_dict[obj] = residual_df
#     ...
#
# trace_files = {obj: SpecSy_folder / f'{obj}_{label_fit}_trace' for obj in target_list}
# ssp_line_grid(trace_files, 'C4_1548A', residual_dict, fname=SpecSy_folder / f'C4_1548A_grid{label_fit}.png')


def plot_bad_pixels(wave_arr, flux_arr, err_arr, mask_arr, nonfinite, neg_err):

    fig, (ax_f, ax_e) = plt.subplots(2, 1, sharex=True, figsize=(10, 6))

    # Only pixels with a finite wavelength can be placed on the x-axis
    wave_ok = np.isfinite(wave_arr)
    used = mask_arr & wave_ok
    masked = ~mask_arr & wave_ok

    for ax, arr, label in ((ax_f, flux_arr, 'Flux'), (ax_e, err_arr, 'Uncertainty')):
        # Valid pixels and masked pixels
        ax.step(wave_arr[used], arr[used], where='mid', color='tab:blue', lw=0.8, label='Unmasked')
        ax.plot(wave_arr[masked], arr[masked], '.', color='lightgrey', ms=2, label='Masked')

        # Non-finite pixels: no y-value to plot, so mark them with ticks along the bottom axis
        bad_nf = nonfinite[label] & wave_ok
        ax.plot(wave_arr[bad_nf], np.zeros(bad_nf.sum()), '|', color='red', ms=15, mew=1.5,
                transform=ax.get_xaxis_transform(), label=f'Non-finite ({nonfinite[label].sum()})')
        ax.set_ylabel(label)

    # Negative uncertainties
    bad_neg = neg_err & wave_ok
    ax_e.plot(wave_arr[bad_neg], err_arr[bad_neg], 'o', color='orange', mfc='none', ms=6,
              label=f'Negative error ({neg_err.sum()})')
    ax_e.axhline(0, color='k', ls=':', lw=0.8)

    # Pixels that cannot be drawn because their wavelength is not finite
    n_no_wave = nonfinite['Wavelength'].sum()
    if n_no_wave:
        ax_f.set_title(f'{n_no_wave} unmasked pixel(s) with non-finite wavelength cannot be shown')

    ax_e.set_xlabel('Wavelength')
    ax_f.legend(loc='upper right', fontsize=8)
    ax_e.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    plt.show()
    return


def _bk_color(key):
    """Theme color as hex, so matplotlib names (e.g. 'tab:blue') also work in Bokeh."""
    return colors.to_hex(theme.colors[key])


def _bk_style(fig_cfg):
    """Bokeh theme configuration without the size, layout, axis type and toolbar entries, which are set per panel."""
    cfg = theme.fig_defaults(fig_cfg, plot_lib='bokeh')
    return {k: v for k, v in cfg.items() if (k not in plot_vars['bokeh']['_BK_STYLE_SKIP']) and not k.startswith(('frame_', 'min_', 'max_'))}


def _box(width, height, left, right, top, bottom):
    """Fixed outer size of a panel and the minimum borders around its frame."""
    return dict(width=int(width), height=int(height), sizing_mode='fixed', margin=(0, 0, 0, 0),
                min_border_left=left, min_border_right=right, min_border_top=top, min_border_bottom=bottom)


def _fixed(layout_fn, *children, width, height):
    """bokeh row / column with a fixed size, so Bokeh does not reflow or stretch its children."""
    return layout_fn(*children, sizing_mode='fixed', width=int(width), height=int(height))


def _spacer(width, height):
    return Spacer(width=int(width), height=int(height), sizing_mode='fixed', margin=(0, 0, 0, 0))


def _bk_figure(style, toolbar=None, **kwargs):
    """Figure with the default tools and the theme style applied."""
    fig = figure(tools=style.get('tools', plot_vars['bokeh']['BOKEH_TOOLS']), toolbar_location=toolbar, **kwargs)
    update_bokeh_figure(fig, style)
    fig.update(**{k: v for k, v in kwargs.items() if k in plot_vars['bokeh']['_BK_LAYOUT_KEYS']})   # the theme must not undo the layout
    return fig


def _bk_span(fig, value, color, dimension='height', dash='solid', width=1.5, alpha=1.0, level='glyph', visible=True):
    """Vertical (dimension='height') or horizontal (dimension='width') line, as axvline/axhline."""
    span = Span(location=value, dimension=dimension, line_color=color, line_dash=dash, line_width=width,
                line_alpha=alpha, level=level, visible=visible)
    fig.add_layout(span)
    return span


def _bk_band(fig, low, high, color, alpha, vertical=True):
    """Shaded interval behind the data, as axvspan (vertical) or axhspan."""
    limits = dict(left=low, right=high) if vertical else dict(bottom=low, top=high)
    fig.add_layout(BoxAnnotation(**limits, fill_color=color, fill_alpha=alpha, line_color=None, level='underlay'))
    return


def _bk_mask_shading(figs, x, mask):
    """Masked-out regions shaded at the pixel edges on every spectral panel."""
    c_mask = _bk_color('masked_region')
    for w_lo, w_hi in _mask_spans(x, mask):
        for fig in figs:
            _bk_band(fig, w_lo, w_hi, c_mask, 0.3)
    return


def _in_view(values, data, pad=0.1):
    """Entries of `values` inside the range of `data` widened by `pad` of its span (grid nodes worth drawing)."""
    if values is None:
        return []
    lo, hi = np.nanmin(data), np.nanmax(data)
    margin = pad * (hi - lo)
    values = np.asarray(values, dtype=float)
    return values[(values >= lo - margin) & (values <= hi + margin)]


def _zoom_ylim(y, model, free_bottom=0.30, free_top=0.20):
    """Flux range of the input and the fit inside a line interval, as in line_zoom_axes (None if nothing finite)."""
    flux = np.concatenate([y, model])
    flux = flux[np.isfinite(flux)]
    if flux.size == 0:
        return None
    span = np.ptp(flux) or abs(flux.max()) or 1.0
    full = span / (1 - free_bottom - free_top)
    return flux.min() - free_bottom * full, flux.max() + free_top * full


def _ssp_trace_data(trace, truth_dict=None, solar_Z=False, z_sun=Z_SUN, show_grid=False):

    """
    Everything the SSP infograph plots, unpacked from the trace inputs/outputs groups (same steps as in
    ssp_fitting_infograph). With `solar_Z` the plotted quantities use Z/Z_sun instead of the sampled log10(Z), while
    `raw_chains`, `best_theta` and `grid` keep the sampled values for the model lookups.
    """

    inputs, out = trace.inputs, trace.outputs
    names = list(out.attrs['param_names'])
    as_dict = lambda key: out[key].to_series().to_dict() if key in out else None

    z_sun = z_sun if solar_Z else None

    def plotted(d):
        """{name: value or array} of sampled values -> plotted values (None stays None)."""
        if (d is None) or (z_sun is None):
            return d
        return {k: (10.0 ** np.asarray(v) / z_sun if k == 'log_Z' else v) for k, v in d.items()}

    grid = SSPGrid.from_idata(trace)
    raw_chains = {k: trace.posterior[k].values for k in names}  # (chain, draw), as sampled
    chains = plotted(raw_chains)

    # True values: the argument takes precedence over those stored in the trace
    ref = truth_to_params(truth_dict) if truth_dict is not None else as_dict('truth')
    p0_arr = dict(zip(names, out['p0_arr'].values.T)) if 'p0_arr' in out else None

    return SimpleNamespace(
        inputs=inputs, names=names, grid=grid, z_sun=z_sun,
        x=inputs['wave'].values, y=inputs['lum'].values, mask=inputs['mask'].values.astype(bool),
        best_model=out['best_model'].values, best_theta=out['theta'].sel(param=names).values,
        band=tuple(out[f'band_{q}'].values for q in ('p16', 'p50', 'p84')) if 'band_p16' in out else None,
        pct={q: plotted(as_dict(f'marginal_{q}')) for q in ('p16', 'p50', 'p84')},
        p0=plotted(as_dict('p0')), p0_arr=plotted(p0_arr), ref=plotted(ref),
        raw_chains=raw_chains, chains=chains, samples={k: v.ravel() for k, v in chains.items()},
        labels=_labels(z_sun),
        nodes=plotted({'log_age': grid.age_log, 'log_Z': grid.z_log}) if show_grid else None)


def _condition_parts(fitted, chi2, z_sun=None):

    """
    HTML pieces of the conditions title (one per parameter plus the reduced chi2) with the fitted value in the color of
    the fitted spectrum, the decimals of each value and the metallicity normalization. Bokeh version of the pieces
    built in _conditions_title: the draw picker appends the values of the selected draw in parentheses.
    """

    labels, formats = plot_vars['bokeh']['CONDITION_LABELS'], plot_vars['bokeh']['CONDITION_FORMATS']
    if z_sun is not None:
        labels[1], formats[1] = 'Z/Z<sub>⊙</sub>', '.3f'
    labels.append('χ²<sub>red</sub>')
    formats.append('.2f')

    z_norm = 1.0 if z_sun is None else z_sun
    values = (fitted[0], 10.0 ** fitted[1] / z_norm, fitted[2], fitted[3], chi2)
    color = _bk_color('fitted_SSP')
    parts = [f'{label} = <span style="color:{color}">{format(v, fmt)}</span>'
             for label, v, fmt in zip(labels, values, formats)]

    return parts, [int(fmt[1:-1]) for fmt in formats], z_norm


def scatter_matrix_bokeh(src, samples, names, pct, style, p0=None, p0_arr=None, ref=None, labels=None, solar=False,
                         nodes=None, cell=120, pick_color='magenta'):

    """
    Bokeh version of scatter_matrix. Returns the (n, n) list of figures (None above the diagonal) and the hidden
    vertical lines of the diagonal panels that mark the selected draw.

    src     : ColumnDataSource with the (plotted) draws shown in the scatter panels, one column per parameter. It is
              shared by every panel, so a selected draw is highlighted in all of them
    samples : {name: 1D posterior samples} for the histograms
    cell    : side of each panel frame in pixels
    The other arguments behave as in scatter_matrix.
    """

    labels = plot_vars['latex'] if labels is None else labels     # NOTE: name clash, see the message above
    n = len(names)
    c_sample, c_band, c_fit = _bk_color('sample_value'), _bk_color('predictive_band'), _bk_color('fitted_SSP')
    c_p0, c_ref, c_fg = _bk_color('initial_value'), _bk_color('true_value'), _bk_color('fg')

    widths, heights = _matrix_sizes(n, cell)
    mat, markers = [[None] * n for _ in range(n)], []
    for i in range(n):
        for j in range(i + 1):

            x_name, y_name = names[j], names[i]

            # Each column shares x with its histogram, each row of scatter panels shares y. The fixed sizes and
            # borders keep the frames aligned whether a panel shows its axes and title or not
            left = plot_vars['bokeh']['BORDER_LEFT'] if j == 0 else plot_vars['bokeh']['BORDER_INNER']
            bottom = plot_vars['bokeh']['BORDER_BOTTOM'] if i == n - 1 else plot_vars['bokeh']['BORDER_MID']
            fig_kw = dict(output_backend='webgl', **_box(widths[j], heights[i], left, plot_vars['bokeh']['BORDER_RIGHT'],
                                                         plot_vars['bokeh']['BORDER_TOP'], bottom))
            if i > j:
                fig_kw['x_range'] = mat[j][j].x_range
            if 0 < j < i:
                fig_kw['y_range'] = mat[i][0].y_range

            if i == j:
                # _value_label gives two lines for matplotlib: one line here
                value_label = _value_label(x_name, pct['p16'][x_name], pct['p50'][x_name], pct['p84'][x_name],
                                           solar=solar)
                fig_kw['title'] = '$$' + value_label.split('\n')[-1].strip('$') + '$$'

            fig = _bk_figure(style, **fig_kw)
            mat[i][j] = fig

            if i == j:
                counts, edges = np.histogram(samples[x_name], bins=30)
                fig.quad(left=edges[:-1], right=edges[1:], bottom=0, top=counts, fill_color=c_sample, fill_alpha=0.7,
                         line_color=None)
                fig.y_range.start = 0
                _bk_band(fig, pct['p16'][x_name], pct['p84'][x_name], c_band, 0.3)
                _bk_span(fig, pct['p50'][x_name], c_fit)
                if p0 is not None: _bk_span(fig, p0[x_name], c_p0, dash='dashed')
                if ref is not None: _bk_span(fig, ref[x_name], c_ref)
                markers.append(_bk_span(fig, pct['p50'][x_name], pick_color, width=1.2, level='overlay', visible=False))
                fig.title.text_font_size = '8pt'
                fig.title.align = 'center'

            else:
                # The selected draw is drawn as an open circle; the other draws keep their look
                rend = fig.scatter(x_name, y_name, source=src, size=3, color=c_sample, alpha=0.1,
                                   nonselection_color=c_sample, nonselection_alpha=0.1)
                rend.selection_glyph = Scatter(x=x_name, y=y_name, size=9, fill_alpha=0, line_color=pick_color,
                                               line_width=1.5)
                if p0_arr is not None:
                    fig.scatter(p0_arr[x_name], p0_arr[y_name], size=4, color=c_p0, alpha=0.3, line_color=None)
                _bk_band(fig, pct['p16'][x_name], pct['p84'][x_name], c_band, 0.2)
                _bk_band(fig, pct['p16'][y_name], pct['p84'][y_name], c_band, 0.2, vertical=False)
                _bk_span(fig, pct['p50'][x_name], c_fit)
                _bk_span(fig, pct['p50'][y_name], c_fit, dimension='width')
                if p0 is not None:
                    _bk_span(fig, p0[x_name], c_p0, dash='dashed', width=0.8)
                    _bk_span(fig, p0[y_name], c_p0, dimension='width', dash='dashed', width=0.8)
                if ref is not None:
                    _bk_span(fig, ref[x_name], c_ref, width=0.8)
                    _bk_span(fig, ref[y_name], c_ref, dimension='width', width=0.8)
                    fig.scatter([ref[x_name]], [ref[y_name]], marker='star', size=14, color=c_ref, line_color='white')

            # Grid nodes inside the posterior range (Spans never rescale a panel)
            if nodes:
                for v in _in_view(nodes.get(x_name), samples[x_name]):
                    _bk_span(fig, v, c_fg, width=0.5, alpha=0.45, level='underlay')
                if i != j:
                    for v in _in_view(nodes.get(y_name), samples[y_name]):
                        _bk_span(fig, v, c_fg, dimension='width', width=0.5, alpha=0.45, level='underlay')

            # Axes only on the outer edges
            fig.xaxis.visible = i == n - 1
            fig.yaxis.visible = (j == 0) and (i > 0)
            fig.xaxis.axis_label = f'$${labels.get(x_name, x_name)}$$'
            fig.yaxis.axis_label = f'$${labels.get(y_name, y_name)}$$'
            for axis in (fig.xaxis, fig.yaxis):
                axis.major_label_text_font_size = '7pt'
                axis.axis_label_text_font_size = '9pt'
            fig.xaxis[0].ticker.desired_num_ticks = 3
            fig.yaxis[0].ticker.desired_num_ticks = 3

    return mat, markers


def _matrix_sizes(n, cell):
    """Outer widths of the scatter-matrix columns and outer heights of its rows (frames of `cell` pixels)."""
    widths = [(plot_vars['bokeh']['BORDER_LEFT'] if j == 0 else plot_vars['bokeh']['BORDER_INNER']) + cell + plot_vars['bokeh']['BORDER_RIGHT']
              for j in range(n)]
    heights = [plot_vars['bokeh']['BORDER_TOP'] + cell + (plot_vars['bokeh']['BORDER_BOTTOM'] if i == n - 1 else plot_vars['bokeh']['BORDER_MID'])
               for i in range(n)]
    return widths, heights


def _matrix_layout(mat, cell):
    """
    Scatter matrix as fixed-size rows (spacers above the diagonal) instead of gridplot, which can redistribute the
    height between its cells. Returns the layout and its outer width and height.
    """
    n = len(mat)
    widths, heights = _matrix_sizes(n, cell)
    rows = [_fixed(row, *[mat[i][j] if mat[i][j] is not None else _spacer(widths[j], heights[i]) for j in range(n)],
                   width=sum(widths), height=heights[i]) for i in range(n)]
    return _fixed(column, *rows, width=sum(widths), height=sum(heights)), sum(widths), sum(heights)


def trace_plots_bokeh(chains, names, pct, style, p0=None, ref=None, labels=None, size=96, height=None, max_draws=500):

    """
    Bokeh version of trace_plots: one figure per parameter with all chains overlaid. Each chain is thinned to at most
    `max_draws` draws to keep the page light. `size` is the frame width in pixels and `height` the outer height of
    each figure (a square frame by default). Returns the figure list.
    """

    height = plot_vars['bokeh']['BORDER_TOP'] + size + plot_vars['bokeh']['BORDER_BOTTOM'] if height is None else height

    labels = plot_vars['latex'] if labels is None else labels
    c_sample, c_band, c_fit = _bk_color('sample_value'), _bk_color('predictive_band'), _bk_color('fitted_SSP')
    c_p0, c_ref = _bk_color('initial_value'), _bk_color('true_value')

    figs = []
    for name in names:
        chain = chains[name]
        step = max(1, int(np.ceil(chain.shape[1] / max_draws)))
        draws = np.arange(chain.shape[1])[::step]

        fig = _bk_figure(style, output_backend='webgl', **_box(size + 55, height, 50, plot_vars['bokeh']['BORDER_RIGHT'],
                                                              plot_vars['bokeh']['BORDER_TOP'], plot_vars['bokeh']['BORDER_BOTTOM']))
        fig.multi_line([draws] * chain.shape[0], list(chain[:, ::step]), color=c_sample, alpha=0.25, line_width=0.5)
        _bk_band(fig, pct['p16'][name], pct['p84'][name], c_band, 0.3, vertical=False)
        _bk_span(fig, pct['p50'][name], c_fit, dimension='width')
        if p0 is not None: _bk_span(fig, p0[name], c_p0, dimension='width', dash='dashed')
        if ref is not None: _bk_span(fig, ref[name], c_ref, dimension='width')

        fig.xaxis.axis_label = f'$${labels.get(name, name)}$$'
        fig.xaxis.major_label_text_font_size = '0pt'
        fig.xaxis.major_tick_line_color = None
        fig.xaxis.minor_tick_line_color = None
        fig.yaxis.major_label_text_font_size = '7pt'
        fig.yaxis[0].ticker.desired_num_ticks = 3
        figs.append(fig)

    return figs


def residuals_plot_bokeh(fig, x, y, best_model, band=None):
    """Bokeh version of residuals_plot: fractional residuals (y - model) / y, optionally with the predictive band."""

    if band is not None:
        rend = fig.varea_step(x=x, y1=_relative_residuals(y, band[0]), y2=_relative_residuals(y, band[2]),
                              step_mode='center', fill_color=_bk_color('predictive_band'), fill_alpha=0.25)
        rend.level = 'underlay'

    fig.step(x, _relative_residuals(y, best_model), mode='center', color=_bk_color('fg'), line_width=1)
    _bk_span(fig, 0, _bk_color('fg'), dimension='width', dash='dotted', width=1)

    fig.yaxis.axis_label = 'Residuals'
    fig.xaxis.axis_label = 'Wavelength (Å)'

    return fig


def line_zoom_bokeh(x, y, best_model, lines, grid, style, width, height, sel_src=None, pick_color='magenta'):

    """
    Bokeh version of line_zoom_axes: one figure per line, zoomed on its [w3, w4] interval with the input and the fitted
    spectrum (and the selected draw from `sel_src`, if provided). The line name and the reduced chi2 of the fit on
    the interval go in the title. `width` is the total width of the figures and `height` their outer height. Returns
    the figures, their titles and their (w3, w4) intervals.
    """

    c_fg, c_fit = _bk_color('fg'), _bk_color('fitted_SSP')
    fig_widths = np.diff(np.linspace(0, width, len(lines) + 1).round()).astype(int)   # adding up to `width`

    figs, titles, windows = [], [], []
    for fig_width, (name, w3, w4) in zip(fig_widths, lines[['w3', 'w4']].itertuples(name=None)):

        window = (x >= w3) & (x <= w4)
        title = f'{name}, {plot_vars['latex']['CHI2_LABEL']} = {_interval_chi2(grid, best_model, w3, w4):.2f}'
        ylim = _zoom_ylim(y[window], best_model[window])

        fig_kw = dict(title=title, x_range=Range1d(w3, w4),
                      **_box(fig_width, height, plot_vars['bokeh']['BORDER_INNER'], plot_vars['bokeh']['BORDER_RIGHT'],
                                                plot_vars['bokeh']['BORDER_TOP'], plot_vars['bokeh']['BORDER_BOTTOM']))
        if ylim is not None:
            fig_kw['y_range'] = Range1d(*ylim)
        fig = _bk_figure(style, **fig_kw)

        fig.step(x, y, mode='center', color=c_fg, line_width=1)
        fig.step(x, best_model, mode='center', color=c_fit, line_width=1)
        if sel_src is not None:
            fig.step('x', 'model', source=sel_src, mode='center', color=pick_color, line_width=1, linestyle='--')

        fig.yaxis.visible = False
        fig.xaxis[0].ticker.desired_num_ticks = 3
        fig.xaxis[0].formatter.use_scientific = False
        fig.title.text_font_size = '9pt'

        figs.append(fig)
        titles.append(title)
        windows.append([float(w3), float(w4)])

    return figs, titles, windows


# JavaScript of the draw picker (runs in the browser, so it also works in the saved HTML). Same steps as
# _connect_draw_picker: nearest draw in screen pixels, SSPGrid.model of its node, reduced chi2 as _reduced_chi2 and
# _interval_chi2, conditions of the node the draw maps to.
_PICK_JS = """
const X = src.data[xcol], Y = src.data[ycol]
const kx = (fig.inner_width || frame) / Math.abs(fig.x_range.end - fig.x_range.start)
const ky = (fig.inner_height || frame) / Math.abs(fig.y_range.end - fig.y_range.start)
let k = -1, d2min = Infinity
for (let i = 0; i < X.length; i++) {
    const dx = (X[i] - cb_obj.x) * kx, dy = (Y[i] - cb_obj.y) * ky
    const d2 = dx * dx + dy * dy
    if (d2 < d2min) { d2min = d2; k = i }
}
if ((k < 0) || (d2min > radius * radius)) { return }

const row = src.data.row[k]
if (row < 0) {
    console.warn(`Draw ${k} falls on an (age, Z) combination missing from the grid`)
    return
}

// Model of the draw: 10^log(A) * node flux * transmission^E(B-V)
const d = spec.data, n = d.wave.length, flux = d['flux_' + row]
const ebv = src.data[names[2]][k], log_amp = src.data[names[3]][k], amp = Math.pow(10, log_amp)
const model = new Float64Array(n), resid = new Float64Array(n)
const win = zwin.map(() => [0, 0])
let chi2 = 0, n_good = 0
for (let i = 0; i < n; i++) {
    model[i] = amp * flux[i] * Math.pow(d.red_corr[i], ebv)
    const res = (d.lum[i] - model[i]) / d.lum[i]
    resid[i] = Number.isFinite(res) ? res : NaN
    if (d.good[i]) {
        const term = ((d.lum[i] - model[i]) / d.lum_err[i]) ** 2
        chi2 += term
        n_good += 1
        zwin.forEach((w, j) => {
            if ((d.wave[i] >= w[0]) && (d.wave[i] <= w[1])) { win[j][0] += term; win[j][1] += 1 }
        })
    }
}
sel.data = {x: d.wave, model: model, resid: resid}

// Draw location: open circle on the scatter panels (shared selection) and line on the histograms
src.selected.indices = [k]
spans.forEach((span, p) => { span.location = src.data[names[p]][k]; span.visible = true })

// Conditions and reduced chi2 of the draw, in parentheses after those of the fit
const fmt = (v, nd) => Number.isFinite(v) ? v.toFixed(nd) : 'nan'
const chi2_red = (n_good > n_par) ? chi2 / (n_good - n_par) : NaN
const vals = [node_age[row], Math.pow(10, node_z[row]) / z_norm, ebv, log_amp, chi2_red]
const html = fit_parts.map((part, p) => `${part} (<span style="color:${color}">${fmt(vals[p], dec[p])}</span>)`)
div.text = div_tmpl.replace('{}', () => html.join(sep))
zfigs.forEach((f, j) => { f.title.text = `${ztitles[j]} (${fmt(win[j][0] / win[j][1], 2)})` })
"""

_CLEAR_JS = """
src.selected.indices = []
const nan = new Float64Array(sel.data.x.length).fill(NaN)
sel.data = {x: sel.data.x, model: nan, resid: nan.slice()}
spans.forEach((span) => { span.visible = false })
div.text = div_tmpl.replace('{}', () => fit_parts.join(sep))
zfigs.forEach((f, j) => { f.title.text = ztitles[j] })
"""


def _connect_draw_picker_bokeh(mat, src, grid, sel_src, markers, div, div_tmpl, sep, fit_parts, decimals, z_norm,
                               zoom_figs=(), zoom_titles=(), zoom_windows=(), names=None, color='magenta', frame=120,
                               pick_radius=10):

    """
    Bokeh version of _connect_draw_picker. A click on an off-diagonal panel of `mat` selects the nearest draw of `src`
    (within `pick_radius` pixels): its model and residuals fill `sel_src`, the draw is marked on every panel and its
    conditions and reduced chi2 are added to `div` and to the line-zoom titles. A double click on any panel clears
    the selection. The grid nodes and the observation travel to the browser in a ColumnDataSource.
    """

    names = list(src.data)[:4] if names is None else list(names)
    good = (grid.mask & np.isfinite(grid.lum) & (grid.lum_err > 0)).astype(np.int32)
    spec_src = ColumnDataSource({'wave': grid.wave, 'lum': grid.lum, 'lum_err': grid.lum_err,
                                 'red_corr': grid.red_corr, 'good': good,
                                 **{f'flux_{k}': flux for k, flux in enumerate(grid.grid_flux)}})

    args = dict(src=src, spec=spec_src, sel=sel_src, spans=list(markers), div=div, div_tmpl=div_tmpl, sep=sep,
                fit_parts=list(fit_parts), dec=list(decimals), z_norm=float(z_norm),
                node_age=grid.age_node.tolist(), node_z=grid.z_node.tolist(),
                zfigs=list(zoom_figs), ztitles=list(zoom_titles), zwin=[list(w) for w in zoom_windows],
                names=names, color=color, n_par=len(names), radius=pick_radius, frame=frame)

    clear = CustomJS(args=args, code=_CLEAR_JS)
    n = len(mat)
    for i in range(n):
        for j in range(i + 1):
            mat[i][j].js_on_event('doubletap', clear)
            if i != j:
                pick = CustomJS(args=dict(args, fig=mat[i][j], xcol=names[j], ycol=names[i]), code=_PICK_JS)
                mat[i][j].js_on_event('tap', pick)

    return


def ssp_fitting_infograph_bokeh(trace, fname=None, truth_dict=None, title=None, show_stellar=True, fig_cfg=None,
                                pick_draws=True, residual_df=None, solar_Z=False, z_sun=Z_SUN, show_grid=False,
                                width=850, cell=120, max_samples=20000, max_trace_draws=500, pick_color='magenta',
                                return_fig=False, **kwargs):

    """
    Bokeh version of ssp_fitting_infograph: spectrum and fit, residuals (or line zooms), scatter matrix and traces.
    The arguments shared with ssp_fitting_infograph behave as there. With `pick_draws`, a click on a draw in an
    off-diagonal scatter panel overplots its model and a double click clears it; this runs in JavaScript, so it also
    works in the saved HTML file.

    Every panel has a fixed outer size: the spectrum is as tall as the scatter-matrix rows and shares the top and
    bottom borders of the first and last ones (BORDER_* constants), so their frames line up.

    width           : outer width of the spectrum panel in pixels
    cell            : side of each scatter-matrix frame in pixels
    kwargs          : log_scale=True for a logarithmic flux axis
    max_samples     : posterior draws shown (and pickable) in the scatter panels; the histograms use all of them
    max_trace_draws : draws per chain shown in the trace panels
    pick_color      : color of the selected draw
    return_fig      : return the layout without displaying or saving it
    """

    if not bokeh_check:
        raise SpecSyError('Please install Bokeh to use ssp_fitting_infograph_bokeh')

    d = _ssp_trace_data(trace, truth_dict, solar_Z, z_sun, show_grid)
    names, grid, x, y, mask = d.names, d.grid, d.x, d.y, d.mask
    style = _bk_style(fig_cfg)
    c_fit, c_band = _bk_color('fitted_SSP'), _bk_color('predictive_band')

    # Observed spectrum as a LiMe object (stored mask: True = fitted, LiMe pixel_mask: True = masked)
    spec = lime.Spectrum(x, y, d.inputs['lum_err'].values, redshift=0, pixel_mask=~mask,
                         units_wave=d.inputs.attrs['units_wave'], units_flux=d.inputs.attrs['units_flux'])

    # Line list: up to 4 lines inside the spectrum replace the residuals panel by a zoom on their [w3, w4] intervals
    if residual_df is not None:
        residual_df = spec.retrieve.lines_frame(ref_bands=residual_df, adjust_central_band=False)
        if (residual_df is None) or (len(residual_df) == 0):
            _logger.warning('None of the residual_df lines is inside the spectrum range: the residuals panel is '
                            'shown instead.')
            residual_df = None
        else:
            residual_df = residual_df.iloc[:4]

    # Model and residuals of the selected draw, filled by the picker
    sel_src = ColumnDataSource(dict(x=x, model=np.full(x.size, np.nan), resid=np.full(x.size, np.nan)))

    # Layout: top row = spectrum | scatter matrix, with the spectrum as tall as the matrix rows; bottom row =
    # conditions + residuals (or line zooms) | "Traces" + trace panels, with the same heights on both sides
    n = len(names)
    top_height = sum(_matrix_sizes(n, cell)[1])
    low_height = plot_vars['bokeh']['BORDER_TOP'] + int(1.2 * cell) + plot_vars['bokeh']['BORDER_BOTTOM']
    trace_size = int(0.8 * cell)

    # Spectrum panel, with the top and bottom borders of the first and last matrix rows so the frames line up. The
    # input spectrum is drawn here with LiMe's axis labels: LiMe's spectrum() applies its own figure size and layout
    fig0 = _bk_figure(style, toolbar='right', title=title, y_axis_type='log' if kwargs.get('log_scale') else 'linear',
                      **_box(width, top_height, 80, 10, plot_vars['bokeh']['BORDER_TOP'], plot_vars['bokeh']['BORDER_BOTTOM']))
    ax_cfg = lime_theme.ax_defaults(None, spec, plotting_library='bokeh')
    fig0.xaxis.axis_label, fig0.yaxis.axis_label = ax_cfg['xlabel'], ax_cfg['ylabel']
    fig0.step(x, y, mode='center', color=_bk_color('fg'), line_width=1, legend_label='Input spectrum')

    # Predictive band, stellar / nebular components and best fit (Bokeh draws in creation order)
    if d.band is not None:
        rend = fig0.varea_step(x=x, y1=d.band[0], y2=d.band[2], step_mode='center', fill_color=c_band,
                               fill_alpha=0.6, legend_label='Predictive 16-84%')
        rend.level = 'underlay'

    star, neb = grid.components(d.best_theta)
    if neb is not None:
        if show_stellar:
            fig0.step(x, star, mode='center', color='darkorange', line_width=0.8, line_dash='dotted',
                      legend_label='Stellar')
        fig0.step(x, neb, mode='center', color='crimson', line_width=0.8, line_dash='dotted',
                  legend_label='Nebular continuum')

    chi2_red = _reduced_chi2(grid, d.best_model, len(names))
    fig0.step(x, d.best_model, mode='center', color=c_fit, line_width=1,
              legend_label=f'Fitted spectrum ({plot_vars['latex']["CHI2_LABEL"]} = {chi2_red:.2f})')

    if pick_draws:
        fig0.step('x', 'model', source=sel_src, mode='center', color=pick_color, line_width=1)

    fig0.legend.location = 'top_right'
    fig0.legend.click_policy = 'hide'

    # Lower-left slot: residuals panel (same borders and toolbar as the spectrum, so the frames line up) or one zoom
    # figure per line
    if residual_df is None:
        fig1 = _bk_figure(style, toolbar='right', x_range=fig0.x_range, y_range=Range1d(-0.7, 0.7),
                          **_box(width, low_height, 80, 10, plot_vars['bokeh']['BORDER_TOP'], plot_vars['bokeh']['BORDER_BOTTOM']))
        residuals_plot_bokeh(fig1, x, y, d.best_model, band=d.band)
        if pick_draws:
            fig1.step('x', 'resid', source=sel_src, mode='center', color=pick_color, line_width=1)
        low, spec_figs = fig1, [fig0, fig1]
        zoom_figs, zoom_titles, zoom_windows = [], [], []
    else:
        zoom_figs, zoom_titles, zoom_windows = line_zoom_bokeh(x, y, d.best_model, residual_df, grid, style, width,
                                                               low_height, sel_src if pick_draws else None,
                                                               pick_color)
        low, spec_figs = _fixed(row, *zoom_figs, width=width, height=low_height), [fig0, *zoom_figs]

    _bk_mask_shading(spec_figs, x, mask)

    # Physical conditions of the fit between the spectrum and the lower panels
    fit_parts, decimals, z_norm = _condition_parts(d.best_theta, chi2_red, d.z_sun)
    div_tmpl = f'<div style="width:{width}px; text-align:center; font-size:10pt">{{}}</div>'
    sep = ',&nbsp;&nbsp;'
    cond_div = Div(text=div_tmpl.replace('{}', sep.join(fit_parts)), width=width, height=plot_vars['bokeh']['DIV_HEIGHT'],
                   margin=(0, 0, 0, 0))

    # Draws shown in the scatter panels (plotted coordinates) with the grid row of their node
    theta = np.column_stack([d.raw_chains[k].ravel() for k in names])
    idcs = np.arange(len(theta))
    if len(theta) > max_samples:
        idcs = np.sort(np.random.default_rng(0).choice(len(theta), size=max_samples, replace=False))
        _logger.info(f'The scatter panels show {max_samples} of the {len(theta)} posterior draws (max_samples).')
    src = ColumnDataSource({**{k: d.samples[k][idcs] for k in names},
                            'row': grid.rows(theta[idcs, 0], theta[idcs, 1]).astype(np.int32)})

    # Scatter matrix and traces
    mat, markers = scatter_matrix_bokeh(src, d.samples, names, d.pct, style, p0=d.p0, p0_arr=d.p0_arr, ref=d.ref,
                                        labels=d.labels, solar=d.z_sun is not None, nodes=d.nodes, cell=cell,
                                        pick_color=pick_color)
    trace_figs = trace_plots_bokeh(d.chains, names, d.pct, style, p0=d.p0, ref=d.ref, labels=d.labels,
                                   size=trace_size, height=low_height, max_draws=max_trace_draws)
    traces_width = sum(fig.width for fig in trace_figs)
    traces_div = Div(text='<div style="text-align:center; font-size:10pt">Traces</div>', width=traces_width,
                     height=plot_vars['bokeh']['DIV_HEIGHT'], margin=(0, 0, 0, 0))

    if pick_draws:
        _connect_draw_picker_bokeh(mat, src, grid, sel_src, markers, cond_div, div_tmpl, sep, fit_parts, decimals,
                                   z_norm, zoom_figs, zoom_titles, zoom_windows, names=names, color=pick_color,
                                   frame=cell)

    # Fixed-size rows and columns all the way up, so Bokeh keeps every box at the size computed here
    matrix, matrix_width, _ = _matrix_layout(mat, cell)
    low_total = plot_vars['bokeh']['DIV_HEIGHT'] + low_height
    top_row = _fixed(row, fig0, _spacer(plot_vars['bokeh']['GAP'], top_height), matrix,
                     width=width + plot_vars['bokeh']['GAP'] + matrix_width, height=top_height)
    bottom_row = _fixed(row,
                        _fixed(column, cond_div, low, width=width, height=low_total),
                        _spacer(plot_vars['bokeh']['GAP'], low_total),
                        _fixed(column, traces_div, _fixed(row, *trace_figs, width=traces_width, height=low_height),
                               width=traces_width, height=low_total),
                        width=width + plot_vars['bokeh']['GAP'] + traces_width, height=low_total)
    layout = _fixed(column, top_row, bottom_row,
                    width=width + plot_vars['bokeh']['GAP'] + max(matrix_width, traces_width),
                    height=top_height + low_total)

    # Save or display the plot
    if return_fig:
        return layout

    if fname is not None:
        save(layout, filename=fname, resources=CDN, title=title or 'SSP fit')
    else:
        show(layout)

    return layout