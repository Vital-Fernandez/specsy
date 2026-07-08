import numpy as np
import arviz as az

from pathlib import Path
from matplotlib import pyplot as plt, gridspec, rc_context, cm, colors

import lime
from innate import load_dataset
from lime import label_decomposition
from lime.plotting.plots import save_close_fig_swicth
from lime.plotting.format import Themer
from specsy import _setup_cfg
from specsy.tools import linear_regression
from specsy.io import SpecSyError, specsy_cfg

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


# Formatting object
theme_file = Path(__file__).resolve().parent/'specsy_theme.toml'
theme = Themer.from_toml(theme_file)


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


def az_trace(sampling_result, fname=None, in_fig=None, var_list=None, exclude=['theo_flux'], var_latex=specsy_cfg['latex_param_notation'],
             fig_cfg=None, maximize=False, display_check=True):

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


def az_scatter_matrix(sampling_result, fname=None, in_fig=None, var_list=None, exclude=['theo_flux'], var_latex=specsy_cfg['latex_param_notation'],
             fig_cfg=None, maximize=False, display_check=True, n_cols=5):

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

    # Set the plot format where the user's overwrites the default
    # size_conf = {'figure.figsize': (6, len(var_list) * 2),
    #              'axes.prop_cycle': plt.cycler(color=[theme.colors['fg']]),
    #              'axes.titlesize': 16, 'axes.titlepad': 10}
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

    # # Set the plot format where the user's overwrites the default
    # size_conf = {'figure.figsize': (6, len(var_list) * 2),
    #              'axes.prop_cycle': plt.cycler(color=[theme.colors['fg']]),
    #              'axes.titlesize': 16, 'axes.titlepad': 10}
    # size_conf = size_conf if fig_cfg is None else {**size_conf, **fig_cfg}

