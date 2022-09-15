import numpy as np
import astropy
from lime.plots import save_close_fig_swicth, STANDARD_PLOT

from matplotlib import pyplot as plt, gridspec, patches, rc_context, cm, colors


def extinction_gradient(cHbeta_array, n_array, x, y_array, idcs_valid=None, line_labels=None, ref_label='ref',
                        save_address=None,  title=None, fig_cfg={}, ax_cfg={}):

    # Adjust default theme
    PLOT_CONF = STANDARD_PLOT.copy()

    # Adjust the axis labels to include the reference line
    x_label = r'$f_{\lambda} - $' + f'$f_{{{ref_label.replace("$","")}}}$'
    y_label = r'$\left(\frac{I_{\lambda}}{I_{ref}}\right)_{Theo}-\left(\frac{F_{\lambda}}{I_{ref}}\right)_{Obs}$'
    y_label = y_label.replace('ref', ref_label.replace("$", ""))
    AXES_CONF = {'xlabel': x_label, 'ylabel': y_label, 'title': title}

    # User configuration overrites user
    PLT_CONF = {**PLOT_CONF, **fig_cfg}
    AXES_CONF = {**AXES_CONF, **ax_cfg}

    # Draw the figure
    with rc_context(PLT_CONF):

        cHbeta, cHbeta_err = cHbeta_array
        n, n_err = n_array
        y, y_err = y_array

        fig, ax = plt.subplots()
        ax.set(**AXES_CONF)

        # Plot valid entries
        idcs_valid = np.ones(x.size).astype(bool) if idcs_valid is None else idcs_valid
        valid_scatter = ax.errorbar(x[idcs_valid], y[idcs_valid], y_err[idcs_valid], fmt='o')

        # Plot excluded entries
        if np.any(~idcs_valid):
            ax.errorbar(x[~idcs_valid], y[~idcs_valid], y_err[~idcs_valid], fmt='o',
                        color='tab:red', label='excluded lines')

        # Linear fitting
        linear_fit = cHbeta * x + n
        linear_label = r'$c(H\beta)={:.2f}\,\pm\,{:.2f}$'.format(cHbeta, cHbeta_err)
        ax.plot(x, linear_fit, linestyle='--', label=linear_label)

        # Labels for the lines
        for i, lineWave in enumerate(line_labels):
            ax.annotate(lineWave,
                        xy=(x[i], y[i]),
                        xytext=(x[i], y[i] + 1.25 * y_err[i]),
                        horizontalalignment="center",
                        rotation=90,
                        xycoords='data', textcoords=("data", "data"))

        # Legend
        ax.legend(loc=3, ncol=2)

        # Increase upper limit
        y_lims = ax.get_ylim()
        ax.set_ylim(y_lims[0], y_lims[1] * 2)

        # Display/save the figure
        save_close_fig_swicth(save_address, 'tight', fig_obj=fig)

    return