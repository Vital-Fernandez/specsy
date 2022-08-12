import numpy as np
import pyneb as pn
import lime

from matplotlib import pyplot as plt, rcParams
from lime.plots import STANDARD_PLOT
from uncertainties import unumpy, ufloat
from lmfit.models import LinearModel


# Function to compute and plot cHbeta
def cHbeta_from_log(line_df, line_labels='all', temp=10000.0, den=100.0, ref_wave='H1_4861A',
                    comp_mode='auto', plot_address=False, input_flux_labels=('line_flux', 'line_err'), title=''):

    # Use all hydrogen lines if none are defined
    if np.isscalar(line_labels):
        if line_labels == 'all':
            idcs_H1 = line_df.ion == 'H1'
            line_labels = line_df.loc[idcs_H1].index.values

    if len(line_df.index) > 1 and np.sum(line_df.ion == 'H1') > 1:

        # Loop through the input lines
        assert ref_wave in line_df.index, f'- ERROR: {ref_wave} not found in input lines log dataframe for c(Hbeta) calculation'

        # Label the lines which are found in the lines log
        # TODO this function fails if we have negative fluxes leaving only one line
        idcs_lines = line_df.index.isin(line_labels) & (line_df.intg_flux > 0) & (line_df.gauss_flux > 0)
        line_labels = line_df.loc[idcs_lines].index.values
        ion_ref, waves_ref, latexLabels_ref = lime.label_decomposition(ref_wave, scalar_output=True)
        ion_array, waves_array, latexLabels_array = lime.label_decomposition(line_labels)


        # Observed ratios
        if comp_mode == 'auto':
            if ('line_flux' in line_df.columns) and ('line_err' in line_df.columns):
                obsFlux = line_df.loc[idcs_lines, input_flux_labels[0]].values
                obsErr = line_df.loc[idcs_lines, input_flux_labels[1]].values
                Href_flux, Href_err = line_df.loc[ref_wave, 'line_flux'], line_df.loc[ref_wave, 'line_err']

            else:
                Href_flux, Href_err = line_df.loc[ref_wave, 'intg_flux'], line_df.loc[ref_wave, 'intg_err']  # TODO what if Hbeta is blended
                obsFlux, obsErr = np.empty(line_labels.size), np.empty(line_labels.size)
                slice_df = line_df.loc[idcs_lines]
                idcs_intg = slice_df.blended_label == 'None'
                obsFlux[idcs_intg] = slice_df.loc[idcs_intg, 'intg_flux'].values
                obsErr[idcs_intg] = slice_df.loc[idcs_intg, 'intg_err'].values
                obsFlux[~idcs_intg] = slice_df.loc[~idcs_intg, 'gauss_flux'].values
                obsErr[~idcs_intg] = slice_df.loc[~idcs_intg, 'gauss_err'].values

            obsRatio_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(Href_flux,
                                                                      Href_err)  # TODO unumpy this with your own model

        elif comp_mode == 'gauss':
            Href_flux, Href_err = line_df.loc[ref_wave, 'gauss_flux'], line_df.loc[ref_wave, 'gauss_err']
            obsFlux, obsErr = line_df.loc[idcs_lines, 'gauss_flux'], line_df.loc[idcs_lines, 'gauss_err']
            obsRatio_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(Href_flux, Href_err)

        else:
            Href_flux, Href_err = line_df.loc[ref_wave, 'intg_flux'], line_df.loc[ref_wave, 'intg_err']
            obsFlux, obsErr = line_df.loc[idcs_lines, 'intg_flux'], line_df.loc[idcs_lines, 'intg_err']
            obsRatio_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(Href_flux, Href_err)

        assert not np.any(np.isnan(obsFlux)) in obsFlux, '- ERROR: nan entry in input fluxes for c(Hbeta) calculation'
        assert not np.any(
            np.isnan(obsErr)) in obsErr, '- ERROR: nan entry in input uncertainties for c(Hbeta) calculation'

        # Theoretical ratios
        H1 = pn.RecAtom('H', 1)
        refEmis = H1.getEmissivity(tem=temp, den=den, wave=waves_ref)
        emisIterable = (H1.getEmissivity(tem=temp, den=den, wave=wave) for wave in waves_array)
        linesEmis = np.fromiter(emisIterable, float)
        theoRatios = linesEmis / refEmis

        # Reddening law
        rc = pn.RedCorr(R_V=3.1, law='CCM89')
        Xx_ref, Xx = rc.X(waves_ref), rc.X(waves_array)
        f_lines = Xx / Xx_ref - 1
        f_ref = Xx_ref / Xx_ref - 1

        # cHbeta linear fit values
        x_values = f_lines - f_ref
        y_values = np.log10(theoRatios) - unumpy.log10(obsRatio_uarray)

        # Perform fit
        lineModel = LinearModel()
        y_nom, y_std = unumpy.nominal_values(y_values), unumpy.std_devs(y_values)
        pars = lineModel.make_params(intercept=y_nom.min(), slope=0)
        output = lineModel.fit(y_nom, pars, x=x_values, weights=1 / y_std)
        cHbeta, cHbeta_err = output.params['slope'].value, output.params['slope'].stderr
        intercept, intercept_err = output.params['intercept'].value, output.params['intercept'].stderr

        if plot_address:

            STANDARD_PLOT = {'figure.figsize': (14, 7), 'axes.titlesize': 12, 'axes.labelsize': 14,
                             'legend.fontsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10}

            axes_dict = {'xlabel': r'$f_{\lambda} - f_{H\beta}$',
                         'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
                         'title':  title + f'\n'+ r'$c(H\beta)$ calculation with ' + ", ".join(list(latexLabels_array))}

            rcParams.update(STANDARD_PLOT)

            fig, ax = plt.subplots(figsize=(8, 4))
            fig.subplots_adjust(bottom=-0.7)

            # Data ratios
            err_points = ax.errorbar(x_values, y_nom, y_std, fmt='o')

            # Linear fitting
            linear_fit = cHbeta * x_values + intercept
            linear_label = r'$c(H\beta)={:.2f}\pm{:.2f}$'.format(cHbeta, cHbeta_err)
            ax.plot(x_values, linear_fit, linestyle='--', label=linear_label)
            ax.update(axes_dict)

            # Legend
            ax.legend(loc='best')
            ax.set_ylim(-0.5, 0.5)

            # Generate plot
            plt.tight_layout()
            # if isinstance(plot_address, (str, pathlib.WindowsPath, pathlib.PosixPath)):
            #     # crs = mplcursors.cursor(_ax, hover=True)
            #     # crs.connect("add", lambda sel: sel.annotation.set_text(sel.annotation))
            #     plt.savefig(plot_address, dpi=200, bbox_inches='tight')
            # else:
            #     mplcursors.cursor(_ax).connect("add",
            #                                   lambda sel: sel.annotation.set_text(latexLabels_array[sel.target.index]))
            # plt.show()
            plt.savefig(plot_address)

    else:

        cHbeta, cHbeta_err = '-', '-'
        print('- Only one hydrogen line, cHbeta could not be calculated')

    return cHbeta, cHbeta_err
