import logging
import numpy as np
import pyneb as pn

from pandas import DataFrame
from lime import label_decomposition, load_lines_log
from pathlib import Path

from uncertainties import unumpy, ufloat
from lmfit.models import LinearModel

from ..tools import blended_label_from_log, get_mixed_fluxes
from ..plots import extinction_gradient


_logger = logging.getLogger('SpecSy')


# Function to compute and plot cHbeta
def cHbeta_from_log(log, line_list='all', R_V=3.1, law='G03 LMC', temp=10000.0, den=100.0, ref_line='auto',
                    flux_entry='gauss', lines_ignore=None, show_plot=False, plot_address=None, plot_title=r'$c(H\beta)$ calculation',
                    fig_cfg={}, ax_cfg={}):

    '''

    This function computes the logarithmic extinction coefficient using the hydrogen lines on the input logs.

    The user can provide a list with the lines to use in the coefficient calculation.

    Moreover, the user can also provide a list of lines to exclude in the coefficient calculation.

    The user can provide the normalization line. If none is provided, the function will try to use Hbeta (H1_4861A). If
    H1_4861A is not in the input log, the library will use the second most intense hydrogen line for the normalization.

    The user can select the flux type ("intg" or "gauss") for the calculation. The default type is "gauss".

    The function also returns the coefficient uncertainty. This value is close to zero if there are only two Hydrogen lines.
    If there aren't hydrogen lines in the log or there are conflicts in the calculation, the function returns "None" for both variables.

    The user can also request the plot with the coefficient calculation. If a file address is provided this plot will be
    stored at the location. In this plot, the "lines_ignore" will be included in the plot even though they are not used
    in the coefficient calculation.

    Logs with hydrogen lines with multiple kinematic components can cause issues in the calculation. The user should index
    the input dataframe lines log to make sure it only includes hydrogen lines from the same kinematic component.

    The emissivities are calculated with the input temperature and density using PyNeb.

    :param log: Lines log with the input fluxes. The pandas dataframe must adhere to LiMe formatting
    :type log: pd.DataFrame

    :param line_list: Array with the lines to use for the cHbeta calculation. If none provided, all lines will be used.
    :type line_list: list, optional

    :param R_V: Total-to-selective extinction ratio. The default value is 3.1
    :type R_V: float, optional

    :param law: Extinction law. The default value is "G03 LMC" from the Gordon et al. (2003, ApJ, 594, 279). The reddening law name should follow the pyneb notation.
    :type law: str, optional

    :param temp: Temperature for the emissivity calculation in degrees Kelvin. The default value is 10000 K.
    :type temp: float, optional

    :param den: Density for the emissivity calculation in particles per centimeter cube. The default value is 100 cm^-3.
    :type den: float, optional

    :param ref_line: Line label of the normalization flux. The default value is "auto" for the automatic selection.
    :type ref_line: str, optional

    :param flux_entry: Flux type for the cHbeta calculation. The default value is "gauss" for a Gaussian flux selection.
    :type flux_entry: str, optional

    :param lines_ignore: List of lines to exclude in the cHbeta calculation. The default value is None.
    :type lines_ignore: list, optional

    :param show_plot: Check to display the cHbeta calculation regression. The default value is False.
    :type show_plot: bool, optional

    :param plot_address: Address for the output image with the cHbeta calculation regression. The default value is None.
    :type plot_address: str, optional

    :param plot_title: Title for the cHbeta calculation regression plot.
    :type plot_title: str, optional

    :param fig_cfg: Configuration for the cHbeta plot figure.
    :type fig_cfg: dict, optional

    :param ax_cfg: Configuration for the cHbeta plot axes.
    :type ax_cfg: dict, optional

    :return: cHbeta value and uncertainty.
    :rtype: float, float

    '''

    # Check if input file is a log or an address to the log file
    if not isinstance(log, DataFrame):
        log_path = Path(log)
        if log_path.is_file():
            log = load_lines_log(log_path)
        else:
            _logger.warning(f'- The file {log} could not be found')
            raise TypeError()

    # Use all hydrogen lines if a list is not specified
    if line_list == 'all':

        # Check for the ion column:
        if 'ion' in log.columns:
            idcs_lines = log.ion == 'H1'
        else:
            ion_array, wave_array, latex_array = label_decomposition(log.index.values)
            idcs_lines = ion_array == 'H1'

        line_list = log.loc[idcs_lines].index.values

    # Proceed if there are enough lines to compute the extinction
    if line_list.size > 1:

        # Use the second most intense line to normalize if non provided
        if ref_line == 'auto':

            # First try to use Hbeta.
            if 'H1_4861A' in log.index:
                ref_line = 'H1_4861A'

            else:
                ref_excluded = lines_ignore if lines_ignore is not None else []
                idcs_candidates = log.index.isin(line_list) & ~log.index.isin(ref_excluded)

                He_cand, fluxes_cand = log.loc[idcs_candidates].index.values, log.loc[idcs_candidates].gauss_flux.values
                ref_line = He_cand[np.argsort(fluxes_cand)[-2]]

        # Check the reference line is there
        if ref_line in log.index:

            # Label the lines which are found in the lines log
            ion_ref, waves_ref, latexLabels_ref = label_decomposition(ref_line, scalar_output=True)
            ion_array, waves_array, latex_array = label_decomposition(line_list)

            # Get the latex labels from the dataframe
            latexLabels_ref = log.loc[ref_line].latex_label
            latex_array = log.loc[idcs_lines].latex_label.values

            # Mixed fluxes ratios
            if flux_entry == 'auto':

                # Integrated fluxes for single lines and gaussian for blended
                obsFlux, obsErr = get_mixed_fluxes(log)
                obsFlux, obsErr = obsFlux[idcs_lines], obsErr[idcs_lines]

                # Check if reference line is blended
                if (log.loc[ref_line, 'profile_label'] == 'no') | (ref_line.endswith('_m')):
                    ref_flux_type = 'intg'
                else:
                    ref_flux_type = 'gauss'

                # Same for the reference line
                Href_flux = log.loc[ref_line, f'{ref_flux_type}_flux']
                Href_err = log.loc[ref_line, f'{ref_flux_type}_err']

            # Use the user param
            else:
                obsFlux = log.loc[idcs_lines, f'{flux_entry}_flux'].values
                obsErr = log.loc[idcs_lines, f'{flux_entry}_flux'].values

                Href_flux = log.loc[ref_line, f'{flux_entry}_flux']
                Href_err = log.loc[ref_line, f'{flux_entry}_flux']

            # Check for negative or nan entries in Href
            if not np.isnan(Href_flux) and not (Href_flux < 0):

                idcs_flux_invalid = np.isnan(obsFlux) | (obsFlux < 0)
                idcs_err_invalid = np.isnan(obsErr) | (obsErr < 0)

                # Check for negative or nan entries in Href to remove them
                if np.any(idcs_flux_invalid):
                    _logger.warning(f'Lines with bad flux entries: {line_list[idcs_flux_invalid]} ='
                                    f' {obsFlux[idcs_flux_invalid]}')

                if np.any(idcs_err_invalid):
                    _logger.warning(f'Lines with bad error entries: {line_list[idcs_err_invalid]} ='
                                    f' {obsErr[idcs_err_invalid]}')

                idcs_valid = ~idcs_flux_invalid & ~idcs_err_invalid
                line_list = line_list[idcs_valid]
                obsFlux = obsFlux[idcs_valid]
                obsErr = obsErr[idcs_valid]
                waves_array = waves_array[idcs_valid]
                latex_array = latex_array[idcs_valid]

                if line_list.size > 1:

                    # Check if there are repeated entries
                    unique_array, counts = np.unique(waves_array, return_counts=True)
                    if np.any(counts > 1):
                        _logger.warning(f'These lines wavelengths are repeated: {unique_array[counts > 1]}\n'
                                        f'Check for repeated transitions or multiple kinematic components.\n')

                    # Array to compute the uncertainty # TODO need own method to propagate the uncertainty
                    obsRatio_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(Href_flux, Href_err)

                    # Theoretical ratios
                    H1 = pn.RecAtom('H', 1)
                    refEmis = H1.getEmissivity(tem=temp, den=den, wave=waves_ref)
                    emisIterable = (H1.getEmissivity(tem=temp, den=den, wave=wave) for wave in waves_array)
                    linesEmis = np.fromiter(emisIterable, float)
                    theoRatios = linesEmis / refEmis

                    # Reddening law
                    rc = pn.RedCorr(R_V=R_V, law=law)
                    Xx_ref, Xx = rc.X(waves_ref), rc.X(waves_array)
                    f_lines = Xx/Xx_ref - 1
                    f_ref = Xx_ref/Xx_ref - 1

                    # cHbeta linear fit values
                    x_values = f_lines - f_ref
                    y_values = np.log10(theoRatios) - unumpy.log10(obsRatio_uarray)

                    # Exclude from the linear fitting the lines requested by the user
                    if lines_ignore is not None:
                        idcs_valid = ~np.in1d(line_list, lines_ignore)
                    else:
                        idcs_valid = np.ones(line_list.size).astype(bool)

                    # Perform fit
                    lineModel = LinearModel()
                    y_nom, y_std = unumpy.nominal_values(y_values), unumpy.std_devs(y_values)

                    pars = lineModel.make_params(intercept=y_nom[idcs_valid].min(), slope=0)
                    output = lineModel.fit(y_nom[idcs_valid], pars, x=x_values[idcs_valid], weights=1/y_std[idcs_valid])

                    cHbeta, cHbeta_err = output.params['slope'].value, output.params['slope'].stderr
                    intercept, intercept_err = output.params['intercept'].value, output.params['intercept'].stderr

                    if x_values[idcs_valid].size == 2:
                        cHbeta_err, intercept_err = 0.0, 0.0

                    # Case lmfit cannot fit the error bars, switch none by nan
                    if not output.errorbars:
                        cHbeta_err, intercept_err = np.nan, np.nan

                    if show_plot:
                        extinction_gradient((cHbeta, cHbeta_err), (intercept, intercept_err),
                                            x_values, (y_nom, y_std),
                                            line_labels=latex_array, ref_label=latexLabels_ref,
                                            idcs_valid=idcs_valid,
                                            save_address=plot_address, title=plot_title,
                                            fig_cfg=fig_cfg, ax_cfg=ax_cfg)

                else:
                    _logger.info(f'{"Zero H1 lines" if line_list.size == 0 else "Just one H1 line"} in the input log, '
                                 f' extinction coefficient could not be calculated')

            else:
                _logger.warning(f'Reference line {ref_line} had an invalid flux value of {Href_flux}')
                cHbeta, cHbeta_err = None, None

        else:
            _logger.info(f'The normalization line {ref_line} could not be found in input log')
            raise IndexError()

    else:
        _logger.info(f'{"Zero H1 lines" if line_list.size == 0 else "Just one H1 line"} in the input log, extinction coefficient '
                     f'could not be calculated')
        cHbeta, cHbeta_err = None, None

    return cHbeta, cHbeta_err
