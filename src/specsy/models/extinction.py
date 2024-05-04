import logging
import numpy as np
try:
    import pyneb as pn
    pyneb_check = True
except ImportError:
    pyneb_check = False

from pandas import DataFrame
from lime import label_decomposition
from pathlib import Path
from lime.io import load_frame
from uncertainties import unumpy, ufloat
from lmfit.models import LinearModel

from ..tools import get_mixed_fluxes
from ..io import SpecSyError, check_file_dataframe

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
            log = load_log(log_path)
        else:
            _logger.warning(f'- The file {log} could not be found')
            raise TypeError()

    # Use all hydrogen lines if a list is not specified
    if isinstance(line_list, str):
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

                enougth_lines = True
                if np.sum(idcs_candidates) > 1:
                    He_cand, fluxes_cand = log.loc[idcs_candidates].index.values, log.loc[idcs_candidates].gauss_flux.values
                    ref_line = He_cand[np.argsort(fluxes_cand)[-2]]
                    idcs_lines = idcs_candidates
                else:
                    enougth_lines = False

        else:
            enougth_lines = True



        # Check if it is possible
        if enougth_lines:

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
                    Href_err = log.loc[ref_line, f'{ref_flux_type}_flux_err']

                # Use the user param
                else:
                    obsFlux = log.loc[idcs_lines, f'{flux_entry}_flux'].values
                    obsErr = log.loc[idcs_lines, f'{flux_entry}_flux_err'].values

                    Href_flux = log.loc[ref_line, f'{flux_entry}_flux']
                    Href_err = log.loc[ref_line, f'{flux_entry}_flux_err']

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

                        # rc.setCorr(obs_over_theo=5.34/2.86, wave1=6563., wave2=4861.)
                        # ratio_dis = np.random.normal(obsRatio_uarray[-1].nominal_value, obsRatio_uarray[-1].std_dev, size=1000)/2.86
                        # rc.setCorr(obs_over_theo=ratio_dis, wave1=6563., wave2=4861.)
                        # cHb, cHb_err = rc.cHbeta.mean(), rc.cHbeta.std() #(0.8395045358309076, 0.15112567990954212)
                        # eBV, eBV_err = rc.EbvFromCHbeta(cHb), rc.EbvFromCHbeta(cHb_err)
                        # print(f'E(B-V) = {eBV:0.2f} +/- {eBV_err:0.2f} || c(Hbeta) = {cHb:0.2f} +/- {cHb_err:0.2f}')


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
            _logger.info(f'Given the excluded lines the extinction coefficient could not be calculated')
            cHbeta, cHbeta_err = None, None

    else:
        _logger.info(f'{"Zero H1 lines" if line_list.size == 0 else "Just one H1 line"} in the input log, extinction coefficient '
                     f'could not be calculated')
        cHbeta, cHbeta_err = None, None

    # eBV, eBV_err = rc.EbvFromCHbeta(cHbeta), rc.EbvFromCHbeta(cHbeta_err)
    # print(f'E(B-V) = {eBV:0.2f} +/- {eBV_err:0.2f} || c(Hbeta) = {cHbeta:0.2f} +/- {cHbeta_err:0.2f}')

    return cHbeta, cHbeta_err


def flambda_calc(wavelength_array, R_v, red_curve, norm_wavelength):

    # Call pyneb
    rcGas = pn.RedCorr(R_V=R_v, law=red_curve)

    # Compute Xx parametrisation
    HbetaXx = rcGas.X(norm_wavelength)
    lineXx = rcGas.X(wavelength_array)

    # Flambda array
    f_lambda = lineXx/HbetaXx - 1.0

    return f_lambda


def reddening_correction(cHbeta, cHbeta_err, log, R_v=3.1, red_curve='G03 LMC', norm_wavelength=None, flux_column='gauss_flux',
                         n_points=1000, intensity_column='line_int'):

    # TODO log must be df only read from the log, Get normalization from log, add new column at front
    #

    # log = check_file_dataframe(log, DataFrame)

    line_wavelengths = log.wavelength.to_numpy()

    log['f_lambda'] = flambda_calc(line_wavelengths, R_v, red_curve, norm_wavelength)

    # Recover the parameters
    flux_array = log[f'{flux_column}'].to_numpy()
    err_array = log[f'{flux_column}_err'].to_numpy()
    f_lambda_array = log['f_lambda'].to_numpy()

    # Prepare distributions
    dist_size = (n_points, len(flux_array))
    flux_dist = np.random.normal(loc=flux_array, scale=err_array, size=dist_size)
    cHbeta_dist = np.random.normal(loc=cHbeta, scale=cHbeta_err, size=dist_size)

    # Compute the line intensities
    int_dist = flux_dist * np.power(10, cHbeta_dist * f_lambda_array)

    log[f'{intensity_column}'] = int_dist.mean(axis=0)
    log[f'{intensity_column}_err'] = int_dist.std(axis=0)

    # rc = pn.RedCorr(R_V=R_v, law=red_curve, cHbeta=cHbeta)
    # e_corr = rc.getCorr(log.wavelength.to_numpy(), rel_wave=norm_wavelength)
    # log['pyneb_int'] = log['line_flux'] * e_corr

    return


class ExtinctionModel:

    def __init__(self, Rv=None, red_curve=None, data_folder=None):

        self.R_v = Rv
        self.red_curve = red_curve

        # Dictionary with the reddening curves
        self.reddening_curves_calc = {'MM72': self.f_Miller_Mathews1972,
                                      'CCM89': self.X_x_Cardelli1989,
                                      'G03_bar': self.X_x_Gordon2003_bar,
                                      'G03_average': self.X_x_Gordon2003_average,
                                      'G03_supershell': self.X_x_Gordon2003_supershell}

        self.literatureDataFolder = data_folder

    def reddening_correction(self, wave, flux, err_flux=None, reddening_curve=None, cHbeta=None, E_BV=None, R_v=None, normWave=4861.331):

        # By default we perform the calculation using the colour excess
        if E_BV is not None:

            E_BV = E_BV if E_BV is not None else self.Ebv_from_cHbeta(cHbeta, reddening_curve, R_v)

            # Perform reddening correction
            wavelength_range_Xx = self.reddening_Xx(wave, reddening_curve, R_v)
            int_array = flux * np.power(10, 0.4 * wavelength_range_Xx * E_BV)

        else:
            lines_flambda = self.gasExtincParams(wave, R_v=R_v, red_curve=reddening_curve, normWave=normWave)

            if np.isscalar(cHbeta):
                int_array = flux * np.pow(10, cHbeta * lines_flambda)

            else:
                cHbeta = ufloat(cHbeta[0], cHbeta[1]),
                obsFlux_uarray = unumpy.uarray(flux, err_flux)

                int_uarray = obsFlux_uarray * unumpy.pow(10, cHbeta * lines_flambda)
                int_array = (unumpy.nominal_values(int_uarray), unumpy.std_devs(int_uarray))

        return int_array

    def Ebv_from_cHbeta(self, cHbeta, reddening_curve, R_v):

        E_BV = cHbeta * 2.5 / self.reddening_Xx(np.array([self.Hbeta_wavelength]), reddening_curve, R_v)[0]
        return E_BV

    def flambda_from_Xx(self, Xx, reddening_curve, R_v):

        X_Hbeta = self.reddening_Xx(np.array([self.Hbeta_wavelength]), reddening_curve, R_v)[0]

        f_lines = Xx / X_Hbeta - 1

        return f_lines

    def reddening_Xx(self, waves, curve_methodology, R_v):

        self.R_v = R_v
        self.wavelength_rc = waves
        return self.reddening_curves_calc[curve_methodology]()

    def f_Miller_Mathews1972(self):

        if isinstance(self.wavelength_rc, np.ndarray):
            y = 1.0 / (self.wavelength_rc / 10000.0)
            y_beta = 1.0 / (4862.683 / 10000.0)

            ind_low = np.where(y <= 2.29)[0]
            ind_high = np.where(y > 2.29)[0]

            dm_lam_low = 0.74 * y[ind_low] - 0.34 + 0.341 * self.R_v - 1.014
            dm_lam_high = 0.43 * y[ind_high] + 0.37 + 0.341 * self.R_v - 1.014
            dm_beta = 0.74 * y_beta - 0.34 + 0.341 * self.R_v - 1.014

            dm_lam = np.concatenate((dm_lam_low, dm_lam_high))

            f = dm_lam / dm_beta - 1

        else:

            y = 1.0 / (self.wavelength_rc / 10000.0)
            y_beta = 1.0 / (4862.683 / 10000.0)

            if y <= 2.29:
                dm_lam = 0.74 * y - 0.34 + 0.341 * self.R_v - 1.014
            else:
                dm_lam = 0.43 * y + 0.37 + 0.341 * self.R_v - 1.014

            dm_beta = 0.74 * y_beta - 0.34 + 0.341 * self.R_v - 1.014

            f = dm_lam / dm_beta - 1

        return f

    def X_x_Cardelli1989(self):

        x_true = 1.0 / (self.wavelength_rc / 10000.0)
        y = x_true - 1.82

        y_coeffs = np.array(
            [np.ones(len(y)), y, np.power(y, 2), np.power(y, 3), np.power(y, 4), np.power(y, 5), np.power(y, 6),
             np.power(y, 7)])
        a_coeffs = np.array([1, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999])
        b_coeffs = np.array([0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002])

        a_x = np.dot(a_coeffs, y_coeffs)
        b_x = np.dot(b_coeffs, y_coeffs)

        X_x = a_x + b_x / self.R_v

        return X_x

    def X_x_Gordon2003_bar(self):

        # Default R_V is 3.4
        R_v = self.R_v if self.R_v != None else 3.4  # This is not very nice
        x = 1.0 / (self.wavelength_rc / 10000.0)

        # This file format has 1/um in column 0 and A_x/A_V in column 1
        curve_address = os.path.join(self.literatureDataFolder, 'gordon_2003_SMC_bar.txt')
        file_data = np.loadtxt(curve_address)

        # This file has column
        Xx_interpolator = interp1d(file_data[:, 0], file_data[:, 1])
        X_x = R_v * Xx_interpolator(x)
        return X_x

    def X_x_Gordon2003_average(self):

        # Default R_V is 3.4
        R_v = self.R_v if self.R_v != None else 3.4  # This is not very nice
        x = 1.0 / (self.wavelength_rc / 10000.0)

        # This file format has 1/um in column 0 and A_x/A_V in column 1
        curve_address = os.path.join(self.literatureDataFolder, 'gordon_2003_LMC_average.txt')
        file_data = np.loadtxt(curve_address)

        # This file has column
        Xx_interpolator = interp1d(file_data[:, 0], file_data[:, 1])
        X_x = R_v * Xx_interpolator(x)
        return X_x

    def X_x_Gordon2003_supershell(self):

        # Default R_V is 3.4
        R_v = self.R_v if self.R_v != None else 3.4  # This is not very nice
        x = 1.0 / (self.wavelength_rc / 10000.0)

        # This file format has 1/um in column 0 and A_x/A_V in column 1
        curve_address = os.path.join(self.literatureDataFolder, 'gordon_2003_LMC2_supershell.txt')
        file_data = np.loadtxt(curve_address)

        # This file has column
        Xx_interpolator = interp1d(file_data[:, 0], file_data[:, 1])
        X_x = R_v * Xx_interpolator(x)
        return X_x

    def Epm_ReddeningPoints(self):

        x_true = np.arange(1.0, 2.8, 0.1)  # in microns -1
        X_Angs = 1 / x_true * 1e4

        Xx = np.array(
            [1.36, 1.44, 1.84, 2.04, 2.24, 2.44, 2.66, 2.88, 3.14, 3.36, 3.56, 3.77, 3.96, 4.15, 4.26, 4.40, 4.52,
             4.64])
        f_lambda = np.array(
            [-0.63, -0.61, -0.5, -0.45, -0.39, -0.34, -0.28, -0.22, -0.15, -0.09, -0.03, 0.02, 0.08, 0.13, 0.16, 0.20,
             0.23, 0.26])

        return x_true, X_Angs, Xx, f_lambda

    def gasExtincParams(self, wave, R_v = None, red_curve = None, normWave = 4861.331):

        if R_v is None:
            R_v = self.R_v
        if red_curve is None:
            red_curve = self.red_curve

        self.rcGas = pn.RedCorr(R_V=R_v, law=red_curve)

        HbetaXx = self.rcGas.X(normWave)
        lineXx = self.rcGas.X(wave)

        lineFlambda = lineXx / HbetaXx - 1.0

        return lineFlambda

    def contExtincParams(self, wave, Rv, reddening_law):

        self.rcCont = pn.RedCorr(R_V=Rv, law=reddening_law)

        lineXx = self.rcGas.X(wave)

        return lineXx