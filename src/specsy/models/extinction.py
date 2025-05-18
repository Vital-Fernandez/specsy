import logging
import numpy as np
from scipy.optimize import curve_fit

try:
    import pyneb as pn
    pyneb_check = True
except ImportError:
    pyneb_check = False

from pandas import DataFrame
from lime import label_decomposition
from lime.io import check_file_dataframe
from pathlib import Path
from lime.io import load_frame
from uncertainties import unumpy, ufloat
from lmfit.models import LinearModel
from specsy.plotting.plots import extinction_gradient
from specsy.tools import get_mixed_fluxes, linear_regression
from specsy.io import SpecSyError, check_file_dataframe

_logger = logging.getLogger('SpecSy')



# Function to compute and plot cHbeta
def extinction_coeff_calc(log, norm_line, R_V=3.1, law='G03 LMC', tem=10000.0, den=100.0, rel_Hbeta=True, flux_column='profile_flux',
                          line_list=None, exclude_list=None, exclude_nonphys=True, auto_ref_line=False, fname=None, plot_results=False,
                          fig_cfg=None, ax_cfg=None):

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

    :param tem: Temperature for the emissivity calculation in degrees Kelvin. The default value is 10000 K.
    :type tem: float, optional

    :param den: Density for the emissivity calculation in particles per centimeter cube. The default value is 100 cm^-3.
    :type den: float, optional

    :param norm_line: Line label of the normalization flux. The default value is "auto" for the automatic selection.
    :type norm_line: str, optional

    :param flux_column: Flux type for the cHbeta calculation. The default value is "gauss" for a Gaussian flux selection.
    :type flux_column: str, optional

    :param exclude_list: List of lines to exclude in the cHbeta calculation. The default value is None.
    :type exclude_list: list, optional

    :param plot_results: Check to display the cHbeta calculation regression. The default value is False.
    :type plot_results: bool, optional

    :param fname: Address for the output image with the cHbeta calculation regression. The default value is None.
    :type fname: str, optional

    :param plot_title: Title for the cHbeta calculation regression plot.
    :type plot_title: str, optional

    :param fig_cfg: Configuration for the cHbeta plot figure.
    :type fig_cfg: dict, optional

    :param ax_cfg: Configuration for the cHbeta plot axes.
    :type ax_cfg: dict, optional

    :return: cHbeta value and uncertainty.
    :rtype: float, float

    '''

    # By default return None values if calculation is not possible
    cHbeta, cHbeta_err = None, None

    # Read input log if it is a file
    log = check_file_dataframe(log)

    # If flux entry is there
    if flux_column not in log.columns:
        raise SpecSyError(f'The input lines table does not have a "{flux_column}" column')
    if f'{flux_column}_err' not in log.columns:
        raise SpecSyError(f'The input lines table does not have a "{flux_column}_err" column')

    # Read the hydrogen lines
    if line_list is None:
        if 'particle' in log.columns:
            idcs_H1 = (log.particle == 'H1') & ~log.index.str.contains('_k-')
        else:
            idcs_H1 = log.index.str.startswith('H1') & ~log.index.str.contains('_k-')
        line_list = log.loc[idcs_H1].index.to_numpy()
    else:
        line_list = np.ravel(np.asarray(line_list))

    # Exclude non-physical entries
    if exclude_nonphys:
        idcs_valid = (log.loc[line_list, flux_column].notna() & (log.loc[line_list, flux_column] > 0) & np.isfinite(log.loc[line_list, flux_column]) &
                      log.loc[line_list, f'{flux_column}_err'].notna() & np.isfinite(log.loc[line_list, f'{flux_column}_err']))

        line_list = line_list[idcs_valid]

    # Exclude lines if requested by the user
    idcs_valid = np.ones(line_list.size).astype(bool) if exclude_list is None else ~np.isin(line_list, exclude_list)

    # Proceed if there are enough lines to compute the extinction
    if line_list.size > 1:

        # Establish the reference line
        if norm_line not in line_list:
            if auto_ref_line is True:
                norm_line = log.loc[line_list, flux_column].nlargest(2).index[1]
            else:
                raise SpecSyError(f'The input normalization line "{norm_line}" was not found in the input lines list or frame. '
                                  f'You may set the argument "auto_ref_line=True" to automatically choose the second most '
                                  f'intense hydrogen line as the normalization flux.')

        # Recover the normalization flux
        norm_flux, norm_err = log.loc[norm_line, [flux_column, f'{flux_column}_err']]
        norm_flux = None if np.isnan(norm_flux) else norm_flux
        norm_err = None if np.isnan(norm_err) else norm_err

        if norm_flux is not None or norm_err is not None:

            # Get line properties
            params = ['wavelength', 'latex_label']
            waves_norm, latex_norm = label_decomposition(norm_line, params_list=params, scalar_output=True)
            waves_array, latex_array = label_decomposition(line_list, params_list=params)

            flux_array = log.loc[line_list, flux_column].to_numpy()
            err_array = log.loc[line_list, f'{flux_column}_err'].to_numpy()

            # Compute flux ratios
            obs_ratios = flux_array / norm_flux
            obs_errs = obs_ratios * np.sqrt(np.square(err_array / flux_array) + np.square(norm_err / norm_flux))

            # Compute the theoretical ratios
            H1 = pn.RecAtom('H', 1)
            emis_norm = H1.getEmissivity(tem=tem, den=den, wave=waves_norm)
            emis_arr = np.fromiter((H1.getEmissivity(tem=tem, den=den, wave=wave) for wave in waves_array), float)
            theo_ratios = emis_arr / emis_norm

            # Reddening law
            rc = pn.RedCorr(R_V=R_V, law=law)
            Xx_ref, Xx = rc.X(waves_norm), rc.X(waves_array)
            f_lines, f_ref = Xx/Xx_ref - 1, Xx_ref/Xx_ref - 1

            # Linear fit for the extinction calculation
            x_arr = f_lines - f_ref
            y_arr = np.log10(theo_ratios) - np.log10(obs_ratios)
            y_err = (1 / np.log(10)) * (obs_errs / obs_ratios)
            c_HI, c_HI_err, intercept, intercept_err = linear_regression(x_arr[idcs_valid], y_arr[idcs_valid], y_err[idcs_valid])

            # Correction for non-Hbeta normalization
            if rel_Hbeta:
                conv_factor =  (1 - f_ref) * (rc.X(4861.25) / Xx_ref)
                cHbeta = conv_factor * c_HI
                cHbeta_err = conv_factor * c_HI_err

            else:
                cHbeta, cHbeta_err = c_HI, c_HI_err

            # Store the results in the lines frame
            log.loc[line_list, 'f_lambda'] = f_lines
            log.loc[line_list, 'theo_ratio'] = theo_ratios
            log.loc[line_list, 'obs_ratio'] = obs_ratios
            log.loc[line_list, 'obs_ratio_err'] = obs_errs

            log.loc[line_list[idcs_valid], 'extinction_idcs'] = 1
            log.loc[line_list[~idcs_valid], 'extinction_idcs'] = 2
            log.loc[norm_line, 'extinction_idcs'] = 0

            # Generate the plot
            if plot_results:
                extinction_gradient(cHbeta, cHbeta_err, log, rel_Hbeta=rel_Hbeta, fname=fname, fig_cfg=fig_cfg, ax_cfg=ax_cfg)

        else:
            _logger.info(f'No valid "{flux_column}" values/uncertainties for the normalization line "{norm_line}".')

    else:
        _logger.info(f'{"No H1 lines" if line_list.size == 0 else "Only one H1 line"} available in the input lines list or measurements frame. '
                     f'extinction coefficient could not be calculated')

    return cHbeta, cHbeta_err, log


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
    # log[f'{intensity_column}'] = int_dist.mean(axis=0)
    # log[f'{intensity_column}_err'] = int_dist.std(axis=0)
    log.insert(0, f'{intensity_column}', int_dist.mean(axis=0))
    log.insert(1, f'{intensity_column}_err', int_dist.std(axis=0))

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