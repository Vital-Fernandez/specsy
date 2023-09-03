import logging
import numpy as np
from scipy.stats import truncnorm, norm
from lime.tools import extract_fluxes, normalize_fluxes

# TODO some of these could go to lime

_logger = logging.getLogger('SpecSy')

# Percentile notation for uncertainty
def percentile_latex_uncertainty(median, superscript, subscript, sig_fig=2):
    return r'${}^{{{}}}_{{{}}}$'.format(np.round(median, sig_fig), np.round(superscript, sig_fig), np.round(subscript, sig_fig))


# Scipy formula for truncation coefficient
def truncation_limits(mu, sigma, lower_limit, upper_limit):
    return (lower_limit - mu) / sigma, (upper_limit - mu) / sigma


# Function to generate a truncated normal function
def truncated_gaussian(diag_int, diag_err, n_steps, low_limit=-np.infty, up_limit=np.infty):
    a, b = truncation_limits(diag_int, diag_err, low_limit, up_limit)
    output_dist = truncnorm.rvs(a, b, loc=diag_int, scale=diag_err, size=n_steps)
    return output_dist


# Function to get the lines which are blended form LiMe log
def blended_label_from_log(log):

    idcs_blended = (log['profile_label'] != 'no') & (~log.index.str.endswith('_m'))

    return idcs_blended.values


# Favoured method to get line fluxes according to resolution
def get_mixed_fluxes(log):

    # Get indeces of blended lines
    idcs_blended = blended_label_from_log(log)

    # First create full arrays with integrated fluxes
    obsFlux = log['intg_flux'].values
    obsErr = log['intg_err'].values

    # Then assign gaussian fluxes to blended
    if np.any(idcs_blended):
        obsFlux[idcs_blended] = log.loc[idcs_blended, 'gauss_flux'].values
        obsErr[idcs_blended] = log.loc[idcs_blended, 'gauss_err'].values

    return obsFlux, obsErr


# Flux_distribution for Monte-Carlo error propagation
def flux_distribution(log, flux_type='auto', n_steps=1000):

    if flux_type == 'auto':
        obsFlux, obsErr = get_mixed_fluxes(log)
    else:
        if flux_type in ['intg', 'gauss']:
            obsFlux, obsErr = log[f'{flux_type}_flux'].values, log[f'{flux_type}_err'].values
        else:
            _logger.warning(f'The flux type {flux_type} is not recognized. Please use "intg" or "gauss" for integrated '
                            f'or gaussian fluxes respectively')
            raise ValueError

    # Generate a normal distribution for every line
    output_dict = {}
    for i, line in enumerate(log.index.values):

        if not np.isnan(obsFlux[i]) and obsFlux[i] > 0:

            if not np.isnan(obsErr[i]) and obsErr[i] > 0:

                output_dict[line] = np.random.normal(obsFlux[i], obsErr[i], n_steps)

            else:
                _logger.info(f'Invalid {line} err ({obsErr[i]}). It is excluded from the distribution dictionary')

        else:
            _logger.info(f'Invalid {line} flux ({obsFlux[i]}). It is excluded from the distribution dictionary')

    return output_dict