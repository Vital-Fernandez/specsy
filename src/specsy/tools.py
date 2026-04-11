import logging
from dataclasses import dataclass, field
from typing import Callable, Optional, Union

import numpy as np
from pandas import DataFrame
from scipy.stats import truncnorm, norm
from scipy.optimize import curve_fit
from lime import Line

# TODO some of these could go to lime

_logger = logging.getLogger('SpecSy')


@dataclass
class RegionParam:
    """
    Defines how a single variable (temperature or density) is specified.

    - "free":  provide distribution and kwargs. The variable is sampled.
    - "tied":  provide ref (name or index) to link to another region.
               Optionally provide an equation key (string) to transform the
               parent value. If no equation is provided, the value is copied as-is.
    """

    label: str

    mode: str  # "free" | "tied"

    # Distribution for free variable and its parameters
    distr: Optional[Callable] = None
    kwargs: dict = field(default_factory=dict)

    # Parent region for the parameter and key for the empirical/theoretical relation if necessary
    ref: Optional[Union[str, int]] = None
    eq: Optional[str] = None  # key into equations dict


class Region:

    def __init__(self, name, species, temp_mode, temp_dist=None, temp_kwargs=None, temp_ref=None, temp_eq=None,
                 den_mode=None, den_dist=None, den_kwargs=None, den_ref=None, den_eq=None, ):
        self.name = name
        self.species = species
        self.temp = RegionParam('temp', temp_mode, distr=temp_dist, kwargs=temp_kwargs, ref=temp_ref, eq=temp_eq)
        self.den = RegionParam('den', den_mode, distr=den_dist, kwargs=den_kwargs, ref=den_ref, eq=den_eq)

        return


class MultiRegionModel:

    def __init__(self, regions: list[Region],
                 temp_equations_dict: dict[str, Callable] = None,
                 den_equations_dict: dict[str, Callable] = None):

        assert 1 <= len(regions) <= 4, "Between 1 and 4 regions allowed"
        self.regions = regions
        self.size = len(regions)
        # self.temp_eq_dict = temp_equations_dict or _TEMP_FUNC
        # self.den_eq_dict = den_equations_dict or _DEN_FUNC
        self.region_map: dict = {r.name: r for r in regions}
        self.sampled: dict = {}

    def build(self):

        # Sample the distribution value if free
        for region in self.regions:
            self.sampled[region.name] = {}

            if region.temp.mode == "free":
                self.sampled[region.name]["temp"] = region.temp.distr(**region.temp.kwargs)

            if region.den.mode == "free":
                self.sampled[region.name]["den"] = region.den.distr(**region.den.kwargs)

        return self

    def get(self, region_name: str, var_type: str):

        # Recover the requested parameter for the physical region
        region = self.region_map[region_name]
        var_spec = getattr(region, var_type)
        equations = getattr(self, f"{var_type}_eq_dict")

        match var_spec.mode:

            # Free variable
            case "free":
                return self.sampled[region_name][var_type]

            # Tied to another example
            case 'tied':
                ref_name = self.regions[var_spec.ref].name if isinstance(var_spec.ref, int) else var_spec.ref
                parent_value = self.get(ref_name, var_type)  # recursive call

                if var_spec.eq is not None:
                    return equations[var_spec.eq](parent_value)
                else:
                    return parent_value

            case _:
                raise ValueError(f"Unknown mode '{var_spec.mode}'. Use 'free' or 'tied'.")

    def map_line_structure(self, line_list, temp_label='temp', den_label='den'):

        # Container to tabulate the data and conditions of the table
        df = DataFrame(index=line_list,
                          columns=["particle", "region", temp_label, den_label, "eq_temp", "eq_den"])

        # Add particles to the table
        line_list = Line.from_list(line_list)
        df['particle'] = [line.particle.label for line in line_list]

        # Map the regions
        for region in self.regions:
            idcs = df.particle.isin(region.species)
            df.loc[idcs, 'region'] = region.name

        # Recover the temperature and density structure
        for region in self.regions:

            # Assign the free params
            idcs = df.region == region.name
            for param_label in [temp_label, den_label]:
                param = getattr(region, param_label)
                if param.mode == "free":
                    df.loc[idcs, param_label] = f'{param.label}_{region.name}'
                else:
                    df.loc[idcs, param_label] = f'{param.label}_{param.ref}'
                    if param.eq is not None:
                        df.loc[idcs, f'eq_temp'] = param.eq

        return df



# Percentile notation for uncertainty
def percentile_latex_uncertainty(median, superscript, subscript, sig_fig=2):
    return r'${}^{{{}}}_{{{}}}$'.format(np.round(median, sig_fig), np.round(superscript, sig_fig), np.round(subscript, sig_fig))


# Scipy formula for truncation coefficient
def truncation_limits(mu, sigma, lower_limit, upper_limit):
    return (lower_limit - mu) / sigma, (upper_limit - mu) / sigma


# Function to generate a truncated normal function
def truncated_gaussian(diag_int, diag_err, n_steps, low_limit=-np.inf, up_limit=np.inf):
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
        obsFlux, obsErr = log[f'{flux_type}'].values, log[f'{flux_type}_err'].values
        # if flux_type in ['intg', 'profile']:
        #     obsFlux, obsErr = log[f'{flux_type}_flux'].values, log[f'{flux_type}_err'].values
        # else:
        #     _logger.warning(f'The flux type {flux_type} is not recognized. Please use "intg" or "profile" for integrated '
        #                     f'or gaussian fluxes respectively')
        #     raise ValueError

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


def linear_model(x, m_cont, n_cont):
    return m_cont * x + n_cont


def linear_regression(x_values, y_values, y_error):

    params, covariance = curve_fit(linear_model, xdata=x_values, ydata=y_values, sigma=y_error, absolute_sigma=True,
                                   check_finite=False)

    m, n = params
    m_err, n_err = np.sqrt(np.diag(covariance))

    return m, m_err, n, n_err