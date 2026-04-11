from dataclasses import dataclass, field
from typing import Callable, Optional, Union
import numpy as np
import pandas as pd
from lime import Line
from pytensor import tensor, function as pt_function
from innate import load_dataset
from scipy.interpolate import RegularGridInterpolator
import pyneb as pn
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
import innate

def double_temp(temp: float) -> float:
    return temp * 2

def half_temp(temp: float) -> float:
    return temp / 2

def square_density(den: float) -> float:
    return den ** 2

# --- Empirical relations ---
def temp_OIII(temp: float) -> float:
    return temp - 1500


def den_OII(den: float) -> float:
    return den * 2


# store them in a dictionary
_TEMP_FUNC = {"double_temp": double_temp, 'TSIII_hagele': temp_OIII, 'half_temp': half_temp}
_DEN_FUNC = {"square_density": square_density, 'nOII_stasinska': den_OII}
_PRIOR_PARAM = {'temp_low': pm.Normal, 'temp_med': pm.Normal, 'temp_high': pm.Normal, 'temp_vhigh': pm.Normal,
                'den_low': pm.HalfNormal, 'den_med': pm.HalfNormal, 'den_high': pm.HalfNormal, 'den_vhigh': pm.HalfNormal}
_DIST_PARAM = {'temp_low': dict(mu=15000, sigma=2000), 'temp_med': dict(mu=15000, sigma=2000), 'temp_high': dict(mu=15000, sigma=2000), 'temp_vhigh': dict(mu=15000, sigma=2000),
                'den_low': dict(sigma=300), 'den_med':  dict(sigma=300), 'den_high':  dict(sigma=300), 'den_vhigh':  dict(sigma=300)}

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
                                      den_mode=None, den_dist=None, den_kwargs=None, den_ref=None, den_eq=None,):

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
        self.temp_eq_dict = temp_equations_dict or _TEMP_FUNC
        self.den_eq_dict = den_equations_dict or _DEN_FUNC
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
        df = pd.DataFrame(index=line_list,
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


def generate_emissivity_file(fname, lines_df):

    # Compute the ranges:
    temp_range = np.linspace(9000, 20000, 251)
    den_range = np.linspace(1, 600, 101)

    # Normalization for the emissivities
    H1 = pn.RecAtom('H', 1)
    norm_emiss = H1.getEmissivity(temp_range, den_range, wave=4861.0)

    # Container output grids
    atom_dict, emiss_dict = {}, {}

    # Loop through the lines and compute the emissivities
    wavelength, particle, t_type = lines_df[['wavelength', 'particle', 'trans']].to_numpy().T
    for i, line_name in enumerate(lines_df.index):

        # Get transition atom pyneb object
        print(f'-- {line_name}')
        elem, ionization = particle[i][:-1], particle[i][-1]
        atom = pn.RecAtom(elem, ionization) if t_type[i] == 'rec' else pn.Atom(elem, ionization)

        # Compute and normalize the emissivities
        grid_i = atom.getEmissivity(temp_range, den_range, wave=np.round(wavelength[i]))
        emiss_dict[line_name] = grid_i/norm_emiss

    # Data attributes
    data_conf = {'parameter': 'emissivity', 'approximation': ('rgi', 'eqn'), 'axes': ('temp', 'den'),
                 'temp_range': (9000, 20000, 251), 'den_range': (1, 600, 101)}

    # Save the data into a dictionary
    innate.save_dataset(fname, emiss_dict, data_conf, custom_cfg=None)

    return

def make_bilinear_interp(x_grid, y_grid, z_grid):

    """
    Returns a pytensor-differentiable bilinear interpolation function.

    x_grid: 1D np.array, shape (N,)
    y_grid: 1D np.array, shape (M,)
    z_grid: 2D np.array, shape (N, M)

    """

    dx = x_grid[1] - x_grid[0]
    dy = y_grid[1] - y_grid[0]
    x0 = x_grid[0]
    y0 = y_grid[0]
    N = len(x_grid) - 2
    M = len(y_grid) - 2
    z = tensor.as_tensor_variable(z_grid)  # baked into the graph, not sampled

    def interp(x, y):
        i = tensor.clip(tensor.cast(tensor.floor((x - x0) / dx), "int32"), 0, N)
        j = tensor.clip(tensor.cast(tensor.floor((y - y0) / dy), "int32"), 0, M)

        tx = (x - (x0 + i * dx)) / dx
        ty = (y - (y0 + j * dy)) / dy

        return (z[i, j] * (1 - tx) * (1 - ty) +
                z[i + 1, j] * tx * (1 - ty) +
                z[i, j + 1] * (1 - tx) * ty +
                z[i + 1, j + 1] * tx * ty)

    return interp


def compile_emis_interp(fname, array_mode=False):

    emis_set = load_dataset(fname)

    interp_dict = {}
    for trans, emis_matrix in emis_set[0].items():

        temp_range = np.linspace(*emis_set[1][trans]['temp_range'])
        den_range = np.linspace(*emis_set[1][trans]['den_range'])
        log_emis = np.log10(emis_matrix)

        interp = make_bilinear_interp(temp_range, den_range, log_emis)

        if array_mode:
            x_sym = tensor.dscalar("x")
            y_sym = tensor.dscalar("y")
            interp_dict[trans] = pt_function([x_sym, y_sym], interp(x_sym, y_sym))

        else:
            interp_dict[trans] = interp

    return interp_dict


def compile_line_fluxes(line_list, f_lambda, emis_interp, params, verbose=False):

    flux_arr = np.zeros(len(line_list))
    flux_log_arr = np.zeros(len(line_list))
    for i, label in enumerate(line_list):
        line = Line.from_transition(label)
        ion = line.particle.label

        # Return parameter
        abund_log = params.get(ion, 0)
        abund = np.power(10, abund_log - 12) if abund_log != 0 else 1
        cHbeta = params['cHbeta']
        temp = params['temp']
        den = params['den']

        # Compute emis
        log_emis = emis_interp[label](temp, den)
        emis = np.power(10, log_emis)

        # Flux computation
        if ion == 'H1':
            flux_i = abund * emis * np.power(10, -f_lambda[i] * cHbeta)
            flux_log = log_emis - f_lambda[i] * cHbeta

        else:
            flux_i = abund * emis * np.power(10, -f_lambda[i] * cHbeta)
            flux_log = abund_log + log_emis - f_lambda[i] * cHbeta - 12

        flux_arr[i] = flux_i
        flux_log_arr[i] = flux_log

        if verbose:
            # O2 = pn.Atom('O', 2)
            # H1 = pn.RecAtom('H', 1)
            # emis = O2.getEmissivity(temp, den, wave=3725.974)/H1.getEmissivity(temp, den, wave=4861)
            print(f'{i}) {label}: {ion} = {abund_log}, emis = {emis:0.3f},'
                  f' f_lambda = {f_lambda[i]:0.3f}, '
                  f'flux={flux_log_arr[i]:0.3f} '
                  f'(match {np.isclose(flux_i, np.power(10, flux_log_arr[i]))})')

    return flux_arr


def emission_flux_model(line_list, chemical_params, struc_df, emis_interp, red_law='CCM89', R_V=3.1, f_lambda_ref='H1_4861A',
                        verbose=False):

    # Container to store the line parameters
    cond_dict = {}

    # Get the abundances
    lime_list = Line.from_list(line_list)
    for line in lime_list:
        cond_dict[line.particle.label] = chemical_params[line.particle.label]

    # Compile extinctions
    rc = pn.RedCorr(law=red_law, R_V=R_V)
    wave_arr = np.array([line.wavelength for line in lime_list])
    f_lambda = rc.X(wave_arr) / rc.X(Line.from_transition(f_lambda_ref).wavelength) - 1.0
    cond_dict['cHBeta'] = chemical_params['cHBeta']

    # Get free temperatures and densities
    print(struc_df.to_string())
    for param in pd.unique(struc_df[['temp', 'den']].values.ravel()):
        cond_dict[param] = chemical_params[param]

    # Get the temperature and density equation list with nan for false values
    line_arr = struc_df.index.to_numpy()
    Tem_label_arr = struc_df.temp.to_numpy()
    den_label_arr = struc_df.den.to_numpy()
    tem_eq_arr = struc_df.eq_temp.to_numpy()
    den_eq_arr = struc_df.eq_den.to_numpy()
    particle_arr = struc_df.particle.to_numpy()

    flux_arr = np.zeros(len(line_list))

    # Convenience arrays for the workflow
    range_arr = np.arange(len(line_list))
    temp_eq_check = struc_df.eq_temp.isnull().to_numpy()
    den_eq_check = struc_df.eq_den.isnull().to_numpy()

    for i in range_arr:

        cHbeta = cond_dict['cHBeta']

        # Compute the emissivity
        tem = cond_dict[Tem_label_arr[i]] if temp_eq_check[i] else _TEMP_FUNC[tem_eq_arr[i]](cond_dict[Tem_label_arr[i]])
        den = cond_dict[den_label_arr[i]] if den_eq_check[i] else _TEMP_FUNC[den_eq_arr[i]](cond_dict[den_label_arr[i]])
        emis = emis_interp[line_arr[i]](tem, den)

        # Compute the flux
        if particle_arr[i] == 'H1':
            flux = emis - f_lambda[i] * cHbeta
        else:
            # Use the local Python variable directly, not dm_model[]
            flux = cond_dict[particle_arr[i]] + emis - f_lambda[i] * cHbeta - 12

        flux_arr[i] = flux

    # Convert to linear scale
    flux_arr = np.power(10, flux_arr)

    return flux_arr


def direct_method_multi_region(lines_df, emis_interp, true_params=None):

    # Convert the fluxes to the log scale
    flux_arr, err_arr = lines_df[['line_flux', 'line_err']].to_numpy().T
    input_obs = np.log10(flux_arr)
    input_err = np.log10(1 + err_arr / flux_arr)

    # Unpack physical parameters
    flambda_arr = lines_df.f_lambda.to_numpy()
    particle_arr = lines_df.particle.to_numpy()
    ion_species = pd.unique(particle_arr)

    # Unpack the temp/den structure arrays
    Tem_label_arr = lines_df.temp.to_numpy()
    den_label_arr = lines_df.den.to_numpy()
    tem_eq_arr = lines_df.eq_temp.to_numpy()
    den_eq_arr = lines_df.eq_den.to_numpy()
    unique_params = pd.unique(lines_df[['temp', 'den']].values.ravel())

    # Convenience arrays for the workflow
    line_arr = lines_df.index.to_numpy()
    range_arr = np.arange(len(line_arr))
    temp_eq_check = lines_df.eq_temp.isnull().to_numpy()
    den_eq_check = lines_df.eq_den.isnull().to_numpy()

    with pm.Model(coords={"lines": line_arr}) as dm_model:

        # Container to store the models
        theo_flux = tensor.zeros(line_arr.size)

        # Compile the abundances
        for ion in ion_species:
            pm.Normal(name=ion, mu=5, sigma=5)

        # Extinction
        cHbeta = pm.HalfCauchy(name="cHBeta", beta=2)

        # Generate the free temperatures and densities
        for param in unique_params:
            _PRIOR_PARAM[param](param, **_DIST_PARAM[param])

        for i in range_arr:

            # Compute the emissivity
            tem = dm_model[Tem_label_arr[i]] if temp_eq_check[i] else _TEMP_FUNC[tem_eq_arr[i]](dm_model[Tem_label_arr[i]])
            den = dm_model[den_label_arr[i]] if den_eq_check[i] else _DEN_FUNC[den_eq_arr[i]](dm_model[den_label_arr[i]])
            emis = emis_interp[line_arr[i]](tem, den)

            if particle_arr[i] == 'H1':
                flux = emis - flambda_arr[i] * cHbeta
            else:
                flux = dm_model[particle_arr[i]] + emis - flambda_arr[i] * cHbeta - 12

            theo_flux = tensor.inc_subtensor(theo_flux[i], flux)

        # Stored the fluxes
        pm.Deterministic('theo_flux', theo_flux)

        # Likelihood
        pm.Normal("likelihood", mu=theo_flux, sigma=input_err, observed=input_obs, dims='lines')


    # Run the model
    with dm_model:
        trace = pm.sample(draws=1000, tune=2000, target_accept=0.9, chains=8, cores=8, nuts_sampler='numpyro')

    # Output the results
    var_names = ["O2", "O3", "S2", "S3", "N2", "Ar3", "Ar4", "Ne3", "cHBeta"] + list(unique_params)
    az.plot_pair(trace, var_names=var_names, divergences=True)
    az.plot_posterior(trace, var_names=var_names)
    summary = az.summary(trace, var_names=var_names)
    print(summary)
    plt.show()

    return




# @dataclass
# class RegionParam:
#
#     """
#     Defines how a single variable (temperature or density) is specified.
#
#     - "free":  provide distribution and kwargs. The variable is sampled.
#     - "tied":  provide ref (name or index) to link to another region.
#                Optionally provide an equation key (string) to transform the
#                parent value. If no equation is provided, the value is copied as-is.
#     """
#
#     mode: str  # "free" | "tied"
#
#     # Distribution for free variable and its parameters
#     distr: Optional[Callable] = None
#     kwargs: dict = field(default_factory=dict)
#
#     # Parent region for the parameter and key for the empirical/theoretical relation if necessary
#     ref: Optional[Union[str, int]] = None
#     eq: Optional[str] = None  # key into equations dict
#
#
#
# class Region:
#
#     def __init__(self, name, species, temp_mode, temp_dist=None, temp_kwargs=None, temp_ref=None, temp_eq=None,
#                                       den_mode=None, den_dist=None, den_kwargs=None, den_ref=None, den_eq=None,):
#
#         self.name = name
#         self.species = species
#         self.temp = RegionParam(temp_mode, distr=temp_dist, kwargs=temp_kwargs, ref=temp_ref, eq=temp_eq)
#         self.den = RegionParam(den_mode, distr=den_dist, kwargs=den_kwargs, ref=den_ref, eq=den_eq)
#
#         return
#
#
#
# class MultiRegionModel:
#
#     def __init__(self, regions: list[Region],
#                  temp_equations_dict: dict[str, Callable] = None,
#                  den_equations_dict: dict[str, Callable] = None):
#
#         assert 1 <= len(regions) <= 4, "Between 1 and 4 regions allowed"
#         self.regions = regions
#         self.size = len(regions)
#         self.temp_eq_dict = temp_equations_dict or _TEMP_FUNC
#         self.den_eq_dict = den_equations_dict or _DEN_FUNC
#         self.region_map: dict = {r.name: r for r in regions}
#         self.sampled: dict = {}
#
#     def build(self):
#
#         # Sample the distribution value if free
#         for region in self.regions:
#             self.sampled[region.name] = {}
#
#             if region.temp.mode == "free":
#                 self.sampled[region.name]["temp"] = region.temp.distr(**region.temp.kwargs)
#
#             if region.den.mode == "free":
#                 self.sampled[region.name]["den"] = region.den.distr(**region.den.kwargs)
#
#         return self
#
#     def get(self, region_name: str, var_type: str):
#
#         # Recover the requested parameter for the physical region
#         region = self.region_map[region_name]
#         var_spec = getattr(region, var_type)
#         equations = getattr(self, f"{var_type}_eq_dict")
#
#         match var_spec.mode:
#
#             # Free variable
#             case "free":
#                 return self.sampled[region_name][var_type]
#
#             # Tied to another example
#             case 'tied':
#                 ref_name = self.regions[var_spec.ref].name if isinstance(var_spec.ref, int) else var_spec.ref
#                 parent_value = self.get(ref_name, var_type)  # recursive call
#
#                 if var_spec.eq is not None:
#                     return equations[var_spec.eq](parent_value)
#                 else:
#                     return parent_value
#
#             case _:
#                 raise ValueError(f"Unknown mode '{var_spec.mode}'. Use 'free' or 'tied'.")
#
#     def map_line_structure(self, line_list, temp_prefix='temp_', den_prefix='den_'):
#
#         # Container to tabulate the data and conditions of the table
#         df = pd.DataFrame(index=line_list,
#                           columns=["line", "particle", "region", "temp", "den", "eq_temp", "eq_den"])
#
#         # Add particles to the table
#         line_list = Line.from_list(line_list)
#         df['particle'] = [line.particle.label for line in line_list]
#
#         # Map the regions
#         for region in self.regions:
#             idcs = df.particle.isin(region.species)
#             df.loc[idcs, 'region'] = region.name
#
#         # Recover the temperature and density structure
#         for region in self.regions:
#
#             # Assign the free params
#             idcs = df.region == region.name
#             for param in ['temp', 'den']:
#                 if getattr(region, param).mode == "free":
#                     df.loc[idcs, param] = f'{param}_{region.name}'
#                 else:
#                     df.loc[idcs, 'temp'] = f'{param}_{region.ref}'
#                     df.loc[idcs, f'eq_temp'] = region.eq
#
#         return df
