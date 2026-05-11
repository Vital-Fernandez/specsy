import logging
import warnings
from dataclasses import dataclass

import lime
import numpy as np
import arviz as az
import xarray as xr
from lime import label_decomposition, Line, lines_frame, normalize_fluxes
from pathlib import Path
from innate import load_dataset

from specsy.io import SpecSyError, specsy_cfg
from specsy.tools import truncated_gaussian, flux_distribution
from specsy.operations.interpolation import compile_bilinear_interp
from specsy.models.extinction import flambda_calc
from specsy.models.literature import _TEM_FUNC_DICT, _DEN_FUNC_DICT
from specsy.models.fluxes_line import DEFAULT_PARTICLE_EQUATIONS_KEYS, FLUX_EQUATION_DICT
from specsy.sampler import direct_method_multi_region, run_model
from matplotlib import pyplot as plt

try:
    import pyneb as pn
    pyneb_check = True
except ImportError:
    pyneb_check = False


_logger = logging.getLogger('SpecSy')


def truncated_SII_density_dist(log=None, SII_lines=('S2_6716A', 'S2_6731A'), temp=10000, S2_pyneb=None, flux_dict=None,
                               n_steps=1000):

    '''

    This function computes the electron density from the [SII]6716,6731A doublet. The line label must adhere to the
    LiMe format.

    The user can input a pandas dataframe lines log. This log should adhere to LiMe formatting. Alternatively, the
    user can provide a dictionary with the lines flux distributions. The keys in this dictionary should be the same as
    in the "SII_lines" argument.

    The emissivity calculation is done using PyNeb. The user can provide its own "S2" Atom object. Otherwise, one is created
    with the default PyNeb atomic data.

    The output density distribution is truncated to avoid values outside the physical emissivity ratios.

    :param log: Lines log with the input fluxes. The pandas dataframe must adhere to LiMe formatting
    :type log: pd.DataFrame

    :param SII_lines: Tupple with the label for the [SII] lines. The default values are ('S2_6716A','S2_6731A')
    :type SII_lines: tuple, optional

    :param temp: Temperature for the emissivity calculation in degrees Kelvin. The default value is 10000 K.
    :type temp: float, optional

    :param S2_pyneb: Pyneb Atom, Atom for the S^+ ion.
    :type S2_pyneb: pyneb.Atom, optional

    :param flux_dict: Dictionary with the flux distribution for the [SII] lines.
    :type flux_dict: dict, optional

    :param n_steps: Number of steps in the Monte-Carlo sampling (only if flux_dict is not provided). The default value is 1000.
    :type n_steps: float, optional

    :return: [SII] electron density distribution.
    :rtype: np.array

    '''

    if flux_dict is None:
        flux_dict = flux_distribution(log, 'auto')

    # Compute the densities
    if (SII_lines[0] in flux_dict) and (SII_lines[1] in flux_dict):

        S2_ratio = flux_dict[SII_lines[0]]/flux_dict[SII_lines[1]]

        RSII, RSII_err = np.mean(S2_ratio), np.std(S2_ratio)
        RSII_dist = truncated_gaussian(RSII, RSII_err, n_steps, low_limit=0.28, up_limit=1.42)

        S2 = S2_pyneb if S2_pyneb is not None else pn.Atom('S', 2)
        neSII_dist = S2.getTemDen(RSII_dist, tem=temp, to_eval='L(6716)/L(6731)')

        if np.any(np.isnan(neSII_dist)):
            _logger.warning(f'ne_[SII] distribution contains nan entries')

    else:
        _logger.info('Both [SII] doublet not found in log, the density was not be calculated')
        neSII_dist = None

    return neSII_dist


def ratio_S23(flux_dict, S2_lines=('S2_6716A', 'S2_6731A'), S3_lines=('S3_9068A', 'S3_9530A'),
              norm_lines=('H1_6563A', 'H1_9546A'), H1_pyneb=None, temp=10000, den=100):

    S_23 = None
    if (S2_lines[0] in flux_dict) and (S2_lines[1] in flux_dict):
        if (S3_lines[0] in flux_dict) and (S3_lines[1] in flux_dict):
            if (norm_lines[0] in flux_dict) and (norm_lines[1] in flux_dict):

                H1_pyneb = H1_pyneb if H1_pyneb is not None else pn.RecAtom('H', 1)

                ion_norm1, norm_lines1, latex_norm1 = label_decomposition(norm_lines[0])
                ion_norm2, norm_lines2, latex_norm2 = label_decomposition(norm_lines[1])

                Hbeta_emis = H1_pyneb.getEmissivity(temp, den, wave=4861)
                S2_norm = H1_pyneb.getEmissivity(temp, den, wave=norm_lines1[0])/Hbeta_emis
                S3_norm = H1_pyneb.getEmissivity(temp, den, wave=norm_lines2[0])/Hbeta_emis

                S_2 = (flux_dict[S2_lines[0]]+flux_dict[S2_lines[1]])/flux_dict[norm_lines[0]] * S2_norm
                S_3 = (flux_dict[S3_lines[0]]+flux_dict[S3_lines[1]])/flux_dict[norm_lines[0]] * S3_norm

                S_23 = S_2 + S_3

            else:
                warn_text = f'{norm_lines[0]} ' if norm_lines[0] not in flux_dict else ""
                warn_text += f'{norm_lines[1]}' if norm_lines[1] not in flux_dict else ""
                _logger.info(f'The normalization lines {warn_text} are missing. Please provide others for S23 calculation')
        else:
            _logger.info('[SIII] lines missing. S_23 could not be calculated')
    else:
        _logger.info('[SII] lines missing for S_23 could not be calculated')

    return S_23



def sufur_diaz_2022(lines_log, S2_lines=('S2_6717A', 'S2_6731A'), S3_lines=('S3_9069A', 'S3_9532A'),
                    S2_norm="H1_6563A", S3_norm="H1_6563A", flux_column=f'line_int', n_steps=5000, temp=10000, den=100):


    # Line distribution container
    dist_dict = {}

    # Check the lines are observed and create Monte-Carlo arrays
    for line in list(S2_lines) + list(S3_lines) + [S2_norm] + [S3_norm]:
        int_flux, err_flux = lines_log.loc[line, flux_column], lines_log.loc[line, f'{flux_column}_err']
        dist_dict[line] = np.random.normal(int_flux, err_flux, size=n_steps)

    # Compute the normalization emissivity
    H1 = pn.RecAtom('H', 1)
    line_H_S2 = Line(S2_norm)
    line_H_S3 = Line(S3_norm)

    # Theoretical emissivity
    Hbeta_emis = H1.getEmissivity(temp, den, wave=4861)
    S2_norm_emis = H1.getEmissivity(temp, den, wave=line_H_S2.wavelength[0])
    S3_norm_emis = H1.getEmissivity(temp, den, wave=line_H_S3.wavelength[0])

    # Calculation S_23
    S_2 = (dist_dict[S2_lines[0]] + dist_dict[S2_lines[1]]) / dist_dict[S2_norm]
    S_3 = (dist_dict[S3_lines[0]] + dist_dict[S3_lines[1]]) / dist_dict[S3_norm]
    S23 = S_2 * (S2_norm_emis/Hbeta_emis) + S_3 * (S3_norm_emis/Hbeta_emis)

    # Calculation abundance coefficients
    k1 = np.random.normal(6.636, 0.011, n_steps)
    k2 = np.random.normal(2.202, 0.050, n_steps)
    k3 = np.random.normal(1.060, 0.098, n_steps)

    # Abundance calculation
    S_H = k1 + k2 * np.log10(S23) + k3 * np.power(np.log10(S23), 2)

    return np.mean(S_H), np.std(S_H)



def assign_temperature_diagnostic(line_list, diagnostic, temp_assign = True):

    # Check if the mathod is in the library:
    if diagnostic[0] in METHODS_DICT.keys():
        diag_func = METHODS_DICT[diagnostic[0]]

    # Confirm is a line:
    else:
        try:
            line = lime.Line(diagnostic[0])
            diag_func = diagnostic[0]
        except:
            raise KeyError(f'Input diagnostic {diagnostic[0]} is not recognized by Specsy')

        if np.any(np.isin(diagnostic[0], line_list)):
            diag_func = diagnostic[0]
        else:
            raise KeyError(f'Input diagnostic line {diagnostic[0]} is not available in observation')

    return diag_func


def review_model(emissivity_grid, prior_dict, temp_zones, verbose=True):

    # Check that the lines are present in the dataframe
    if emissivity_grid is not None:
        lines_cand = np.array(list(emissivity_grid.keys()))
        idcs = np.isin(self.lines, lines_cand)
        if not np.all(idcs):
            raise KeyError(f'- Missing lines from emissivity grid database: {self.lines[~idcs]}')

    # Show the priors configuration
    if verbose:
        print(f'\n- Prior configuration')
        for key, value in prior_dict.items():
            print(f'-- {key.split("_prior")[0]}: {value}')

    # High ionization ions
    if verbose:
        print(f'\n- {len(temp_zones)} temperature zones model with:')
        for temp, ions_array in temp_zones.items():
            print(f'-- {temp}: {ions_array}')

    if verbose:
        print(f'\n- Input fluxes: ')
        for i, line in enumerate(self.lines):
            print(f'-- {line} ({self.particles[i]}):'
                  f'flux = {self.fluxes[i]:.4f} +/- {self.errs[i]:.4f} '
                  f'|| err/flux = {100 * self.errs[i] / self.fluxes[i]:.2f} % '
                  f'|| flambda = {self.f_lambda[i]:.3f}')

    return


def package_results(fname, inference_data, inputs=None, prior_dict=None, true_values=None):

    # First save it just in case
    az.to_netcdf(inference_data, fname)

    # Recalibrate the fluxes
    if "calcFluxes_Op" in inference_data.posterior:
        inference_data.posterior['calcFluxes_Op'] = np.power(10, inference_data.posterior['calcFluxes_Op'])

    # Remove the parametrization
    if prior_dict is not None:
        parameter_list = list(inference_data.posterior.data_vars)
        for param in parameter_list:
            if param in prior_dict:

                # Recover the trace and parametrization
                pos_xarr = inference_data.posterior[param]
                reparam0, reparam1 = prior_dict[param][3], prior_dict[param][4]

                if 'logParams_list' in prior_dict:
                    if param not in prior_dict['logParams_list']:
                        pos_xarr = pos_xarr * reparam0 + reparam1
                    else:
                        pos_xarr = np.power(10, pos_xarr * reparam0 + reparam1)
                else:
                    pos_xarr = pos_xarr * reparam0 + reparam1

                # Reset the data
                inference_data.posterior[param] = pos_xarr

    # Store the inputs in a custom group
    if inputs is not None:
        inputs_dict = {'fluxes': xr.DataArray(data=inputs.fluxes, dims=['labels'],
                                              coords={'labels': inputs.lines}, name='fluxes'),
                       'errs': xr.DataArray(data=inputs.errs, dims=['labels'],
                                            coords={'labels': inputs.lines}, name='errs')}
    else:
        inputs_dict = None

    # Add the true values if provided
    if true_values is not None:
        true_values_dict = {'magnitude': xr.DataArray(data=list(true_values.values()), dims=['parameters'],
                                         coords={'parameters': list(true_values.keys())}, name='magnitude')}
    else:
        true_values_dict = None

    # Save to a file
    save_inference_data(fname, inference_data, inputs=inputs_dict, true_values=true_values_dict)

    return


class InputsDirectMethod:

    def __init__(self, input_lines, input_flux, input_err, flambda_arr, ion_arr,
                 temp_id_arr, den_id_arr, tem_eq_arr, den_eq_arr, eq_flux_arr,):

        self.labels  = input_lines
        self.flux_arr   = input_flux
        self.err_arr    = input_err
        self.flambda_arr  = flambda_arr
        self.ion_arr      = ion_arr
        self.temp_id_arr  = temp_id_arr
        self.den_id_arr   = den_id_arr
        self.eq_tem_arr   = tem_eq_arr
        self.eq_den_arr   = den_eq_arr
        self.eq_flux_arr   = eq_flux_arr

    @classmethod
    def from_dataframe(cls, lines_structure):

        # Remove normalization line
        if 'norm_line' in lines_structure.columns:
            idcs = ~lines_structure.index.isin(lines_structure['norm_line'].unique())
        else:
            idcs = np.ones(lines_structure.index.size).astype(bool)

        # Line inputs the data
        input_lines = lines_structure.loc[idcs].index.to_numpy()
        input_flux  = lines_structure.loc[idcs].line_flux.to_numpy()
        input_err   = lines_structure.loc[idcs].line_flux_err.to_numpy()

        # Unpack physical parameters
        flambda_arr = lines_structure.loc[idcs].f_lambda.to_numpy()
        ion_arr     = lines_structure.loc[idcs].particle.to_numpy()

        # Unpack the temp/den structure arrays
        temp_id_arr = lines_structure.loc[idcs].temp.to_numpy()
        den_id_arr  = lines_structure.loc[idcs].den.to_numpy()
        tem_eq_arr  = lines_structure.loc[idcs].eq_temp.to_numpy()
        den_eq_arr  = lines_structure.loc[idcs].eq_den.to_numpy()
        eq_flux_arr  = lines_structure.loc[idcs].eq_flux.to_numpy()

        print(f'Multi-region direct method sampler')
        print(f'- Readying inputs:')

        return cls(input_lines, input_flux, input_err, flambda_arr, ion_arr,
                   temp_id_arr, den_id_arr, tem_eq_arr, den_eq_arr, eq_flux_arr)


class DirectMethod:

    def __init__(self, lines_df, ion_struct):

        # Default prior cfg
        self.prior_cfg = specsy_cfg['direct_method_priors']

        # Default emissivity
        self.emis_interp = None

        # Containers for the emission lines
        self._lines_frame = lines_df
        self.lines_structure = None

        # Load lines data
        self.norm_line = None
        self.ion_struct = ion_struct

        # Data attributes
        self.model = None
        self.trace = None
        self.inputs = None

        return

    def prepare_inputs(self, line_list=None, emissivity_source=None, prior_cfg=None, kinematic_component=0,
                       R_V=3.1, law='G03 LMC', norm_line='H1_4861A', normalize_flux=True, flux_column='profile_flux',
                       review_model=True):

        # Check the lines frame and normalization line
        if self._lines_frame is None:
            raise SpecSyError(f'The object does not have a lines_frame declared')

        self.norm_line = norm_line
        if self.norm_line not in self._lines_frame.index:
            raise SpecSyError(f'The normalization line "{self.norm_line}" is not in the input lines frame')

        # Prepare the emissivity interpolator
        if isinstance(emissivity_source, (str, Path)):
            emis_dataset = load_dataset(emissivity_source)
        else:
            emis_dataset = emissivity_source
        self.emis_interp = compile_bilinear_interp(emis_dataset)

        # Prepare the prior cfg
        self.prior_cfg = self.prior_cfg if prior_cfg is None else prior_cfg

        # Select the lines
        if line_list is None:
            self.lines_structure = self._lines_frame.copy()
        else:
            self.lines_structure = self._lines_frame.loc[self._lines_frame.index.isin(line_list)].copy()

        # Remove the kinematic components
        if kinematic_component == 0:
            idcs_kinem = ~self.lines_structure.index.str.contains('_k-')
        else:
            idcs_kinem = self.lines_structure.index.contains(f'_k-{kinematic_component}')
        self.lines_structure = self.lines_structure.loc[idcs_kinem]

        # Make a copy of the frame and normalize the fluxes
        if normalize_flux:
            self.lines_structure = normalize_fluxes(self.lines_structure, norm_list=norm_line, flux_column=flux_column,
                                                    clear_empty=True)
        else:
            self.lines_structure.insert(3, 'norm_line', self.norm_line)

        # Compute the reddening law
        flambda_arr = flambda_calc(self.lines_structure.wavelength, R_V, law, self.lines_structure.loc['H1_4861A'].wavelength)
        if 'f_lambda' not in self.lines_structure.columns:
            self.lines_structure.insert(4, 'f_lambda', flambda_arr)

        # Map the target lines to the ionization structure
        self.lines_structure = self.ion_struct.map_line_structure(self.lines_structure)

        # Add the equations per ion
        if 'eq_flux' not in self.lines_structure.columns:
            self.lines_structure.insert(9, 'eq_flux', '-')
        for line_label in self.lines_structure.index:
            eq_name = DEFAULT_PARTICLE_EQUATIONS_KEYS.get(self.lines_structure.at[line_label, 'particle'], 'metals')
            self.lines_structure.loc[line_label, 'eq_flux'] = eq_name

        # Check for issues
        if review_model:
            self._review_inputs()

        return

    def _review_inputs(self, return_message=False):

        errors = []
        no_dash = self.lines_structure['region'] != '-'

        # 1) region must not be '-'
        bad = self.lines_structure.index[~no_dash].tolist()
        if bad:
            errors.append(f"The following lines do not have a region defined: {bad}")

        # 2) for non-'-' region rows, 'temp' and 'den' must not be '-'
        sub = self.lines_structure[no_dash]
        for col in ('temp', 'den'):
            bad = sub.index[sub[col] == '-'].tolist()
            if bad:
                errors.append(f"'{col}' is '-' for non-'-' region rows at lines: {bad}")

        # 3) eq_temp / eq_den values must be keys in their respective dicts
        for col, d in (('eq_temp', _TEM_FUNC_DICT), ('eq_den', _DEN_FUNC_DICT)):
            mask = self.lines_structure[col] != '-'
            bad = self.lines_structure.index[mask & ~self.lines_structure[col].isin(d)].tolist()
            if bad:
                errors.append(f"'{col}' has unrecognised values at lines: {bad}")

        # 4) each particle must map to a unique region
        dup = self.lines_structure.groupby('particle')['region'].nunique()
        bad = dup.index[dup > 1].tolist()
        if bad:
            detail = {p: self.lines_structure.loc[self.lines_structure['particle'] == p, 'region'].unique().tolist() for p in bad}
            errors.append(f"'particle' maps to multiple regions: {detail}")

        # 5) Check if we have emissivity data
        bad = [idx for idx in self.lines_structure.index if idx not in self.emis_interp]
        if bad:
            errors.append(f"Missing emissivity data for transitions: {bad}")

        # 6) Check the flux equation is recognized
        eq_flux_names = self.lines_structure.eq_flux.unique()
        bad = [eq_name for eq_name in eq_flux_names if eq_name not in FLUX_EQUATION_DICT]
        if bad:
            errors.append(f"Flux equation not available in database for: {bad}")

        if not return_message:
            if errors:
                msg = "lines_structure validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
                warnings.warn(msg)
                raise ValueError(msg)
            else:
                return None
        else:
            return errors


    def run(self, draws=1000, tune=2000, chains=8, cores=8, target_accept=0.8, nuts_sampler='numpyro', callback=None,
            linear_scale_results=True):

        # Input data
        self.inputs = InputsDirectMethod.from_dataframe(self.lines_structure)

        # Create the model
        print(f'- Compiling model:')
        self.model = direct_method_multi_region(inputs=self.inputs, emis_interp=self.emis_interp, prior_dict=self.prior_cfg,
                                                tem_EQDB=_TEM_FUNC_DICT, den_EQDB=_DEN_FUNC_DICT)

        # Run the model
        print(f'- Launching sampler:')
        self.trace = run_model(self.model,  draws=draws, tune=tune, chains=chains, cores=cores, target_accept=target_accept,
                               nuts_sampler=nuts_sampler, callback=callback)

        # Remove the normalization from the fluxes
        if linear_scale_results:
            self.trace.posterior['theo_flux'] = np.power(10, self.trace.posterior['theo_flux'])
            self.trace.observed_data['likelihood'] = np.power(10, self.trace.observed_data['likelihood'])

        # Remove the log scale for the helium abundaces
        for helium in ['He1', 'He2']:
            if helium in self.inputs.ion_arr:
                self.trace.posterior[helium] = np.power(10, self.trace.posterior[helium])

        return

    def save_line_structure(self, fname):
        lime.save_frame(fname, self.lines_structure)

        return

    def save_trace(self, fname):
        az.to_netcdf(self.trace, fname)

        return


