import logging
import numpy as np
import lime
import arviz as az
import xarray as xr

from pathlib import Path
from pandas import DataFrame
from lime import label_decomposition

from .. import _setup_cfg
from ..operations.pytensors import EmissionFluxModel
from ..tools import truncated_gaussian, flux_distribution
from .chemistry_inference import direct_method_inference
from .extinction import flambda_calc
from ..innate import save_inference_data

try:
    import pyneb as pn
    pyneb_check = True
except ImportError:
    pyneb_check = False

_logger = logging.getLogger('SpecSy')


def TOIII_from_TSIII_relation(T_low):
    # From Hagele et al 2006
    return (0.8403 * T_low / 10000.0 + 0.2689) * 10000.0


def TOII_from_TOIII_relation(T_high, n_e):
    # From Epm and Cotini 2009
    return ((1.2 + 0.002*n_e + 4.2/n_e) / (10000.0/T_high + 0.08 + 0.003*n_e + 2.5/n_e)) * 10000.0


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


def sulfur_diaz_2020(S_23):

    n_steps = S_23.size
    a_dist = np.random.normal(6.636, 0.010, size=n_steps)
    b_dist = np.random.normal(2.202, 0.050, size=n_steps)
    c_dist = np.random.normal(1.060, 0.098, size=n_steps)

    SH = a_dist + b_dist * np.log10(S_23) + c_dist * np.square(np.log10(S_23))

    return SH


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
    line_H_S2 = lime.Line(S2_norm)
    line_H_S3 = lime.Line(S3_norm)

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



class DmInputs:

    def __init__(self, lines_frame, lines_list=None, ext_frame='FRAME', R_v=None, extinction_law=None):

        # Attributes
        self.frame = None
        self.lines = None
        self.fluxes = None
        self.errs = None
        self.f_lambda = None
        self.wave = None
        self.particles = None

        # Review the inputs
        self.frame = lime.io.check_file_dataframe(lines_frame, DataFrame, ext=ext_frame)

        # Crop to the target lines
        if lines_list is not None:
            self.frame = self.frame.index.isin(lines_list)

        # Declare the inputs
        self.lines = self.frame.index.to_numpy()
        self.particles = self.frame.particle.to_numpy()
        self.fluxes = self.frame.line_flux.to_numpy()
        self.errs = self.frame.line_flux_err.to_numpy()
        self.wave = self.frame.wavelength.to_numpy()

        # Compute the extinction # TODO normalization for multiple lines
        norm_line = np.unique(self.frame.norm_line.to_numpy())[0]
        line = lime.Line(norm_line)
        self.f_lambda = flambda_calc(self.wave, R_v, extinction_law, line.wavelength)

        return

    def review_model(self, emissivity_grid, prior_dict, temp_zones, verbose=True):

        # Check that the lines are present in the dataframe
        if emissivity_grid is not None:
            lines_cand = np.array(list(emissivity_grid.interpl.keys()))
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

class DmFunctions():

    def __init__(self, model):

        # Instantiate the dependencies
        # LineFitting.__init__(self)

        # Lime spectrum object with the scientific data
        self._model = model
        self._i_iter = 0
        self._n_iter = 0

        self.output_path = None
        self.label_fit = None

        return

    def frame(self, lines_frame, output_folder, results_label, lines_list=None, iter=2000, cores=4,
              plots=('traces', 'posterior', 'sc_matrix'), verbose=True, true_values=None):

        # Check output path
        self.output_path, self.label_fit = Path(output_folder), results_label
        assert self.output_path.is_dir(), f'- Directory "{output_folder}" does not exist'

        # Prepare the input data for the fitting
        inputs = DmInputs(lines_frame, lines_list, R_v=self._model.R_v, extinction_law=self._model.extinction_law)
        inputs.review_model(self._model.emiss_grids, self._model.prior_conf, self._model.temp_zones, verbose)

        # Prepare auxiliary parameters
        idcs_highTemp_ions = np.isin(inputs.particles, self._model.temp_zones['high'])
        lowTemp_check = np.any(np.isin(self._model.temp_low_diag, inputs.lines))
        highTemp_check = np.any(np.isin(self._model.temp_high_diag, inputs.lines))

        # Confirm all data is available
        if lowTemp_check or highTemp_check:

            # Output file
            fname = self.output_path/f'{self.label_fit}_inference_data.nc'
            print(f'\n- Launching direct method inference: ')

            # Run the model
            infer_data = direct_method_inference(fname, inputs,
                                    prior_dict=self._model.prior_conf, idcs_highTemp_ions=idcs_highTemp_ions,
                                    emiss_interp=self._model.emiss_grids.interpl, eq_tt=self._model.eq_tt,
                                    fit_Tlow=lowTemp_check, fit_Thigh=highTemp_check)

            # Save the output data
            output_db = self.output_path / f'{self.label_fit}_infer_db.nc'
            print(f'-- done, saving the results at: {output_db}')
            self.package_results(output_db, infer_data, inputs, self._model.prior_conf, true_values=true_values)

        else:
            msg = (f'\n- The observation does not have temperature diagnostics:'
                   f'\n-- Low temperature diagnostics: {self._model.temp_low_diag}',
                   f'\n-- High temperature diagnostics: {self._model.temp_high_diag}')
            print(msg)

        return

    def package_results(self, fname, inference_data, inputs=None, prior_dict=None, true_values=None):

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


class DirectMethod:

    def __init__(self, emiss_grids, R_v, extinction_law, min_err=0.02, temp_low_diag=None, temp_high_diag=None,
                 temp_zones=None, den_zones=None, prior_cfg=None, tensor_model='pytensor'):

        # Declare the attributes
        self.emiss_grids = None
        self.candidate_lines = None
        self.R_v = None
        self.extinction_law = None
        self.tensor_library = None
        self.eq_tt = None
        self.prior_conf = None

        # Extinction parameters
        self.R_v = R_v
        self.extinction_law = extinction_law

        # Grid dictionary
        self.emiss_grids = emiss_grids

        # Compute flux equations as tensors
        self.tensor_library = tensor_model
        line_array = np.array(list(self.emiss_grids.interpl.keys()))
        particle_array = lime.label_decomposition(line_array, params_list=['particle'])[0]
        self.eq_tt = EmissionFluxModel(line_array, particle_array)

        # Get the prior configuration
        self.prior_conf = prior_cfg if prior_cfg is not None else _setup_cfg['direct_method_priors']

        # If temps zones is None
        self.temp_zones = temp_zones if temp_zones is not None else _setup_cfg["direct_method_cfg"]['temp_zones']

        # If temps zones is None
        self.temp_low_diag = temp_low_diag if temp_low_diag is not None else _setup_cfg["direct_method_cfg"]['temp_low_diag']
        self.temp_high_diag = temp_high_diag if temp_high_diag is not None else _setup_cfg["direct_method_cfg"]['temp_high_diag']
        self.temp_low_diag, self.temp_high_diag = np.atleast_1d(self.temp_low_diag), np.atleast_1d(self.temp_high_diag)

        # Declare the function methods
        self.fit = DmFunctions(self)

        return