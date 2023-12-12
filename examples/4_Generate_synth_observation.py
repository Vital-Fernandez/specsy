import numpy as np
import pandas as pd
import specsy as sy

from specsy.astro.chemistry import TOIII_from_TSIII_relation
from specsy.astro.extinction import flambda_calc
from specsy.operations.tensors import EmissionTensors
from specsy.operations.interpolation import emissivity_grid_calc
import pyneb as pn

# Load the model configuration
model_cfg = sy.load_cfg('./sample_data/default_cfg.toml')

# Set the synthetic observation parameter values
param_dict = {'n_e': 150.0,
              'T_low': 12500.0,
              'T_high': TOIII_from_TSIII_relation(12500.0),
              'tau': 0.60 + 0.15,
              'cHbeta': 0.08 + 0.02,
              'H1': 0.0,
              'He1': np.log10(0.070 + 0.005),
              'He2': np.log10(0.00088 + 0.0002),
              'O2': 7.805 + 0.15,
              'O3': 8.055 + 0.15,
              'N2': 5.845 + 0.15,
              'S2': 5.485 + 0.15,
              'S3': 6.94 + 0.15,
              'Ne3': 7.065 + 0.15,
              'Fe3': 5.055 + 0.15,
              'Ar3': 5.725 + 0.15,
              'Ar4': 5.065 + 0.15}

# Declare lines to simulate
input_lines = ['H1_4341A', 'H1_4861A', 'H1_6563A',
                'He1_4026A', 'He1_4471A', 'He1_5876A', 'He1_6678A', 'He1_7065A', 'He2_4686A',
                'O2_3726A_m', 'O2_7319A_m', 'O3_4363A', 'O3_4959A', 'O3_5007A',
                'Ne3_3968A',
                'N2_6548A', 'N2_6584A',
                'S2_6716A', 'S2_6731A', 'S3_6312A', 'S3_9069A', 'S3_9531A',
                'Ar3_7136A', 'Ar4_4740A',
                'Fe3_4658A']

merged_lines = {'O2_3726A_m': 'O2_3726A+O2_3729A',
                'O2_7319A_m': 'O2_7319A+O2_7330A'}

# Generate dataframe with store line data
particle_array, wave_array, latex_array = sy.label_decomposition(input_lines, fit_conf=merged_lines)
lineLogHeaders = ['wavelength', 'intg_flux', 'intg_err', 'ion', 'blended_label']

line_log = pd.DataFrame(index=input_lines, columns=lineLogHeaders)
line_log = line_log.assign(wavelength=wave_array, ion=particle_array, blended_label='None')
for line, components in merged_lines.items(): line_log.loc[line, 'blended_label'] = components
line_log.sort_values(by=['wavelength'], ascending=True, inplace=True)

# New line parameter arrays sorted by wavelenght
line_array = line_log.index.to_numpy()
wavelength_array = line_log.wavelength.to_numpy()
particle_array = line_log.ion.to_numpy()

# Compute the reddening curve for the input emission lines
f_lambda_array = flambda_calc(wavelength_array, red_curve="G03 LMC", R_v=3.4, norm_wavelength=4861)

# Interpolator functions for the emissivity grids
emis_grid_interp = emissivity_grid_calc(lines_array=line_array, comp_dict=merged_lines)

# We generate an object with the tensor emission functions
emtt = EmissionTensors(line_array, particle_array)

# High ionization region ions
T_High_ions = model_cfg['simulation_properties']['high_temp_ions_list']

# Loop through the lines and compute the fluxes
flux_array = np.empty(line_log.index.size)
for i in range(len(wavelength_array)):

    # Declare line properties
    lineLabel = line_log.iloc[i].name
    lineIon = line_log.iloc[i].ion
    lineFlambda = f_lambda_array[i]

    # Declare grid interpolation parameters
    emisCoord_low = np.stack(([param_dict['T_low']], [param_dict['n_e']]), axis=-1)
    emisCoord_high = np.stack(([param_dict['T_high']], [param_dict['n_e']]), axis=-1)

    # Compute emisivity for the corresponding ion temperature
    T_calc = emisCoord_high if lineIon in T_High_ions else emisCoord_low
    line_emis = emis_grid_interp[lineLabel](T_calc).eval()

    # Compute line flux
    flux_array[i] = emtt.compute_flux(lineLabel,
                                      line_emis[0][0],
                                      param_dict['cHbeta'],
                                      lineFlambda,
                                      param_dict[lineIon],
                                      ftau=0.0,
                                      O3=param_dict['O3'],
                                      T_high=param_dict['T_high'])

    print(f'{i}) {lineLabel}, f_lambda = {f_lambda_array[i]:0.3f}, F_i = {flux_array[i]:0.3f}')

# Convert to a natural scale
lineFluxes = np.power(10, flux_array)
line_log['intg_flux'] = lineFluxes
line_log['intg_err'] = lineFluxes * model_cfg['simulation_properties']['lines_minimum_error']

# # We proceed to safe the synthetic spectrum as if it were a real observation
log_file_address = './sample_data/synth_linesLog.txt'
sy.save_log(line_log, log_file_address)

# Save the parameters in the natural scale
model_cfg['true_values'] = param_dict
for param in model_cfg['priors_configuration']['logParams_list']:
    param_true_ref = param + '_true'
    model_cfg['true_values'][param] = float(np.power(10, model_cfg['true_values'][param]))
model_cfg['inference_model_configuration']['input_lines_list'] = input_lines

cfg_file_address = './sample_data/synth_conf.toml'
sy.save_cfg(model_cfg, cfg_file_address)

