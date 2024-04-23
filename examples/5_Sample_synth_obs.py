import numpy as np
import specsy as sy
from pathlib import Path
from specsy.astro.extinction import flambda_calc
from specsy.operations.interpolation import emissivity_grid_calc

from fastprogress import fastprogress
fastprogress.printing = lambda: True

# Data location
synthConfigPath = Path('./sample_data/synth_conf.toml')
synthLinesLogPath = Path('./sample_data/synth_linesLog.txt')
output_db = Path('./sample_data')

# Load simulation parameters
model_cfg = sy.load_cfg(synthConfigPath)

# Load emission lines
input_lines = model_cfg['inference_model_configuration']['input_lines_list']
merged_lines = {'O2_3726A_m': 'O2_3726A+O2_3729A', 'O2_7319A_m': 'O2_7319A+O2_7330A'}
log = sy.load_frame(synthLinesLogPath)

normLine = 'H1_4861A'
idcs_lines = (log.index != normLine)
lineLabels = log.loc[idcs_lines].index
lineWaves = log.loc[idcs_lines, 'wavelength'].values
lineIons = log.loc[idcs_lines, 'ion'].values
lineFluxes = log.loc[idcs_lines, 'intg_flux'].values
lineErr = log.loc[idcs_lines, 'intg_err'].values

# Compute the reddening curve for the input emission lines
flambda = flambda_calc(lineWaves, red_curve="G03 LMC", R_v=3.4, norm_wavelength=4861)

# Interpolator functions for the emissivity grids
emis_grid_interp = emissivity_grid_calc(lines_array=lineLabels, comp_dict=merged_lines)

# Declare sampler
obj1_model = sy.SpectraSynthesizer(emis_grid_interp)

# Declare region physical model
obj1_model.define_region(lineLabels, lineFluxes, lineErr, flambda, merged_lines)

# Declare sampling properties
obj1_model.simulation_configuration(prior_conf_dict=model_cfg['priors_configuration'],
                                    highTempIons=model_cfg['simulation_properties']['high_temp_ions_list'],)

# Declare simulation inference model
obj1_model.inference_model()

obj1_model.save_fit(output_db, ext_name='synth_sampling')


# Run the simulation
# obj1_model.run_sampler(2000, 2000, nchains=4, njobs=1)

# # Load the results
# fit_pickle = sr.load_fit_results(output_db)
# inLines, inParameters = fit_pickle['inputs']['lines_list'], fit_pickle['inputs']['parameter_list']
# inFlux, inErr = fit_pickle['inputs']['line_fluxes'].astype(float), fit_pickle['inputs']['line_err'].astype(float)
# traces_dict = fit_pickle['outputs']
#
# # Print the results
# print('-- Model parameters table')
# figure_file = user_folder/f'obj_fitted_fluxes'
# sr.table_fluxes(figure_file, inLines, inFlux, inErr, traces_dict, merged_lines)
#
# # Print the results
# print('-- Fitted fluxes table')
# figure_file = user_folder/f'obj_MeanOutputs'
# sr.table_params(figure_file, inParameters, traces_dict, true_values=objParams['true_values'])
#
# print('-- Model parameters posterior diagram')
# figure_file = user_folder/f'obj_traces_plot.png'
# sr.plot_traces(figure_file, inParameters, traces_dict, true_values=objParams['true_values'])
#
# print('-- Line flux posteriors')
# figure_file = user_folder/f'obj_fluxes_grid.png'
# sr.plot_flux_grid(figure_file, inLines, inFlux, inErr, traces_dict)
#
# print('-- Model parameters corner diagram')
# figure_file = user_folder/f'obj_corner.png'
# sr.plot_corner(figure_file, inParameters, traces_dict, true_values=objParams['true_values'])



