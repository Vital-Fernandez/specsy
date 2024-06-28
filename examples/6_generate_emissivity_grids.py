import numpy as np
import lime
import innate
import specsy as sy
from pathlib import Path


# Data location
synthConfigPath = Path('./sample_data/synth_conf.toml')
synthLinesLogPath = Path('./sample_data/synth_linesLog.txt')
emissivity_conf = Path('./sample_data/emissivity_formulae.toml')
emissivity_file = Path('./sample_data/emissivity_grids.nc')
old_emissivity_file = Path('./sample_data/old_emissivity_grids.nc')
output_db = Path('./sample_data/synth_fitting.nc')

# Load the data
cfg = sy.load_cfg(synthConfigPath)
emis = sy.load_cfg(emissivity_conf)
log = sy.load_frame(synthLinesLogPath, flux_type='intg')

# # Create a database with emissivity grids
# lines_db = lime.line_bands((3000, 10000))
# lines_db['norm_line'] = 'H1_4861A'
# lines_db['group_label'] = 'none'
# lines_db.loc['O2_3726A_m'] = lines_db.loc['O2_3726A']
# lines_db.loc['O2_3726A_m', 'wavelength'] = 3726.0300
# lines_db.loc['O2_3726A_m', 'group_label'] = "O2_3726A+O2_3729A"
# lines_db.loc['O2_7319A_m'] = lines_db.loc['O2_7319A']
# lines_db.loc['O2_7319A_m', 'wavelength'] = 7318.8124
# lines_db.loc['O2_7319A_m', 'group_label'] = "O2_7319A+O2_7330A"
# # emiss_dict = sy.models.emissivity.generate_emis_grid_orig(lines_db, norm_header='norm_line')
# # sy.innate.save_grids(old_emissivity_file, emiss_dict)
#
# # Data grids parameters
# params = 'tem', 'den'
# temp_range = (9000, 20000, 251)
# den_range = (1, 600, 101)
#
# data_conf = {'parameter': 'emissivity',
#              'approximation': ('rgi', 'eqn'),
#              'axes': ('temp', 'den'),
#              'temp_range': (9000, 20000, 251),
#              'den_range': (1, 600, 101)}
#
# # Prepare the configuration for the formula regression
# trans_conf = {}
# for key, value in emis['emissivity_equations'].items():
#     if not key.startswith('emisEquation'):
#         string_equ = emis['emissivity_equations'][emis['emissivity_equations'][key]]
#         trans_conf[key] = {'eqn': string_equ,
#                            'eqn_coeffs': emis['emissivity_coefficients'][key]}
#
# # Prepare data grids:
# grid_dict = {}
# temp_array = np.linspace(data_conf['temp_range'][0], data_conf['temp_range'][1], data_conf['temp_range'][2])
# den_array = np.linspace(data_conf['den_range'][0], data_conf['den_range'][1], data_conf['den_range'][2])
# print('Computing emissivities:')
# for line in lines_db.index:
#     print(f'-- {line}')
#     grid_dict[line] = sy.models.emissivity.generate_emis_grid(line, temp_array, den_array, lines_df=log,
#                                                               normalization_line='H1_4861A', log_scale=True)
#
# # Save the data into a dictionary
# innate.save_dataset(emissivity_file, grid_dict, data_conf, trans_conf)

# Compare old and new version
emiss_old = sy.Innate(old_emissivity_file, x_space=cfg['simulation_properties']['temp_grid_array'],
                                       y_space=cfg['simulation_properties']['den_grid_array'])

emiss_new = innate.DataSet.from_file(Path('./sample_data/emissivity_grids_porSi.nc'))

print('Coso')
print(emiss_old.grid['H1_6563A'])
print(emiss_new['H1_6563A'].data)
# print(np.all(np.isclose(emiss_old.grid['H1_6563A'], emiss_new['H1_6563A'].data)))

# data_conf = {'parameter': 'emissivity', 'approximation': ('rgi', 'eqn'), 'axes': ('temp', 'den'),
#              'temp_range': (9000, 20000, 251), 'den_range': (1, 600, 101)}
# trans_conf = {'H1_6563A': {'eqn': 'y~2.3x + 4.3x**2 + 6.3x**3'}}
#
# # Save the data into a dictionary
# output_file = f'../data/emissivity_grids.nc'
# innate.save_dataset(output_file, emiss_dict, data_conf, trans_conf)

# Generate the containers for the data and configuration

# lines_db['f_lambda'] = 0.0
# for index in lines_db.index:
#     f_lambda_i = flambda_calc(lines_db.loc[index, 'wavelength'], 3.4, "G03 LMC", 4861.2582)
#     lines_db.loc[index, 'f_lambda'] = f_lambda_i
#     print(index, f_lambda_i)
#
# lines_db.sort_values(by='wavelength', inplace=True)

# f_lambda_dict = lines_db.f_lambda.to_dict()
# print(f_lambda_dict)

# # Create a database with emissivity grids
# lines_db = lime.line_bands((3000, 10000))
# lines_db['norm_line'] = 'H1_4861A'
# lines_db['group_label'] = 'none'
# lines_db.loc['O2_3726A_m'] = lines_db.loc['O2_3726A']
# lines_db.loc['O2_3726A_m', 'wavelength'] = 3726.0300
# lines_db.loc['O2_3726A_m', 'group_label'] = "O2_3726A+O2_3729A"
# lines_db.loc['O2_7319A_m'] = lines_db.loc['O2_7319A']
# lines_db.loc['O2_7319A_m', 'wavelength'] = 7318.8124
# lines_db.loc['O2_7319A_m', 'group_label'] = "O2_7319A+O2_7330A"
# # emiss_dict = sy.models.emissivity.generate_emis_grid(lines_db, norm_header='norm_line')
# # save_grids(emissivity_file, emiss_dict)
#
# lines_db['f_lambda'] = 0.0
# for index in lines_db.index:
#     f_lambda_i = flambda_calc(lines_db.loc[index, 'wavelength'], 3.4, "G03 LMC", 4861.2582)
#     lines_db.loc[index, 'f_lambda'] = f_lambda_i
#     print(index, f_lambda_i)
#
# lines_db.sort_values(by='wavelength', inplace=True)
#
# f_lambda_dict = lines_db.f_lambda.to_dict()
# print(f_lambda_dict)

# # Load lines observations
# log = sy.load_frame(synthLinesLogPath, flux_type='intg')
#
# # Load emissivity grids and generate the interpolators
# emiss_db = sy.Innate(emissivity_file, x_space=cfg['simulation_properties']['temp_grid_array'],
#                      y_space=cfg['simulation_properties']['den_grid_array'])
#
# # Declare model
# dm_twoTemps = sy.models.DirectMethod(emiss_grids=emiss_db, R_v=3.4, extinction_law="G03 LMC")
#
# # Declare model
# dm_twoTemps.fit.frame(log, './sample_data/', 'synth_fitting', true_values=cfg['true_values'])
#
# # Plot the results
# results_address = f'./sample_data/synth_fitting_infer_db.nc'
# sy.plots.plot_traces(results_address, f'./sample_data/synth_traces.png')
# sy.plots.plot_flux_grid(results_address, f'./sample_data/flux_posteriors.png')
# sy.plots.plot_corner_matrix(results_address, f'./sample_data/corner_matrix.png')
#
# print(f'Finished')




