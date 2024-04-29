import numpy as np
import lime
import specsy as sy
from specsy.innate import save_grids
from pathlib import Path


# Data location
synthConfigPath = Path('./sample_data/synth_conf.toml')
synthLinesLogPath = Path('./sample_data/synth_linesLog.txt')
emissivity_file = Path('./sample_data/emissivity_db.nc')
output_db = Path('./sample_data/synth_fitting.nc')

# Configuration
cfg = sy.load_cfg(synthConfigPath)

# Create a database with emissivity grids
# lines_db = lime.line_bands((3000, 10000))
# lines_db['norm_line'] = 'H1_4861A'
# lines_db['group_label'] = 'none'
# lines_db.loc['O2_3726A_m'] = lines_db.loc['O2_3726A']
# lines_db.loc['O2_3726A_m', 'wavelength'] = 3726.0300
# lines_db.loc['O2_3726A_m', 'group_label'] = "O2_3726A+O2_3729A"
# lines_db.loc['O2_7319A_m'] = lines_db.loc['O2_7319A']
# lines_db.loc['O2_7319A_m', 'wavelength'] = 7318.8124
# lines_db.loc['O2_7319A_m', 'group_label'] = "O2_7319A+O2_7330A"
# emiss_dict = sy.models.emissivity.generate_emis_grid(lines_db, norm_header='norm_line')
# save_grids(emissivity_file, emiss_dict)

# Load lines observations
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

# Plot the results
results_address = f'./sample_data/synth_fitting_infer_db.nc'
sy.plots.plot_traces(results_address, f'./sample_data/synth_traces.png')
sy.plots.plot_flux_grid(results_address, f'./sample_data/flux_posteriors.png')
sy.plots.plot_corner_matrix(results_address, f'./sample_data/corner_matrix.png')

print(f'Finished')




