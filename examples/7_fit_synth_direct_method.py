import numpy as np
import lime
import specsy as sy
import innate
from pathlib import Path
from specsy.models.extinction import flambda_calc


# Data location
synthConfigPath = Path('./sample_data/synth_conf.toml')
emissivity_cfg = Path('./sample_data/emissivity_grids.nc')
synthLinesLogPath = Path('./sample_data/synth_linesLog.txt')
output_db = Path('./sample_data/synth_fitting.nc')

# Configuration
cfg = sy.load_cfg(synthConfigPath)
emiss_data_set = innate.DataSet.from_file(emissivity_cfg)

# Load lines observations
log = sy.load_frame(synthLinesLogPath, flux_type='intg')

# Declare model
dm_twoTemps = sy.models.DirectMethod(emiss_grids=emiss_data_set, R_v=3.4, extinction_law="G03 LMC")

# Declare model
dm_twoTemps.fit.frame(log, './sample_data/', 'synth_fitting', true_values=cfg['true_values'])

# Plot the results
results_address = f'./sample_data/synth_fitting_infer_db.nc'
sy.plots.plot_traces(results_address, f'./sample_data/synth_traces.png')
sy.plots.plot_flux_grid(results_address, f'./sample_data/flux_posteriors.png')
sy.plots.plot_corner_matrix(results_address, f'./sample_data/corner_matrix.png')

print(f'Finished')



# # Load emissivity grids and generate the interpolators
# emiss_db = sy.Innate(emissivity_file, x_space=cfg['simulation_properties']['temp_grid_array'],
#                      y_space=cfg['simulation_properties']['den_grid_array'])
#


