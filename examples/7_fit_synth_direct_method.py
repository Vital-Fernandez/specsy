import specsy as sy
from pathlib import Path
import innate

# Data location
synthConfigPath = Path('./sample_data/synth_conf.toml')
emissivity_cfg = Path('./sample_data/emissivity_grids.nc')
synthLinesLogPath = Path('./sample_data/synth_linesLog.txt')
output_db = Path('./sample_data/synth_fitting.nc')

from specsy.plotting.plots import theme

theme.set_style('dark')

# Configuration
cfg = sy.load_cfg(synthConfigPath)
emiss_data_set = innate.DataSet.from_file(emissivity_cfg)

# Load lines observations
log = sy.load_frame(synthLinesLogPath, flux_type='intg')

# Declare model
dm_twoTemps = sy.models.DirectMethod(emiss_grids=emiss_data_set, R_v=3.4, extinction_law="G03 LMC")

# Declare model
line_list = ['O2_3726A_m','O3_4363A','O3_4959A','O3_5007A','S3_6312A','S3_9068A','S3_9530A','H1_4340A','H1_4861A','H1_6563A']
idcs = log.index.isin(line_list)
dm_twoTemps.fit.frame(log.loc[idcs], './sample_data/', 'synth_fitting', true_values=cfg['true_values'])

# Plot the results
results_address = f'./sample_data/synth_fitting_infer_db.nc'
sy.plot_traces(results_address, f'./sample_data/synth_traces.png')
sy.plot_flux_grid(results_address, f'./sample_data/flux_posteriors.png')
sy.plot_corner_matrix(results_address, f'./sample_data/corner_matrix.png')

print(f'Finished')



# # Load emissivity grids and generate the interpolators
# emiss_db = sy.Innate(emissivity_file, x_space=cfg['simulation_properties']['temp_grid_array'],
#                      y_space=cfg['simulation_properties']['den_grid_array'])
#


