import numpyro

# Number of cores available
numpyro.set_host_device_count(12)

import jax

# Selecting CPU sampling with JAX
jax.config.update("jax_platform_name", "cpu")

import lime
import specsy as sy
from specsy.plotting.arviz_functions import plot_fitted_fluxes, plot_fitted_params, plot_fitted_pairs


cfg_fname = './synthetic_spectrum_region_v0.toml'
cfg = lime.load_cfg(cfg_fname)

linesdf_fname = '/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/Tutorial/SHOC579_sdss_measurements.txt'
lines_df = lime.load_frame(linesdf_fname)

# Emissivity file
emis_fname = '/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/Tutorial/emissivity_grids_pyneb_1.1.30.nc'
structure_cfg = cfg['default_ionization_structure']

obj = sy.Nebula.from_lines_frame(lines_df, structure_cfg)
line_list = ['H1_4861A', 'H1_4340A', 'H1_6563A',
             'O2_3726A', 'O2_3729A',
             'O3_4363A', 'O3_4959A', 'O3_5007A',
             'N2_6548A', 'N2_6583A',
             'Ar3_7136A', 'Ar3_7751A', 'Ar4_4740A',
             'S2_6716A', 'S2_6731A',
             'S3_6312A', 'S3_9068A', 'S3_9530A']

# Region line structure file
structure_fname = f'./real_spectrum_region_structure_v0.txt'
obj.infer.direct_method.prepare_inputs(emissivity_source=emis_fname, line_list=line_list, prior_cfg=cfg['direct_method_priors'])
obj.infer.direct_method.save_line_structure(structure_fname)

# Run the fitting
obj.infer.direct_method.run()

# Save the fitting trace
trace_fname = f'./real_spectrum_trace_v0.nc'
obj.infer.direct_method.save_trace(trace_fname)

# Plots
plot_fitted_fluxes(trace_fname)
plot_fitted_params(trace_fname)
plot_fitted_pairs(trace_fname, var_names=['temp_high', 'temp_low', 'den_low', 'O2', 'O3', 'S2'])
