import numpyro

# Number of cores available
numpyro.set_host_device_count(8)

import jax

# Selecting CPU sampling with JAX
jax.config.update("jax_platform_name", "cpu")

import specsy as sy

# Synthetic region base parameters
cfg_fname = f'./synthetic_spectrum_region_v0.toml'
lines_fname = f'./synthetic_spectrum_lines_region_v5.txt'
emis_fname = f'./emissivity_grids_pyneb_1.1.30.nc'
trace_fname = f'./synthetic_spectrum_trace_v6.nc'
lines_struc_fname = f'./synthetic_spectrum_line_structure_v6.txt'

# Load the data
spec_cfg = sy.load_cfg(cfg_fname)
lines_df = sy.load_frame(lines_fname)

obj = sy.Nebula.from_lines_frame(lines_df, spec_cfg['default_ionization_structure'])

obj.infer.direct_method.prepare_inputs(emissivity_source=emis_fname, normalize_flux=False, prior_cfg=spec_cfg['direct_method_priors'])

print(obj.infer.direct_method.lines_structure.to_string())
obj.infer.direct_method.save_line_structure(lines_struc_fname)

obj.infer.direct_method.run(nuts_sampler='pymc')

obj.infer.direct_method.save_trace(trace_fname)

# obj.infer.direct_method.plot_trace()


# # Create the chemical model
# chem_model = sy.DirectMethod(structure_df, emis_fname, spec_cfg['direct_method_priors'])
#
# # Run the sampler
# chem_model.run()
#
# # Save the results
# # chem_model.save_trace(trace_fname)
#
# # Show the results
# chem_model.plot_trace()
#

# import pymc as pm
# from pymc.progress_bar import ProgressBarManager
#
# old_update = ProgressBarManager.update
#
#
# def new_update(self, chain_idx, is_last, draw, tuning, stats):
#     print(chain_idx, draw)  # Do whatever you want with this info
#     old_update(self, chain_idx, is_last, draw, tuning, stats)
#
#
# ProgressBarManager.update = new_update
#
# with pm.Model() as m:
#     x = pm.Normal("x")
#     pm.sample(tune=5, draws=5, chains=2)
