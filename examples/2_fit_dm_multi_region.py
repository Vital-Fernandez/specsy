import numpyro

# Number of cores available
numpyro.set_host_device_count(12)

import jax

# Selecting CPU sampling with JAX
jax.config.update("jax_platform_name", "cpu")

import specsy as sy

# Synthetic region base parameters
cfg_fname = f'./synthetic_spectrum_region_v0.toml'
lines_fname = f'./synthetic_spectrum_lines_region_v7.txt'
trace_fname = f'./synthetic_spectrum_trace_v7.nc'
lines_struc_fname = f'./synthetic_spectrum_line_structure_v7.txt'

# Load the data
spec_cfg = sy.load_cfg(cfg_fname)
lines_df = sy.load_frame(lines_fname)

# Generate model and fit it
obj = sy.Nebula.from_lines_frame(lines_df, spec_cfg['default_ionization_structure'])
obj.infer.direct_method.prepare_inputs(normalize_flux=False, prior_cfg=spec_cfg['direct_method_priors'])
obj.infer.direct_method.run()

# Save inputs and outputs
obj.infer.direct_method.save_line_structure(lines_struc_fname)
obj.infer.direct_method.save_trace(trace_fname)
