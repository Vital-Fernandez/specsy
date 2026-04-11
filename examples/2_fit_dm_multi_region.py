import numpyro
numpyro.set_host_device_count(8)

import jax
jax.config.update("jax_platform_name", "cpu") # Force JAX to use CPU

import lime
from innate import load_dataset
from specsy.operations.interpolation import compile_bilinear_interp
from specsy.models.chemistry_inference import direct_method_multi_region


# Synthetic region base parameters
cfg_fname = f'synthetic_region_v0.toml'
lines_fname = f'./synthetic_lines_region_v0.txt'
trace_fname = f'./synthetic_trace_v0.nc'

spec_cfg = lime.load_cfg(cfg_fname)
structure_df = lime.load_frame(lines_fname)

# Load the emissivity interpolator
emis_file = load_dataset('./emissivity_grids_pyneb_1.1.30.nc')
interp_dict_np = compile_bilinear_interp(emis_file)

# Run the Bayesian sampler
true_params = spec_cfg['synth_spectrum']['true_params']
priors = spec_cfg['direct_method_priors']
trace = direct_method_multi_region(structure_df, interp_dict_np, priors, fname=trace_fname)


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
