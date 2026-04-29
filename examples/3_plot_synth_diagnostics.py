import arviz as az
import lime
from matplotlib import pyplot as plt
from specsy.plotting.plots import plot_traces

# Synthetic region base parameters
cfg_fname = f'./synthetic_spectrum_region_v0.toml'
lines_fname = f'./synthetic_spectrum_lines_region_v1.txt'
trace_fname = f'./synthetic_spectrum_trace_v1.nc'

# Load the data
spec_cfg = lime.load_cfg(cfg_fname)
structure_df = lime.load_frame(lines_fname)
trace = az.from_netcdf(trace_fname)

# Save the files
var_names = ["den_low", "temp_low", "temp_high", "O2", "O3", "S2", "S3", "N2", "Ar3", "Ar4", "Ne3", "cHBeta", ]

# Make the trace plots
# summary = az.summary(trace, var_names=var_names)
# print(summary)

az.plot_pair(trace, var_names=var_names, divergences=True, gridsize=(len(var_names), len(var_names)))
# az.plot_posterior(trace, var_names=var_names)
plt.show()

# plot_traces(trace_fname)