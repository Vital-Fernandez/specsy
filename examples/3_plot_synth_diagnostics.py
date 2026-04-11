import arviz as az
import lime
from matplotlib import pyplot as plt

# Synthetic region base parameters
cfg_fname = f'synthetic_region_v0.toml'
lines_fname = f'./synthetic_lines_region_v0.txt'
trace_fname = f'./synthetic_trace_v0.nc'

# Load the data
spec_cfg = lime.load_cfg(cfg_fname)
structure_df = lime.load_frame(lines_fname)
trace = az.from_netcdf(trace_fname)

# Make the trace plots
var_names = ["O2", "O3", "S2", "S3", "N2", "Ar3", "Ar4", "Ne3", "cHBeta", "den_low", "temp_low", "temp_high"]
az.plot_pair(trace, var_names=var_names, divergences=True)
az.plot_posterior(trace, var_names=var_names)
summary = az.summary(trace, var_names=var_names)
print(summary)
plt.show()

