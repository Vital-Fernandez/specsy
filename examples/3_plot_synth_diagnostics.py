import arviz as az
import specsy as sy
import specsy.plotting.arviz_functions as sy_plots
from specsy.plotting.arviz_functions import plot_fitted_fluxes, plot_fitted_params, plot_fitted_pairs, plot_prior_posterior, plot_traces


# Theme settings
sy_plots.theme.set_style('dark')

# Inputs
cfg_fname = f'./synthetic_spectrum_region_v0.toml'
trace_fname = f'./synthetic_spectrum_trace_v7.nc'
true_values = sy.load_cfg(cfg_fname)['synth_spectrum']['true_params']

# Load the fitting data
trace_data = az.from_netcdf(trace_fname)

# Plots
plot_fitted_fluxes(trace_fname)

plot_fitted_params(trace_fname, true_values=true_values)

plot_fitted_pairs(trace_fname, var_names=['temp_high', 'temp_low', 'den_low', 'O2', 'O3', 'S2'])

plot_prior_posterior(trace_fname)

plot_traces(trace_fname)

