import arviz as az
import lime
from matplotlib import pyplot as plt
import specsy.plotting.plots as sy_plots
import specsy.plotting.bokeh_functions as sy_bokeh
import specsy as sy
from bokeh.plotting import show as bokeh_show


sy_plots.theme.set_style('dark')

# Synthetic region base parameters
cfg_fname = f'./synthetic_spectrum_region_v0.toml'
lines_fname = f'./synthetic_spectrum_lines_region_v3.txt'
lines_struc_fname = f'./synthetic_spectrum_line_structure_v3.txt'
emis_fname = f'./emissivity_grids_pyneb_1.1.30.nc'
trace_fname = f'./synthetic_spectrum_trace_v4.nc'
plot_trace = f'/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/Tutorial/trace_plot.png'
plot_sm = f'/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/Tutorial/sm_plot.png'
plot_flug_grid = f'/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/Tutorial/sm_plot.png'

# Load the data
spec_cfg = lime.load_cfg(cfg_fname)
structure_df = lime.load_frame(lines_fname)
trace = az.from_netcdf(trace_fname)

# # Trace plot
# fig_cfg = {'width': 200, 'height': 100}
# fig = sy_bokeh.bokeh_trace(trace, in_fig=None, fig_cfg=fig_cfg)
# bokeh_show(fig)

# # Scatter plot matrix
# fig_cfg = {}
# fig = sy_bokeh.bokeh_scatter_matrix(trace, in_fig=None, fig_cfg=fig_cfg)
# bokeh_show(fig)

# Flux grid
fig_cfg = None
fig = sy_bokeh.bokeh_flux_grid(trace, in_fig=None, fig_cfg=fig_cfg)
bokeh_show(fig)

az.summary(az.summary(trace, stat_focus="median"))

# sy_plots.az_scatter_matrix(trace, plot_sm)
# sy_plots.az_flux_grid(trace, plot_flug_grid, structure_df, n_cols=4)
