from pathlib import Path
import specsy
from specsy.plotting.plots import plot_corner_matrix

specsy.theme.set_style('dark')

# # Locate the array
conf_file = Path('./sample_data/synth_conf.toml')
output_db = Path('./sample_data/synth_fitting_true.nc')

synth_cfg = specsy.load_cfg(conf_file)
inference_data = specsy.load_inference_data(output_db)

plot_corner_matrix(output_db, params_list=['T_low', 'T_high', 'S2', 'S3', 'n_e'])

# # Recalibrate the fluxes
# if "calcFluxes_Op" in trace.posterior:
#     trace.posterior['calcFluxes_Op'] = np.power(10, trace.posterior['calcFluxes_Op'])
#
# prior_dict = {}
#
# # Loop through the parameters
# parameter_list = list(trace.posterior.data_vars)
# for param in parameter_list:
#     if param in prior_dict:
#
#         # Recover the trace and parametrization
#         pos_xarr = trace.posterior[param]
#         reparam0, reparam1 = prior_dict[param][3], prior_dict[param][4]
#
#         if 'logParams_list' in prior_dict:
#             if param not in prior_dict['logParams_list']:
#                 pos_xarr = pos_xarr * reparam0 + reparam1
#             else:
#                 pos_xarr = np.power(10, pos_xarr * reparam0 + reparam1)
#         else:
#             pos_xarr = pos_xarr * reparam0 + reparam1
#
#         # Set the umparametrization
#         trace.posterior[param] = pos_xarr

# # Save inputs
# line_list = np.array(['H1_4861A', 'H1_6563A', 'O3_4959A', 'O3_5007A'])
# flux_list = np.array([2., 6., 4, 12.])
# err_list = np.array([0.05, 0.06, 0.04, 0.07])
#
# input_dataset = xr.Dataset({'fluxes': xr.DataArray(data=flux_list, dims=['labels'],
#                                                    coords={'labels': line_list}, name='fluxes'),
#                             'errs': xr.DataArray(data=err_list, dims=['labels'],
#                                                  coords={'labels': line_list}, name='errs')})
#
# with catch_warnings():
#     simplefilter("ignore", UserWarning)
#     trace.add_groups({"specsy_inputs": input_dataset})
#
# # Restore the data
# parameter_list = list(trace.posterior.data_vars)
# line_list = trace.specsy_inputs.labels.values
# line_fluxes = trace.specsy_inputs.fluxes.values
# line_errs = trace.specsy_inputs.errs.values
# trace_arrs = trace.posterior["calcFluxes_Op"].values
