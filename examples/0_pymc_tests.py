import pymc as pm
import arviz as az
from matplotlib import pyplot as plt
from pathlib import Path
from fastprogress import fastprogress
fastprogress.printing = lambda: True
import specsy
import xarray as xr

# Assume 10 trials and 5 successes out of those trials
# Change these numbers to see how the posterior plot changes
trials = 10
successes = 5


# Locate the array
conf_file = Path('./sample_data/synth_conf.toml')
synth_cfg = specsy.load_cfg(conf_file)
true_values = synth_cfg['true_values']


# Store the inputs and true values in a custom group
lineLabels = ['H1_4861', 'H1_6563A', 'O3_4959A',  'O3_5007A']
emissionFluxes = [2, 6, 4, 12]
emissionErr = [0.1, 0.2, 0.3, 0.4]

inputs_dict = {'fluxes': xr.DataArray(data=emissionFluxes, dims=['labels'],
                                      coords={'labels': lineLabels}, name='fluxes'),
               'errs': xr.DataArray(data=emissionErr, dims=['labels'],
                                    coords={'labels': lineLabels}, name='errs')}

inputs_dict = None

true_values = {'magnitude': xr.DataArray(data=list(true_values.values()), dims=['parameters'],
                                          coords={'parameters': list(true_values.keys())},
                                          name='magnitude')}
true_values = None

print(pm.__version__)
print(az.__version__)

# Set up models context
with pm.Model() as coin_flip_model:
    # Probability p of success we want to estimate
    # and assign Beta prior
    p = pm.Beta("p", alpha=1, beta=1)

    # Define likelihood
    obs = pm.Binomial("obs", p=p, n=trials, observed=successes)

    # Hit Inference Button
    trace = pm.sample(cores=4, chains=4)

az.to_netcdf(trace, 'example1.nc')

specsy.save_inference_data('example1_specsy.nc', trace)#, inputs_specsy=inputs_dict, true_values=true_values)
trace2 = specsy.load_inference_data('example1_specsy.nc')



lines = trace2.inputs_specsy.labels.values
fluxes = trace2.inputs_specsy.fluxes.values
errs = trace2.inputs_specsy.errs.values
true_params = trace2.true_values.parameters.values
true_mag = trace2.true_values.magnitude.values