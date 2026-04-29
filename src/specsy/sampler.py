import pymc as pm
import numpy as np
from pytensor import tensor as tt


def set_prior(param, prior_dict, abund_type=False, name_param=None):

    # Read distribution configuration
    dist_name = prior_dict[param][0]
    dist_loc, dist_scale = prior_dict[param][1], prior_dict[param][2]
    dist_norm, dist_reLoc = prior_dict[param][3], prior_dict[param][4]

    # Load the corresponding probability distribution
    probDist = getattr(pm, dist_name)

    if abund_type:
        priorFunc = probDist(name_param, dist_loc, dist_scale) * dist_norm + dist_reLoc

    elif probDist.__name__ in ['HalfNormal']:  # These distributions only have one parameter
        priorFunc = probDist(param, dist_scale) * dist_norm + dist_reLoc

    elif probDist.__name__ in ['HalfCauchy']:  # These distributions only have one parameter
        priorFunc = probDist(param, dist_scale) * dist_norm + dist_reLoc

    elif probDist.__name__ == 'Uniform':

        if param == 'logOH':
            priorFunc = pm.Bound(pm.Normal, lower=7.1, upper=9.1)('logOH', mu=8.0, sigma=1.0, testval=8.1)
        if param == 'logU':
            priorFunc = pm.Bound(pm.Normal, lower=-4.0, upper=-1.5)('logU', mu=-2.75, sigma=1.5, testval=-2.75)
        if param == 'logNO':
            priorFunc = pm.Bound(pm.Normal, lower=-2.0, upper=0.0)('logNO', mu=-1.0, sigma=0.5, testval=-1.0)

    else:
        priorFunc = probDist(param, dist_loc, dist_scale) * dist_norm + dist_reLoc

    return priorFunc


def direct_method_multi_region(inputs, emis_interp, prior_dict, tem_EQDB, den_EQDB):

    # Observables
    labels_arr = inputs.labels
    input_obs = np.log10(inputs.flux_arr)
    input_err = np.log10(1 + inputs.err_arr / inputs.flux_arr)
    ion_arr = inputs.ion_arr
    flambda_arr = inputs.flambda_arr

    # Structure array
    temp_id_arr = inputs.temp_id_arr
    den_id_arr = inputs.den_id_arr
    eq_tem_arr = inputs.eq_tem_arr
    eq_den_arr = inputs.eq_den_arr

    # Convenience arrays
    range_arr = np.arange(labels_arr.size)
    temp_eq_check = eq_tem_arr == '-'
    den_eq_check = eq_den_arr == '-'
    unique_species = np.unique(ion_arr)
    unique_params = np.unique(np.concatenate([temp_id_arr, den_id_arr]))

    # PyMC model
    with (pm.Model(coords={"lines": labels_arr}) as model):

        # Container to store the models
        theo_flux = tt.zeros(labels_arr.size)

        # Compile the abundances
        for ion in unique_species:
            if ion != 'H1':
                set_prior(ion, prior_dict, abund_type=True, name_param=ion)

        # Extinction
        # cHbeta = pm.HalfCauchy(name="cHBeta", beta=2)
        cHbeta = set_prior('cHBeta', prior_dict)

        # Generate the free temperatures and densities
        for param in unique_params:
            set_prior(param, prior_dict)

        # Loop through the lines and compute the fluxes
        for i in range_arr:

            # Compute the emissivity
            tem = model[temp_id_arr[i]] if temp_eq_check[i] else tem_EQDB[eq_tem_arr[i]](model[temp_id_arr[i]])
            den = model[den_id_arr[i]] if den_eq_check[i] else den_EQDB[eq_den_arr[i]](model[den_id_arr[i]])
            emis = emis_interp[labels_arr[i]](tem, den)

            if ion_arr[i] == 'H1':
                flux = emis - flambda_arr[i] * cHbeta
            else:
                flux = model[ion_arr[i]] + emis - flambda_arr[i] * cHbeta - 12

            theo_flux = tt.inc_subtensor(theo_flux[i], flux)

        # Stored the fluxes
        pm.Deterministic('theo_flux', theo_flux)

        # Likelihood
        pm.Normal("likelihood", mu=theo_flux, sigma=input_err, observed=input_obs, dims='lines')

    return model

def run_model(model, draws=1000, tune=2000, target_accept=0.9, chains=8, cores=8, nuts_sampler='numpyro'):

    with model:
        trace = pm.sample(draws=draws, tune=tune, target_accept=target_accept, chains=chains, cores=cores, nuts_sampler=nuts_sampler)

    return trace