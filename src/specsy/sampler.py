import numpy as np
from time import time
from specsy.models.fluxes_line import FLUX_EQUATION_DICT

try:
    from pytensor import tensor as tt
    import pymc as pm
    pytensor_check = True
except ImportError:
    pytensor_check = False


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
        priorFunc = probDist(param, dist_loc, dist_scale) * dist_norm + dist_reLoc

    else:
        priorFunc = probDist(param, dist_loc, dist_scale) * dist_norm + dist_reLoc

    return priorFunc


def direct_method_multi_region(inputs, emis_interp, prior_dict, tem_EQDB, den_EQDB):

    # This is the good one

    # Convenience arrays
    merge_d = inputs.merge_dict

    # PyMC model
    with (pm.Model(coords={"lines": inputs.labels}) as model):

        # Save input fluxes
        pm.Data('input_flux', inputs.flux_arr, dims='lines')
        pm.Data('input_err', inputs.err_arr, dims='lines')

        # Container to store the models
        theo_flux = tt.zeros(inputs.labels.size)

        # Compile the abundances
        for ion in inputs.unique_species:
            if ion != 'H1':
                set_prior(ion, prior_dict, abund_type=True, name_param=ion)
        pm.Data('H1', 1.0)

        # Extinction
        cHbeta = set_prior('cHBeta', prior_dict)

        # Generate the free temperatures and densities
        for param in inputs.unique_params:
            set_prior(param, prior_dict)

        # Loop through the lines and compute the fluxes
        for i in inputs.range_arr:

            # Compute the emissivity
            if inputs.single_arr[i]:
                tem = model[inputs.temp_id_arr[i]] if inputs.temp_eq_check[i] else tem_EQDB[inputs.eq_tem_arr[i]](model[inputs.temp_id_arr[i]])
                den = model[inputs.den_id_arr[i]] if inputs.den_eq_check[i] else den_EQDB[inputs.eq_den_arr[i]](model[inputs.den_id_arr[i]])
                emis = emis_interp[inputs.labels[i]](tem, den)

                # Compute the flux
                flux = FLUX_EQUATION_DICT[inputs.eq_flux_arr[i]](abund=model[inputs.ion_arr[i]], emis=emis,
                                                                 flambda=inputs.flambda_arr[i], cHbeta=cHbeta)


            else:
                flux_terms = []
                merge_in = merge_d[inputs.labels[i]]
                for j in merge_in.range_arr:
                    tem = model[merge_in.temp_id_arr[j]] if merge_in.temp_eq_check[j] else tem_EQDB[merge_in.eq_tem_arr[j]](model[merge_in.temp_id_arr[j]])
                    den = model[merge_in.den_id_arr[j]] if merge_in.den_eq_check[j] else den_EQDB[merge_in.eq_den_arr[j]](model[merge_in.den_id_arr[j]])
                    emis = emis_interp[merge_in.labels[j]](tem, den)
                    flux_terms.append(FLUX_EQUATION_DICT[merge_in.eq_flux_arr[j]](abund=model[merge_in.ion_arr[j]], emis=emis,
                                                                                  flambda=inputs.flambda_arr[i], cHbeta=cHbeta))
                # convert from log to linear, sum, convert back to log
                flux = tt.log10(tt.sum(tt.pow(10, tt.stack(flux_terms))))

            theo_flux = tt.inc_subtensor(theo_flux[i], flux)

        # Stored the fluxes and input fluxes, uncertainty
        pm.Deterministic('theo_flux', theo_flux, dims='lines')

        # Likelihood
        pm.Normal("likelihood", mu=theo_flux, sigma=inputs.log_err_arr, observed=inputs.log_flux_arr, dims='lines')

    return model




def run_model(model, draws=1000, tune=2000, target_accept=0.8, chains=8, cores=8,
              nuts_sampler='numpyro', callback=None):


    '''
        nuts_sampler : str, optional
        Which NUTS implementation to run. One of ["pymc", "nutpie", "blackjax", "numpyro"].
        This requires the chosen sampler to be installed.
        All samplers, except "pymc", require the full model to be continuous.
        
        If ``None`` (default), "nutpie" is used if installed and can be compiled to the desired backend.
            backend: str, optional.
        Which computational backend to use. Recommended to be one of "numba", "c", and "jax".
        May require installing extra dependencies.

    '''

    with model:
        time_it = time()
        trace = pm.sample(draws=draws, tune=tune, target_accept=target_accept, chains=chains, cores=cores,
                          nuts_sampler=nuts_sampler, progressbar='combined')
        time_it = time() - time_it
        print(f'-- Complete: {time_it:.1f} seconds' if time_it < 60 else f'Sampling time: {time_it / 60:.1f} minutes')

        print(f'- Drawing prior samples')
        prior = pm.sample_prior_predictive(draws=1000)

    # Save the prior
    trace['prior'] = prior['prior']
    if 'prior_predictive' in prior.children:
        trace['prior_predictive'] = prior['prior_predictive']

    return trace