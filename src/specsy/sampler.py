import pymc as pm
import numpy as np
from pytensor import tensor as tt
from specsy.models.fluxes_line import FLUX_EQUATION_DICT


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


    # Convenience arrays
    # # range_arr = np.arange(inputs.labels.size)
    # single_arr = np.ones(inputs.labels.size).astype(bool) if inputs.merged_dict is None else
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
                # flux = 0
                # merge_in = merge_d[inputs.labels[i]]
                # for j in merge_in.range_arr:
                #     tem = model[merge_in.temp_id_arr[j]] if merge_in.temp_eq_check[j] else tem_EQDB[merge_in.eq_tem_arr[j]](model[merge_in.temp_id_arr[j]])
                #     den = model[merge_in.den_id_arr[j]] if merge_in.den_eq_check[j] else den_EQDB[merge_in.eq_den_arr[j]](model[merge_in.den_id_arr[j]])
                #     emis = emis_interp[merge_in.labels[j]](tem, den)
                #     flux = flux + FLUX_EQUATION_DICT[merge_in.eq_flux_arr[j]](abund=model[merge_in.ion_arr[j]], emis=emis,
                #                                                               flambda=inputs.flambda_arr[i], cHbeta=cHbeta)
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


def direct_method_multi_region_orig(inputs, emis_interp, prior_dict, tem_EQDB, den_EQDB):

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
    eq_flux_arr = inputs.eq_flux_arr

    # Convenience arrays
    range_arr = np.arange(labels_arr.size)
    temp_eq_check = eq_tem_arr == '-'
    den_eq_check = eq_den_arr == '-'
    unique_species = np.unique(ion_arr)
    unique_params = np.unique(np.concatenate([temp_id_arr, den_id_arr]))

    # PyMC model
    with (pm.Model(coords={"lines": labels_arr}) as model):

        # Save input fluxes
        pm.Data('input_flux', inputs.flux_arr, dims='lines')
        pm.Data('input_err', inputs.err_arr, dims='lines')

        # Container to store the models
        theo_flux = tt.zeros(labels_arr.size)

        # Compile the abundances
        for ion in unique_species:
            if ion != 'H1':
                set_prior(ion, prior_dict, abund_type=True, name_param=ion)
        pm.Data('H1', 1.0)

        # Extinction
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

            # Compute the flux
            flux = FLUX_EQUATION_DICT[eq_flux_arr[i]](abund=model[ion_arr[i]], emis=emis,
                                                      flambda=flambda_arr[i], cHbeta=cHbeta)

            theo_flux = tt.inc_subtensor(theo_flux[i], flux)

        # Stored the fluxes and input fluxes, uncertainty
        pm.Deterministic('theo_flux', theo_flux, dims='lines')

        # Likelihood
        pm.Normal("likelihood", mu=theo_flux, sigma=input_err, observed=input_obs, dims='lines')

    return model


def run_model(model, draws=1000, tune=2000, target_accept=0.8, chains=8, cores=8, nuts_sampler='numpyro',
              callback=None):

    nuts_sampler_kwargs = None if nuts_sampler != 'nutpie' else {"backend": "jax", 'gradient_backend': "jax"}

    with model:
        trace = pm.sample(draws=draws, tune=tune, target_accept=target_accept, chains=chains, cores=cores,
                          nuts_sampler=nuts_sampler, callback=callback, progressbar='combined',
                          nuts_sampler_kwargs=nuts_sampler_kwargs)

    return trace