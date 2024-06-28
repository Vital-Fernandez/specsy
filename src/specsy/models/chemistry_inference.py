import numpy as np
import pymc as pm
from pytensor import tensor as tt


def storeValueInTensor(idx, value, tensor1D):
    return tt.inc_subtensor(tensor1D[idx], value)


def set_prior(param, prior_dict, abund_type=False, name_param=None, total_regions=1):

    # Read distribution configuration
    dist_name = prior_dict[param][0]
    dist_loc, dist_scale = prior_dict[param][1], prior_dict[param][2]
    dist_norm, dist_reLoc = prior_dict[param][3], prior_dict[param][4]

    # Load the corresponding probability distribution
    probDist = getattr(pm, dist_name)

    if abund_type:
        priorFunc = probDist(name_param, dist_loc, dist_scale) * dist_norm + dist_reLoc

    elif probDist.__name__ in ['HalfCauchy']:  # These distributions only have one parameter
        priorFunc = probDist(param, dist_loc, shape=total_regions) * dist_norm + dist_reLoc

    elif probDist.__name__ == 'Uniform':

        if param == 'logOH':
            priorFunc = pm.Bound(pm.Normal, lower=7.1, upper=9.1)('logOH', mu=8.0, sigma=1.0, testval=8.1)
        if param == 'logU':
            priorFunc = pm.Bound(pm.Normal, lower=-4.0, upper=-1.5)('logU', mu=-2.75, sigma=1.5, testval=-2.75)
        if param == 'logNO':
            priorFunc = pm.Bound(pm.Normal, lower=-2.0, upper=0.0)('logNO', mu=-1.0, sigma=0.5, testval=-1.0)

    else:
        priorFunc = probDist(param, dist_norm, dist_scale, shape=total_regions) * dist_norm + dist_reLoc

    return priorFunc


# def temperature_selection(self, fit_T_low=True, fit_T_high=True):
#
#     if self.lowTemp_check and fit_T_low:
#         self.set_prior('T_low')
#
#         if self.highTemp_check:
#             self.set_prior('T_high')
#         else:
#             self.prior_vars['T_high'] = TOIII_from_TSIII_relation(self.prior_vars['T_low'])
#
#     else:
#         if self.highTemp_check and fit_T_high:
#             self.set_prior('T_high')
#
#         self.prior_vars['T_low'] = TOII_from_TOIII_relation(self.prior_vars['T_high'], self.prior_vars['n_e'])
#
#     return

def temperature_selection(Tlow_diag, Thigh_diag, prior_dict):

    # Both diagnostics available
    if (not callable(Tlow_diag)) and (not callable(Thigh_diag)):
        Tlow_prior, Thigh_prior = set_prior('T_low', prior_dict), set_prior('T_high', prior_dict)

    # Only one:
    else:
        if callable(Tlow_diag):
            Thigh_prior = set_prior('T_high', prior_dict)
            Tlow_prior = Tlow_diag(Thigh_prior)

        else:
            Tlow_prior = set_prior('T_low', prior_dict)
            Thigh_prior = Thigh_diag(Tlow_prior)

    return Tlow_prior, Thigh_prior


def direct_method_inference(fname, inputs, prior_dict, idcs_highTemp_ions, emiss_interp, eq_tt,
                            Tlow_diag, Thigh_diag):

    # Container synthetic fluxes # FIXME do I need this one for loop inferences
    prior_vars = {}

    # Unpack the inputs
    line_arr, ion_arr = inputs.lines, inputs.particles
    flux_arr, err_arr = inputs.fluxes, inputs.errs
    lambda_arr = inputs.f_lambda

    # Define observable input
    fluxTensor = tt.zeros(line_arr.size)
    inputFlux = np.log10(flux_arr)
    inputFluxErr = np.log10(1 + err_arr / flux_arr)

    # Define the counters for loops
    linesRangeArray = np.arange(line_arr.size)

    # Unique ions for the analysis
    ions_unique = np.unique(ion_arr)
    ions_unique = ions_unique[ions_unique != 'H1']

    # Assign variable values
    prior_vars['H1'] = 0.0

    with pm.Model() as model:

        # Declare models parameters priors
        prior_vars['n_e'] = set_prior('n_e', prior_dict)
        prior_vars['cHbeta'] = set_prior('cHbeta', prior_dict)

        # Establish models temperature structure
        # prior_vars['T_low'], prior_vars['T_high'] = temperature_selection(fit_Tlow, fit_Thigh, prior_dict)
        # prior_vars['T_low'] = set_prior('T_low', prior_dict)
        # prior_vars['T_high'] = set_prior('T_high', prior_dict)
        # Both diagnostics available
        if (not callable(Tlow_diag)) and (not callable(Thigh_diag)):
            prior_vars['T_low'], prior_vars['T_high'] = set_prior('T_low', prior_dict), set_prior('T_high', prior_dict)

        # Only one:
        else:
            if callable(Tlow_diag):
                prior_vars['T_high'] = set_prior('T_high', prior_dict)
                prior_vars['T_low'] = Tlow_diag(prior_vars['T_high'])

            else:
                prior_vars['T_low'] = set_prior('T_low', prior_dict)
                prior_vars['T_high'] = Thigh_diag(prior_vars['T_low'])

        # Define grid interpolation variables
        emisCoord_low = tt.stack([[prior_vars['T_low'][0]], [prior_vars['n_e'][0]]], axis=-1)
        emisCoord_high = tt.stack([[prior_vars['T_high'][0]], [prior_vars['n_e'][0]]], axis=-1)

        # Establish models composition
        for ion in ions_unique:
            prior_vars[ion] = set_prior(ion, prior_dict, abund_type=True, name_param=ion)

        # Loop through the lines to compute the synthetic fluxes
        for i in linesRangeArray:

            # Declare line properties
            lineLabel = line_arr[i]
            lineIon = ion_arr[i]
            lineFlambda = lambda_arr[i]

            # Compute emisivity for the corresponding ion temperature
            T_calc = emisCoord_high if idcs_highTemp_ions[i] else emisCoord_low
            line_emis = emiss_interp[lineLabel](T_calc)

            # Declare fluorescence correction
            lineftau = 0.0

            # Compute line flux
            lineFlux_i = eq_tt.compute_flux(lineLabel,
                                            line_emis[0][0],
                                            prior_vars['cHbeta'],
                                            lineFlambda,
                                            prior_vars[lineIon],
                                            lineftau,
                                            O3=prior_vars['O3'],
                                            T_high=prior_vars['T_high'])

            # Assign the new value in the tensor
            fluxTensor = storeValueInTensor(i, lineFlux_i[0], fluxTensor)

        # Store computed fluxes
        pm.Deterministic('calcFluxes_Op', fluxTensor)

        # Likelihood gas components
        pm.Normal('Y_emision', mu=fluxTensor, sigma=inputFluxErr, observed=inputFlux)

        # Run the data
        inference_data = pm.sample(1000, tune=2000, chains=8, cores=8, init='auto', progressbar=True)

        #package_results(fname, inference_data, prior_dict, true_values)

    return inference_data


def chemical_model(transition_list, obs_mean_array, obs_std_array, element_array, temperature_array, function_dict):

    unique_elements = np.unique(element_array)
    unique_temperatures = np.unique(temperature_array)

    coords = {'element': unique_elements, 'region': unique_temperatures}
    with pm.Model(coords=coords) as model:

        # Model data
        flux_obs = pm.Data('observed_data', obs_mean_array)
        flux_err = pm.Data('observed_uncertainty', obs_std_array)

        # Priors definition
        abundance_priors = pm.Normal('abundance', 5, 5, dims='element')
        temperature_priors = pm.Normal('temperature', 15000, 5000,  dims='region')
        density = pm.HalfCauchy('density', 2.0, 0.0, shape=1) + 200
        cHbeta = pm.HalfCauchy('density', 2.0, 0.0, shape=1)

        # Loop through the input transitions and compute the likelihoods:
        for i, transition in enumerate(transition_list):

            # Get the transition specific parameters
            transition_temperature = temperature_priors[element_array[i]]
            transition_abundance = abundance_priors[element_array[i]]

            transition_emission = function_dict[transition](transition_temperature, density, cHbeta)
            flux_theo = pm.Deterministic(transition, transition_abundance * transition_emission)

            likelihood_i = pm.Normal(transition, mu=flux_theo, sigma=flux_obs[i], observed=flux_err[i])

    return model
