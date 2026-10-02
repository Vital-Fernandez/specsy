from logging import getLogger
_logger = getLogger('SpecSy')

import astropy.units as u
import numpy as np
import pymc as pm
import pytensor.tensor as pt
import extinction
import emcee
import inspect
import os
from dust_extinction.averages import G03_LMCAvg, G03_SMCBar
from dust_extinction.parameter_averages import G23
from scipy import interpolate
from scipy.optimize import minimize
from scipy.special import logsumexp
from specsy.io import specsy_cfg, SpecSyError


# SSPs variables
BPASS_Z = np.array(specsy_cfg['stellar']['ssp']['BPASS'].values()) #np.array([1e-5, 1e-4, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.014, 0.02, 0.03, 0.04])


# Nebular continuum constants, as in models.nebular_continuum().
# gamma holds emission coefficients for HI (free-free, bound-free, 2-photon)
# and HeI assuming He/H = 0.1, from Aller (1984) and Ferland (1980), for Case B
# conditions at T = 1e4 K with f_esc = 0.
gamma = np.array([0., 2.11e-4, 5.647, 9.35, 9.847, 10.582, 16.101, 24.681, 26.736,
                  24.883, 29.979, 6.519, 8.773, 11.545, 13.585, 6.333, 10.444, 7.023,
                  9.361, 7.59, 9.35, 8.32, 9.53, 8.87]) * 1e-40   # erg cm3 s-1 Hz-1
nebx = np.array([912., 913., 1300., 1500., 1800., 2200., 2855., 3331., 3421., 3422.,
                 3642., 3648., 5700., 7000., 8207., 8209., 14583., 14585., 22787., 22789.,
                 32813., 32815., 44680., 44682.])
alpha_B = 2.6e-13
Qbase = 52
L_SUN = 3.83e33   # erg s-1

sparse_nebcont = (2.998e18 * gamma * 10 ** Qbase) / (alpha_B * nebx ** 2)

# Extinction laws, each evaluated at E(B-V) = 1 so the result is A_lambda per
# unit reddening.  R_V is fixed per law, as in SESAMME v1.0.  The `extinction`
# package returns magnitudes directly; the dust_extinction models return
# A_lambda / A_V, hence the factor of R_V.
EXT_LAWS = {'CCM':           lambda w: extinction.ccm89(w, 3.1, 3.1),
            'ODonnell':      lambda w: extinction.odonnell94(w, 3.1, 3.1),
            'Fitzpatrick99': lambda w: extinction.fitzpatrick99(w, 3.1, 3.1),
            'FitzMassa07':   lambda w: extinction.fm07(w, 3.1),
            'Calzetti':      lambda w: extinction.calzetti00(w, 4.05, 4.05),
            'Gordon23':      lambda w: G23(Rv=3.1)(w * u.AA) * 3.1,
            'LMC':           lambda w: G03_LMCAvg()(w * u.AA) * G03_LMCAvg().Rv,
            'SMC':           lambda w: G03_SMCBar()(w * u.AA) * G03_SMCBar().Rv}


def extinction_curve(x, ext_law='CCM'):

    """A_lambda / E(B-V) in magnitudes, for one of SESAMME's extinction laws.

    Every law in EXT_LAWS is linear in A_V at the fixed R_V that SESAMME
    assumes, so the whole curve reduces to this one array, evaluated once here
    and reused at every sampling step as 10**(-0.4 * ebv * klam).

    Parameters
    ----------
    ext_law : str, optional
        Defaults to models.use_ext_law, so set_ext_law() still controls it.

    Raises
    ------
    ValueError
        If ext_law is not an implemented option.
    """

    if ext_law not in EXT_LAWS:
        raise ValueError("'" + str(ext_law) + "' is not a valid choice of extinction law.\n"
                         "Accepted values are " + ", ".join(EXT_LAWS) + ".")

    wave = np.asarray(x, dtype=np.float64)
    klam = np.asarray(EXT_LAWS[ext_law](wave), dtype=float)

    if np.any(klam < 0) or not np.all(np.isfinite(klam)):
        _logger.warning(f"Extinction curve for {ext_law} has negative or non-finite values; check the wavelength coverage"
                        " against the law's validity range.")

    return klam



def set_prior_bounds(prior_dict=None, prior_lowbounds=None, prior_highbounds=None):

    """
    Set boundaries on the priors, in the order log(age), log(Z), E(B-V), log(A).

    The age and metallicity bounds select which grid points are retained;
    the E(B-V) and log(A) bounds become uniform priors, as in SESAMME v1.0.

    Raises
    ------
    ValueError
        If any boundary pair is not increasing.
    """

    if prior_dict is None:
        prior_dict = dict.fromkeys(specsy_cfg['stellar']['ssp']['synthesis']['params'])
        prior_lowbounds = specsy_cfg['stellar']['ssp']['synthesis']['prior_lowbounds']
        prior_highbounds = specsy_cfg['stellar']['ssp']['synthesis']['prior_highbounds']

    prior_dict['age'] = [prior_lowbounds[0], prior_highbounds[0]]
    prior_dict['met'] = [prior_lowbounds[1], prior_highbounds[1]]
    prior_dict['ebv'] = [prior_lowbounds[2], prior_highbounds[2]]
    prior_dict['amp'] = [prior_lowbounds[3], prior_highbounds[3]]

    for key, value in prior_dict.items():
        if value[0] > value[1]:
            raise ValueError("Prior boundaries are out of order for variable '" + key + "'")

    return prior_dict


def grid_axes(binaries):
    """Return the model grid axes from a StellarBinaries object.

    `ages` are already log(age/yr); `metallicities` are linear mass fractions,
    so `mets` is their log10 for the priors and the posterior, while
    `met_vals` is kept for the get_spectrum() lookup.
    """
    ages = np.asarray(binaries.age_log, dtype=float)
    met_vals = np.asarray(binaries.metallicities, dtype=float)

    if np.any(met_vals <= 0):
        raise SpecSyError(f'Metallicities {met_vals[met_vals <= 0]} have no logarithm; '
                          f'drop them before building the grid.')

    return ages, np.log10(met_vals), met_vals


def nebular_shape(x):
    """Wavelength dependence of the nebular continuum, at unit amplitude and Q = Qbase.

    Identical for every grid point, since nebx and sparse_nebcont are fixed, so
    the interpolator is built once here rather than once per sampling step.
    Units are L_sun / A.
    """
    interp_function = interpolate.interp1d(nebx, sparse_nebcont, fill_value='extrapolate')

    return interp_function(np.asarray(x, dtype=float)) / L_SUN


def nebular_scale(ion_table, age_keys, met_keys):

    """
    Per-grid-point rescaling 10**(Q - Qbase), matching stellar_grid() ordering.

    Q is the log ionizing photon rate of the SSP and is the only part of the
    nebular continuum that varies across the grid.
    """

    n_met = len(met_keys)
    q = np.empty(len(age_keys) * n_met, dtype=float)

    for i, ak in enumerate(age_keys):
        for j, mk in enumerate(met_keys):
            Q_new = float(np.asarray(ion_table[ion_table['Z'] == mk][ak]).ravel()[0])
            q[i * n_met + j] = 10 ** (Q_new - Qbase)

    return q


def extinction_curve(x, ext_law='CCM'):

    """A_lambda / E(B-V) in magnitudes, for one of SESAMME's extinction laws.

    Every law in EXT_LAWS is linear in A_V at the fixed R_V that SESAMME
    assumes, so the whole curve reduces to this one array, evaluated once here
    and reused at every sampling step as 10**(-0.4 * ebv * klam).

    Parameters
    ----------
    ext_law : str, optional
        Defaults to models.use_ext_law, so set_ext_law() still controls it.

    Raises
    ------
    ValueError
        If ext_law is not an implemented option.
    """

    if ext_law not in EXT_LAWS:
        raise ValueError("'" + str(ext_law) + "' is not a valid choice of extinction law.\n"
                         "Accepted values are " + ", ".join(EXT_LAWS) + ".")

    wave = np.asarray(x, dtype=np.float64)
    klam = np.asarray(EXT_LAWS[ext_law](wave), dtype=float)

    if np.any(klam < 0) or not np.all(np.isfinite(klam)):
        _logger.warning(f"Extinction curve for {ext_law} has negative or non-finite values; check the wavelength coverage"
                        " against the law's validity range.")

    return klam


def normalize_grid(flat_m, mask=None):
    """Scale each grid row to unit median over the fitted pixels.

    The SSP rows span several dex in absolute luminosity, so a single shared
    log(A) cannot fit them all at once -- the marginalisation then collapses
    onto whichever cell happens to match that one amplitude. Dividing the
    spread out here keeps log(A) a genuine sampled parameter: it sets the
    overall flux scale, while each cell's intrinsic luminosity is carried in
    row_norm and added back when converting to mass.
    """
    sub = flat_m if mask is None else flat_m[:, mask]
    row_norm = np.median(sub, axis=1)

    if np.any(row_norm <= 0):
        raise SpecSyError(f'{np.count_nonzero(row_norm <= 0)} grid rows have a non-positive '
                          f'median over the fitted pixels; check the mask and the flux units.')

    return flat_m / row_norm[:, None], row_norm


def stellar_grid(binaries):

    """
    Stack every spectrum in `binaries` into a flat (n_age * n_met, n_wave) array.

    Row k corresponds to age index i and metallicity index j via
    k = i * met_vals.size + j, matching nebular_scale() and posterior_grid(). `binaries` is used as
    given — no library/alpha/imf selection is applied, only a warning if it turns out to hold more than
    one value for any of them.
    """

    if not binaries.uniform_dispersion:
        raise SpecSyError('The binaries do not share a common wavelength grid. '
                          'Call to_stellar_binaries(..., pixel_width=...) (or disp_intvl/pixel_number) '
                          'before building the model grid.')

    frame = binaries.frame

    # Flag (but don't act on) a mixed library, alpha or IMF: the grid stacks everything regardless
    for col in ('library', 'alpha', 'imf'):
        n_unique = frame[col].nunique(dropna=False)
        if n_unique > 1:
            _logger.warning(f'The binaries hold {n_unique} different "{col}" values: '
                            f'{sorted(frame[col].dropna().unique().tolist())}. All of them are stacked '
                            f'into the grid.')

    # Entries that share an (age, metallicity) pair cannot both occupy the same grid cell
    dup_mask = frame.duplicated(subset=['age', 'metallicity'], keep='last')
    if dup_mask.any():
        dups = frame.loc[dup_mask, ['age', 'metallicity']].drop_duplicates()
        _logger.warning(f'{int(dup_mask.sum())} of the {len(frame)} spectra share their (age, metallicity) '
                        f'with another entry: {list(dups.itertuples(index=False, name=None))}. Only the last '
                        f'one of each pair is kept.')
        frame = frame.loc[~dup_mask]

    ages = np.unique(frame['age'])
    met_vals = np.unique(frame['metallicity'])

    if len(frame) != ages.size * met_vals.size:
        raise SpecSyError(f'The binaries hold {len(frame)} spectra but a full grid needs '
                          f'{ages.size} x {met_vals.size} = {ages.size * met_vals.size}: some (age, '
                          f'metallicity) combinations are missing.')

    # Row order k = i * met_vals.size + j: sort by (age, metallicity) to match that layout
    ordered_idcs = frame.sort_values(['age', 'metallicity']).index
    rows = np.stack([binaries.flux_series.at[idx] for idx in ordered_idcs])

    return rows, ages, met_vals


def build_grid(x, binaries, ion_table=None, ext_law='CCM', add_nebular=False, age_keys=None, met_keys=None,
               prior_bounds=None):

    """Assemble everything the PyMC model needs, in one numpy pass.

    Parameters
    ----------
    x : array-like
        Wavelength array; normally binaries.wave_series.at[0].
    binaries : StellarBinaries
        Source of the SSP spectra. Every spectrum it holds is stacked into the grid, as is.
    age_keys, met_keys : list of str, optional
        Only needed when add_nebular is True: the ion_table column and row keys, which are a naming
        convention of that table rather than of the model grid.
    prior_bounds : dict, optional
        Prior (low, high) bounds per parameter, e.g. {'age': [6.0, 7.5], 'amp': [...], 'ebv': [...],
        'met': [-5.0, -1.4]}. Only 'age' (log(age/yr)) and 'met' (log10(Z)) are checked here: a warning
        is logged if either prior range extends beyond what `binaries` actually covers, since the model
        grid cannot constrain the prior past its own coverage. Other keys ('amp', 'ebv', ...) are ignored.

    Returns
    -------
    flat : ndarray, shape (n_age * n_met, n_wave)
        Stellar (+ nebular) model at unit amplitude, unreddened.
    klam : ndarray, shape (n_wave,)
    ages, mets : ndarray
        The grid axes, as log(age/yr) and log10(Z).
    """

    flat, ages, mets = stellar_grid(binaries)

    if prior_bounds is not None:

        # met_vals from stellar_grid are the grid's raw metallicity values (Z); prior bounds are log10(Z)
        grid_bounds = {'age': (float(ages.min()), float(ages.max())),
                       'met': (float(np.log10(mets.min())), float(np.log10(mets.max())))}

        for key, (grid_low, grid_high) in grid_bounds.items():
            if key not in prior_bounds:
                continue
            prior_low, prior_high = prior_bounds[key]
            if (prior_low < grid_low) or (prior_high > grid_high):
                _logger.warning(f'The "{key}" prior bounds ({prior_low:.4g}, {prior_high:.4g}) extend beyond the '
                                f'StellarBinaries coverage ({grid_low:.4g}, {grid_high:.4g}). The prior is not '
                                f'fully constrained by the model grid over that range.')

    if add_nebular:
        if (ion_table is None) or (age_keys is None) or (met_keys is None):
            raise SpecSyError('add_nebular requires an ion_table plus its age_keys and '
                              'met_keys; StellarBinaries carries no ionizing photon rate.')
        q = nebular_scale(ion_table, age_keys, met_keys)
        flat = flat + q[:, None] * nebular_shape(x)[None, :]

    if not np.all(np.isfinite(flat)):
        n_bad = np.count_nonzero(~np.isfinite(flat))
        _logger.warning(f'Model grid contains {n_bad} non-finite values; check the model suite coverage.')

    return flat, extinction_curve(x, ext_law), ages, mets



def build_model(x, y, yerr, flat, klam, mask, prior_dict, p_grid=None):

    """Construct the PyMC model with the grid index marginalised out.

    Parameters
    ----------
    flat, klam : ndarray
        Output of build_grid().
    mask : array-like
        Boolean array marking wavelength bins to keep.
    p_grid : ndarray, optional
        Prior weights over the flattened grid.  Defaults to uniform, matching
        SESAMME's flat priors on age and metallicity.

    Returns
    -------
    pymc.Model
    """

    y_m = np.asarray(y)[mask]
    err_m = np.asarray(yerr)[mask]
    klam_m = np.asarray(klam)[mask]

    # Rows scaled to unit median so a single ampl is meaningful across all cells;
    # the intrinsic luminosity spread is carried in row_norm instead
    flat_m, row_norm = normalize_grid(np.asarray(flat)[:, mask])
    _logger.debug(f'Grid rows normalised to unit median; row_norm spans '
                  f'{np.log10(row_norm.max() / row_norm.min()):.2f} dex')

    n_models = flat_m.shape[0]
    if p_grid is None:
        p_grid = np.full(n_models, 1.0 / n_models)
    logp_grid = np.log(p_grid)

    with pm.Model() as model:

        pm.Data("row_norm", row_norm)              # lands in idata.constant_data
        pm.Data("klam_mean", klam_m.mean())        # <-- NEW, next to row_norm

        ebv = pm.Uniform("ebv", *prior_dict['ebv'])
        ampl = pm.Uniform("ampl", *prior_dict['amp'])

        # Error rescaling: the quoted uncertainties understate the model-data
        # mismatch, which otherwise makes the grid posterior far too confident
        lnf = pm.Normal("lnf", 0., 1.)
        err_eff = err_m * pt.exp(lnf)

        # Centring the curve leaves ebv controlling only the slope of the
        # attenuation across the band; the overall level is ampl's job, and
        # sharing it between the two is what correlates them
        klam_c = klam_m - klam_m.mean()            # <-- NEW
        atten = pt.pow(10.0, -0.4 * ebv * klam_c)  # <-- was klam_m

        mu = pt.pow(10.0, ampl) * flat_m * atten         # (n_models, n_wave)
        r = (y_m[None, :] - mu) / err_eff[None, :]
        ll = -0.5 * pt.sum(r ** 2 + 2 * pt.log(err_eff), axis=1)

        logw = ll + logp_grid
        pm.Potential("marginal_likelihood", pt.logsumexp(logw))
        pm.Deterministic("p_k", pt.special.softmax(logw[None, :])[0])

    return model



def run_sesamme_pymc(x, y, yerr, binaries, ion_table, mask, ext_law='CCM', add_nebular=False,
                     draws=2000, tune=2000, chains=4, target_accept=0.9, **sample_kwargs):

    prior_dict = set_prior_bounds()

    # ampl is a flux scale against unit-median rows, not SESAMME's mass-scale
    # log(A), so its bounds come from the data rather than from the config
    centre = np.log10(np.median(np.asarray(y)[mask]))
    prior_dict['amp'] = [centre - 3, centre + 3]

    flat, klam, ages, mets = build_grid(x, binaries, ion_table=ion_table, ext_law=ext_law,
                                        add_nebular=add_nebular,
                                        age_bounds=prior_dict['age'],
                                        met_bounds=prior_dict['met'])
    _logger.info(f'Grid built with {flat.shape} spectra and {np.sum(mask)} pixels')

    model = build_model(x, y, yerr, flat, klam, mask, prior_dict)

    with model:
        idata = pm.sample(draws=draws, tune=tune, chains=chains,
                          target_accept=target_accept, **sample_kwargs)

    return idata, {"ages": ages, "mets": mets, "flat": flat, "klam": klam}


#########################
# Initial values
#########################

def marginal_log_likelihood_profiled(ebv, y, yerr, mask, flat, klam, logp_grid):
    """Marginal log-likelihood with the amplitude profiled out per grid cell.

    For a Gaussian likelihood with a purely multiplicative scale, the optimal
    A_k has a closed form, so the amplitude need not be sampled at all:
    cells are then compared on spectral shape, not absolute normalisation.
    """
    y_m, err_m = y[mask], yerr[mask]
    mod_m = flat[:, mask] * 10 ** (-0.4 * ebv * klam[mask])

    w = 1.0 / err_m ** 2
    A_k = (mod_m @ (y_m * w)) / (mod_m ** 2 @ w)          # (n_models,)

    r = (y_m[None, :] - A_k[:, None] * mod_m) / err_m[None, :]
    ll = -0.5 * np.sum(r ** 2, axis=1)

    return logsumexp(ll + logp_grid), ll, A_k


def negative_marginal_log_likelihood(params, *args):
    """Named negation of the profiled marginal likelihood, for minimize().

    params is length-1: only E(B-V) is searched, since the amplitude is
    profiled out analytically per grid cell.
    """
    logp, _, _ = marginal_log_likelihood_profiled(params[0], *args)

    return -logp


def get_initial_values(y, yerr, mask, flat, klam, ages, mets, p_grid=None,
                       initial=(0.1, 0.0), bounds=((0.01, 1.0), (-20., 1.)),
                       return_grid=True, show_solution=True):
    """Find a starting point for the sampler, in SESAMME's [age, met, ebv, ampl] order.

    Only E(B-V) is searched numerically: for a Gaussian likelihood with a
    multiplicative scale, the optimal amplitude per grid cell has a closed
    form, so it is profiled out rather than optimized. The returned ampl is
    converted into the normalised convention build_model() samples in, and
    the age and metallicity are those of the most probable cell at the
    optimized E(B-V) -- a report rather than a starting position, since NUTS
    never samples a grid index.

    Parameters
    ----------
    flat, klam, ages, mets : ndarray
        Output of build_grid().
    p_grid : ndarray, optional
        Prior weights over the flattened grid; defaults to uniform, matching
        build_model()'s default.
    initial : tuple
        Starting guess (ebv, ampl); only the first entry is used.
    bounds : tuple of tuples
        (ebv_bounds, ampl_bounds); only the first is used by the search.
    return_grid : bool
        If True (default), return [age, met, ebv, ampl] -- matching the order
        and shape of the original emcee-based function. If False, return just
        [ebv, ampl].

    Returns
    -------
    initial_optimized : list
        [age, met, ebv, ampl] if return_grid, else [ebv, ampl].
    """
    y, yerr = np.asarray(y), np.asarray(yerr)

    n_models = flat.shape[0]
    if p_grid is None:
        p_grid = np.full(n_models, 1.0 / n_models)
    logp_grid = np.log(p_grid)

    solution = minimize(negative_marginal_log_likelihood, [initial[0]],
                        method='Nelder-Mead', bounds=[bounds[0]],
                        args=(y, yerr, mask, flat, klam, logp_grid),
                        options={'maxiter': 8000})

    if show_solution:
        print(solution)

    ebv_opt = float(solution.x[0])

    _, ll, A_k = marginal_log_likelihood_profiled(ebv_opt, y, yerr, mask, flat, klam, logp_grid)
    k_map = int(np.argmax(ll + logp_grid))

    if A_k[k_map] <= 0:
        raise SpecSyError(f'Optimal amplitude for the best cell is non-positive '
                          f'({A_k[k_map]:.3e}); check the flux units and the mask.')

    # A_k multiplies the raw row; build_model samples against unit-median rows
    _, row_norm = normalize_grid(flat[:, mask])
    ampl_opt = float(np.log10(A_k[k_map] * row_norm[k_map]))

    if not return_grid:
        return [ebv_opt, ampl_opt]

    i_map, j_map = np.unravel_index(k_map, (ages.size, mets.size))

    return [float(ages[i_map]), float(mets[j_map]), ebv_opt, ampl_opt]


#########################
# Posterior interpretation
#########################

def posterior_grid(idata, grid_info):
    """Collapse the marginalised index posterior into age and metallicity PDFs.

    Returns
    -------
    dict with keys:
        'joint'               : (n_age, n_met) posterior over the grid
        'age_pdf', 'met_pdf'  : marginals over log(age) and log10(Z)
        'age_map', 'met_map'  : most probable grid point
        'age_mean', 'met_mean': probability-weighted means
        'logA'                : cluster amplitude in SESAMME's convention
    """
    ages, mets = grid_info["ages"], grid_info["mets"]

    p_k = idata.posterior["p_k"].mean(dim=("chain", "draw")).values
    joint = p_k.reshape(ages.size, mets.size)

    age_pdf = joint.sum(axis=1)
    met_pdf = joint.sum(axis=0)
    i_map, j_map = np.unravel_index(np.argmax(joint), joint.shape)

    row_norm = idata.constant_data["row_norm"].values
    klam_mean = float(idata.constant_data["klam_mean"].values)

    ampl = idata.posterior["ampl"].values.ravel()
    ebv = idata.posterior["ebv"].values.ravel()

    # SESAMME's log(A). Three corrections to the sampled ampl:
    #   - row_norm      : rows were scaled to unit median before fitting
    #   - klam_mean term: the extinction curve is centred, so ampl is the level
    #                     at the band-mean attenuation rather than unreddened
    logA_k = (ampl[:, None]
              - np.log10(row_norm)[None, :]
              + 0.4 * ebv[:, None] * klam_mean)
    logA = float(np.sum(p_k * logA_k.mean(axis=0)))

    return {"joint": joint,
            "age_pdf": age_pdf,
            "met_pdf": met_pdf,
            "age_map": ages[i_map],
            "met_map": mets[j_map],
            "age_mean": float(np.sum(age_pdf * ages)),
            "met_mean": float(np.sum(met_pdf * mets)),
            "logA": logA}


def credible_cells(joint, level=0.68):
    """Smallest set of grid cells whose summed probability exceeds `level`.

    The discrete analogue of a credible interval.  Returns (i, j) index pairs;
    the set may be disjoint if the posterior is multimodal, which is expected
    behaviour for a genuinely discrete PDF.
    """
    order = np.argsort(joint.ravel())[::-1]
    cum = np.cumsum(joint.ravel()[order])
    keep = order[:np.searchsorted(cum, level) + 1]

    return np.array(np.unravel_index(keep, joint.shape)).T





def _snap_to_grid(values, grid):
    """
    Snap each value to its nearest entry in a 1D grid.

    Uses np.searchsorted to bracket each value between its two neighbouring
    grid points, then picks whichever of the two is closer.
    """

    grid = np.sort(np.asarray(grid, dtype=float))
    idx = np.searchsorted(grid, values, side='left')
    idx = np.clip(idx, 1, grid.size - 1)

    left, right = grid[idx - 1], grid[idx]
    idx = np.where(np.abs(values - left) <= np.abs(right - values), idx - 1, idx)

    return grid[idx]


def _apply_snap_grid(samples, snap_grid):

    if snap_grid:
        for p, grid in snap_grid.items():
            if p in samples:
                samples[p] = _snap_to_grid(samples[p], grid)

    return samples


def load_emcee_run(filename, runname, param_names=('age', 'Z', 'ebv', 'logA'), discard=40, thin=None,
                    autocorr_tol=0, snap_grid=None):
    """
    Read one emcee HDF5 backend and return a samples dict ready for plot_corner/plot_ssp_diag.

    filename:     path to the .h5 file
    runname:      the run's name inside the HDF5 file (the `name` argument to HDFBackend)
    param_names:  names for each chain column, in the order the sampler used them --
                  defaults to your (age, Z, ebv, logA) ordering
    discard:      burn-in steps to discard
    thin:         thinning step; if None, computed as int(tau.max() / 2) from the
                  autocorrelation time, matching your existing script
    autocorr_tol: tol passed to get_autocorr_time (0 disables the convergence check,
                  as in your current usage)
    snap_grid:    optional dict {param: grid_values}, for parameters drawn from a
                  discrete grid (age, Z). Snaps every sample to its nearest grid point
                  so the histograms/contours in plot_corner land exactly on grid bins,
                  matching how the fit itself looks up nearest_t/nearest_z.

    Returns a dict {param_name: 1D array}.
    """

    reader = emcee.backends.HDFBackend(filename, name=runname)

    if thin is None:
        tau = reader.get_autocorr_time(tol=autocorr_tol)
        thin = int(tau.max() / 2.)

    flat_samples = reader.get_chain(discard=discard, thin=thin, flat=True)

    if flat_samples.shape[1] != len(param_names):
        raise ValueError(f'Chain from "{runname}" has {flat_samples.shape[1]} parameters but '
                         f'{len(param_names)} param_names were given: {param_names}')

    samples = {name: flat_samples[:, i] for i, name in enumerate(param_names)}

    return _apply_snap_grid(samples, snap_grid)


def load_pymc_run(filename, param_names=('age', 'Z', 'ebv', 'logA'), var_names=None, group='posterior',
                  snap_grid=None):
    """
    Read a PyMC trace saved with `trace.to_netcdf(filename)` (an arviz InferenceData)
    and return a samples dict ready for plot_corner/plot_ssp_diag.

    filename:    path to the .nc file
    param_names: names to use as the output dict's keys, in the same order as var_names
    var_names:   the variable names as they appear in the InferenceData's `group`;
                 defaults to param_names, i.e. assumes your PyMC model's variable names
                 already match ('age', 'Z', 'ebv', 'logA') -- pass this explicitly if
                 your model used different names (e.g. var_names=['t', 'z', 'ebv', 'amp'])
    group:       which InferenceData group to read ('posterior' by default; use
                 'posterior_predictive' etc. if that's what you need)
    snap_grid:   as in load_emcee_run -- snaps discrete-grid parameters (age, Z) to
                 their nearest grid point after flattening chain/draw

    Returns a dict {param_name: 1D array}, with the chain and draw dimensions flattened
    together. If a variable has extra dimensions (e.g. per-population values), those are
    flattened into the same 1D array too -- pass one param_name per scalar variable if
    you need them kept separate.
    """

    import arviz as az

    idata = az.from_netcdf(filename)

    if not hasattr(idata, group):
        raise ValueError(f'Group "{group}" not found in {filename}; available groups: '
                         f'{list(idata.groups())}')

    data = getattr(idata, group)
    var_names = var_names or param_names

    if len(var_names) != len(param_names):
        raise ValueError('param_names and var_names must have the same length: got '
                         f'{len(param_names)} and {len(var_names)}')

    samples = {}
    for pname, vname in zip(param_names, var_names):
        if vname not in data:
            raise KeyError(f'Variable "{vname}" not found in the "{group}" group of {filename}; '
                           f'available variables: {list(data.data_vars)}')
        samples[pname] = data[vname].values.reshape(-1)

    return _apply_snap_grid(samples, snap_grid)


_EXTENSION_BACKENDS = {'.h5': 'emcee', '.hdf5': 'emcee', '.nc': 'pymc', '.netcdf': 'pymc'}

_LOADERS = {'emcee': load_emcee_run, 'pymc': load_pymc_run}


def _infer_backend(filename):

    ext = os.path.splitext(filename)[1].lower()
    backend = _EXTENSION_BACKENDS.get(ext)

    if backend is None:
        raise ValueError(f'Could not infer a backend from the extension "{ext}" of {filename}; '
                         f'pass backend="emcee" or backend="pymc" explicitly for this run')

    return backend


def load_run(spec, **common_kwargs):
    """
    Read a single run, auto-dispatching to load_emcee_run or load_pymc_run.

    spec: a dict of keyword arguments for the run (may include "backend": "emcee"/"pymc"
          to force the loader; otherwise it is inferred from the filename's extension --
          .h5/.hdf5 -> emcee, .nc/.netcdf -> pymc). Or, for backward compatibility, an
          (filename, runname) tuple, treated as an emcee run.
    common_kwargs: defaults shared across runs (filename excepted), e.g. param_names,
                   snap_grid, discard.

    Returns a samples dict, as from load_emcee_run / load_pymc_run.
    """

    if isinstance(spec, dict):
        kwargs = {**common_kwargs, **spec}
    else:
        filename, runname = spec
        kwargs = {**common_kwargs, 'filename': filename, 'runname': runname}

    backend = kwargs.pop('backend', None) or _infer_backend(kwargs['filename'])
    loader = _LOADERS[backend]

    # Keep only the keywords this particular loader accepts (emcee- and pymc-specific
    # keywords, e.g. runname vs var_names, can otherwise collide in common_kwargs)
    accepted = set(inspect.signature(loader).parameters)
    kwargs = {k: v for k, v in kwargs.items() if k in accepted}

    return loader(**kwargs)


def load_runs(runs, **common_kwargs):
    """
    Read several runs at once, of any mix of backends (emcee HDF5, PyMC/arviz netCDF).

    runs: list of dict specs (or (filename, runname) tuples for emcee, kept for
          backward compatibility) -- see load_run.
    common_kwargs: defaults applied to every run, e.g. param_names, snap_grid.

    Returns a list of samples dicts, in the same order as `runs` -- pass it straight
    as `samples=` to plot_corner / plot_ssp_diag, with `set_labels` giving each run's name.
    """

    return [load_run(run, **common_kwargs) for run in runs]


# Backward-compatible alias for the previous emcee-only entry point
load_emcee_runs = load_runs