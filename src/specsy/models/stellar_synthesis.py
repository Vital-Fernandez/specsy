from logging import getLogger

_logger = getLogger('SpecSy')

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
from pytensor.compile.ops import as_op
import arviz as az
import xarray as xr
from scipy.stats import truncnorm
import matplotlib.pyplot as plt

from specsy.io import SpecSyError, specsy_cfg, load_ssp_trace_results
from specsy.models.nebular_continuum import nebular_continuum, ion_table_log_q, load_ionization_frame, NEB_GAMMA_WAVE, NEB_GAMMA_COEF
from specsy.models.ssp import SSPGrid
from specsy.models.extinction import extinction_curve, DEFAULT_RV

PARAM_NAMES = ('log_age', 'log_Z', 'ebv', 'log_A')
MCMC_SAMPLERS = ('DEMetropolis', 'DEMetropolisZ', 'Metropolis')


def _cell_widths(grid, lo, hi):
    """Width of each node's nearest-neighbour cell inside [lo, hi] = the implicit prior mass of that node."""
    if grid.size == 1:
        return np.array([hi - lo], dtype=float)
    edges = np.clip(np.concatenate([[lo], 0.5 * (grid[1:] + grid[:-1]), [hi]]), lo, hi)
    return np.diff(edges)


def sample_truncnorm(mu, lower, upper, n, stellar_binaries, spread=0.10, ebv_zero=0.5, seed=None, atol=1e-6):

    mu = np.array(mu, dtype=float)
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    mu_in = mu.copy()

    # log(age) and log10(Z): a value on the edge node (inside the bounds) moves one node inwards
    grids = (np.asarray(stellar_binaries.ages, float), np.log10(np.asarray(stellar_binaries.metallicities, float)))
    for i, grid in enumerate(grids):
        nodes = grid[(grid >= lower[i] - atol) & (grid <= upper[i] + atol)]
        if nodes.size > 1:
            idx = SSPGrid._nearest(nodes, mu[i])[0]
            mu[i] = nodes[1] if idx == 0 else nodes[-2] if idx == nodes.size - 1 else mu[i]

    # E(B-V) at zero
    if np.isclose(mu[2], 0.0, atol=atol):
        mu[2] = ebv_zero

    # log(A) on a prior limit: one unit inwards
    if np.isclose(mu[3], lower[3], atol=atol):
        mu[3] = lower[3] + 1
    elif np.isclose(mu[3], upper[3], atol=atol):
        mu[3] = upper[3] - 1

    moved = np.flatnonzero(mu != mu_in)
    if moved.size:
        _logger.info('Initial values moved from the grid/prior edges: '
                     + ', '.join(f'{PARAM_NAMES[i]} {mu_in[i]:.4f} -> {mu[i]:.4f}' for i in moved))

    # Relative spread, absolute fallback for parameters still at zero (e.g. log(A) = 0)
    sigma = spread * np.abs(mu)
    sigma[sigma == 0] = spread

    a, b = (lower - mu) / sigma, (upper - mu) / sigma
    samples = truncnorm.rvs(a, b, loc=mu, scale=sigma, size=(n, mu.size), random_state=seed)

    return [dict(zip(PARAM_NAMES, walker)) for walker in samples], mu


def ssp_inputs(idata):

    """
    SSPGrid rebuilt from the trace's constant_data group: observation, extinction and (visited) grid arrays, plus
    the node lookup built on the full age and Z axes, so nearest-node matching and cell widths match the original
    grid.
    """

    return SSPGrid.from_idata(idata)


def marginal_stats(idata, param_names=PARAM_NAMES):

    theta = flat_samples(idata)
    table = {name: np.percentile(theta[:, i], (16, 50, 84)) for i, name in enumerate(param_names)}

    post = idata.posterior
    n_samples = post.sizes['chain'] * post.sizes['draw']
    ess = az.ess(idata, var_names=list(param_names), method='bulk')
    for name in param_names:
        ess_i = float(ess[name])
        table[name] = np.append(table[name], [ess_i, n_samples / ess_i])

    cols = ['p16', 'p50', 'p84'] + (['ess', 'tau'] if len(next(iter(table.values()))) == 5 else [])

    return pd.DataFrame.from_dict(table, orient='index', columns=cols)


def node_rows(d, log_age, log_met):
    """Flux row of the grid node nearest to each (log_age, log_met), -1 if the combination is missing."""
    return d.rows(log_age, log_met)


def ssp_model(theta, d):
    return d.model(theta)


def ssp_components(theta, d):
    """Stellar and nebular parts of the rescaled, reddened model (nebular is None if it was not fitted)."""
    return d.components(theta)


def flat_samples(idata):
    post = idata.posterior
    return np.column_stack([post[name].values.ravel() for name in PARAM_NAMES])


def node_table(idata, uniform_node_prior=False, d=None):
    """
    Posterior mass of every visited (age, Z) grid node plus the conditional medians of E(B-V) and log(A)
    within it. The likelihood is constant inside each node's cell, so the chain is a discrete mixture.
    uniform_node_prior=True divides out the cell widths, i.e. every node gets the same prior weight instead of
    one proportional to its (age, log Z) cell area.
    """
    d = ssp_inputs(idata) if d is None else d
    theta = flat_samples(idata)
    rows = d.rows(theta[:, 0], theta[:, 1])

    ok = rows >= 0
    if not ok.all():
        _logger.warning(f'{int((~ok).sum())} of the {rows.size} posterior samples fall on (age, Z) combinations '
                        f'missing from the grid. They are excluded from the node table.')

    df = pd.DataFrame({'row': rows[ok], 'ebv': theta[ok, 2], 'log_A': theta[ok, 3]})
    tab = df.groupby('row').agg(n=('ebv', 'size'), ebv_p50=('ebv', 'median'), log_A_p50=('log_A', 'median'))

    weight = tab['n'].to_numpy(dtype=float)
    if uniform_node_prior:
        (t_lo, t_hi), (z_lo, z_hi) = d.bounds[0], d.bounds[1]
        i_age = np.searchsorted(d.age_log, d.age_node[tab.index])  # exact matches: age_log = unique(age_node)
        i_z = np.searchsorted(d.z_log, d.z_node[tab.index])
        weight /= _cell_widths(d.age_log, t_lo, t_hi)[i_age] * _cell_widths(d.z_log, z_lo, z_hi)[i_z]

    tab['prob'] = weight / weight.sum()
    tab.insert(0, 'log_age', d.age_node[tab.index])
    tab.insert(1, 'Z', np.power(10, d.z_node[tab.index]))

    return tab.sort_values('prob', ascending=False)


def predictive_band(idata, n_draws=500, q=(16, 50, 84), seed=None, d=None):
    """Pointwise percentiles of the model spectra of random posterior draws (marginalizes over nodes)."""
    d = ssp_inputs(idata) if d is None else d
    theta = flat_samples(idata)
    idx = np.random.default_rng(seed).choice(len(theta), size=min(n_draws, len(theta)), replace=False)
    models = np.array([m for m in (d.model(t) for t in theta[idx]) if m is not None])
    return np.percentile(models, q, axis=0)


def best_fit(idata, uniform_node_prior=False, d=None):

    """
    Self-consistent optimum without interpolation: the (age, Z) node with the highest posterior mass and, for
    that same node, the E(B-V) and log(A) that maximize the likelihood (refine=True) or their conditional
    posterior medians (refine=False). Returns theta, the model spectrum and the node table.
    """

    d = ssp_inputs(idata) if d is None else d
    tab = node_table(idata, uniform_node_prior, d=d)
    row = int(tab.index[0])

    ebv, log_A = tab.iloc[0][['ebv_p50', 'log_A_p50']]

    theta = [float(d.age_node[row]), float(d.z_node[row]), float(ebv), float(log_A)]

    med = np.median(flat_samples(idata)[:, :2], axis=0)
    row_med = int(d.rows(med[0], med[1])[0])
    if row_med != row:
        p_med = float(tab['prob'].get(row_med, 0.0))
        _logger.warning(f'The marginal-median node (log(age)={med[0]:.2f}, log10(Z)={med[1]:.3f}) differs from '
                        f'the maximum-posterior node (log(age)={theta[0]:.2f}, log10(Z)={theta[1]:.3f}): '
                        f'P = {p_med:.3f} vs {tab["prob"].iloc[0]:.3f}.')

    return theta, d.model(theta), tab


def truth_to_params(truth_dict):
    """Simulation truth (log_age, metallicity, ebv, mass [Msun]) -> {param name: value} in the sampled parameters."""
    return {'log_age': float(truth_dict['log_age']), 'log_Z': float(np.log10(truth_dict['metallicity'])),
            'ebv': float(truth_dict['ebv']), 'log_A': float(np.log10(truth_dict['mass'] / 1e6))}


def plot_best_fit(trace, truth=None, title=None, savefile=None):

    inputs_dict = ssp_inputs(trace)
    theta, best_model, band, nodes, stats = load_ssp_trace_results(trace, PARAM_NAMES)

    names = list(PARAM_NAMES)
    n = len(names)

    # True values: the argument takes precedence over those stored in the trace
    ref = None
    if truth is not None:
        ref = truth_to_params(truth)
    elif 'truth' in trace.fit_results:
        ref = dict(zip(names, trace.fit_results['truth'].values))

    # Initial values of the sampler (absent in traces saved before p0 was stored)
    p0 = None
    if 'p0' in trace.fit_results:
        p0 = dict(zip(names, trace.fit_results['p0'].values))

    # Bootstrap initial values (only present if the p0 bootstrap was requested)
    p0_arr = None
    if 'p0_arr' in trace.fit_results:
        p0_arr = dict(zip(names, trace.fit_results['p0_arr'].values.T))

    cd = trace.constant_data
    x, y, mask = cd['wave'].values, cd['lum'].values, cd['mask'].values.astype(bool)

    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], width_ratios=[2, 1.3])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax0.tick_params(labelbottom=False)

    # Shade masked-out regions (mask True = fitted)
    for lo, hi in np.flatnonzero(np.diff(np.r_[False, ~mask, False])).reshape(-1, 2):
        for ax in (ax0, ax1):
            ax.axvspan(x[lo], x[min(hi, x.size - 1)], alpha=0.5, color='lightgrey', lw=0)

    ax0.step(x, y, color='black', lw=1, label='Input spectrum', zorder=1)

    if band is not None:
        ax0.fill_between(x, band[0], band[2], alpha=0.25, color='teal', lw=0, label='Predictive 16-84%', zorder=0)

    label = f'Optimal Model (log(age)={theta[0]:.2f}, Z={10 ** theta[1]:.4f}, E(B-V)={theta[2]:.2f})'
    ax0.step(x, best_model, color='royalblue', lw=1, label=label, zorder=500)

    star, neb = ssp_components(theta, inputs_dict)
    if neb is not None:
        ax0.step(x, star, color='darkorange', lw=0.8, label='Stellar', zorder=400)
        ax0.step(x, neb, color='crimson', lw=0.8, label='Nebular continuum', zorder=400)

    ax0.set_ylabel(r"L$_{\odot}$ $\AA^{-1}$")
    if title is not None:
        ax0.set_title(title, weight='semibold')
    ax0.legend(loc='best')

    ax1.step(x, (y - best_model) / y, color='royalblue', lw=1)
    ax1.axhline(0, color='black')
    ax1.set_ylabel('Residuals')
    ax1.set_ylim(-0.7, 0.7)
    ax1.set_xlabel(r"Wavelength ($\AA$)")

    # Scatter matrix nested in the ax2 slot: histograms on the diagonal, scatter in the lower triangle

    # Marginal posterior percentiles (p16, p50, p84) of each parameter
    pct = {q: stats[q].to_dict() for q in ('p16', 'p50', 'p84')}
    samples = {k: trace.posterior[k].values.ravel() for k in names}

    sub = gs[0, 1].subgridspec(n, n, wspace=0.08, hspace=0.08)
    mat = np.empty((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            if j > i:
                continue
            ax = fig.add_subplot(sub[i, j], sharex=mat[0, j] if i > 0 else None)
            mat[i, j] = ax
            if i == j:
                ax.hist(samples[names[i]], bins=30, color='teal', alpha=0.7)
                ax.axvspan(pct['p16'][names[i]], pct['p84'][names[i]], color='gray', alpha=0.3)
                ax.axvline(pct['p50'][names[i]], color='black', lw=1.5)
                if p0 is not None:
                    ax.axvline(p0[names[i]], color='blue', ls='--', lw=1.5)
                if ref is not None:
                    ax.axvline(ref[names[i]], color='red', ls='-', lw=1.5)
                ax.set_yticks([])
            else:
                ax.scatter(samples[names[j]], samples[names[i]], s=3, alpha=0.1, color='teal', rasterized=True)
                if p0_arr is not None:
                    ax.scatter(p0_arr[names[j]], p0_arr[names[i]], s=6, alpha=0.3, color='blue', lw=0, zorder=1.5,
                               rasterized=True)
                ax.axvspan(pct['p16'][names[j]], pct['p84'][names[j]], color='gray', alpha=0.2)
                ax.axhspan(pct['p16'][names[i]], pct['p84'][names[i]], color='gray', alpha=0.2)
                ax.axvline(pct['p50'][names[j]], color='black', lw=0.8)
                ax.axhline(pct['p50'][names[i]], color='black', lw=0.8)
                if p0 is not None:
                    ax.axvline(p0[names[j]], color='blue', ls='--', lw=0.8)
                    ax.axhline(p0[names[i]], color='blue', ls='--', lw=0.8)
                if ref is not None:
                    ax.axvline(ref[names[j]], color='red', ls='-', lw=0.8)
                    ax.axhline(ref[names[i]], color='red', ls='-', lw=0.8)
                    ax.scatter(ref[names[j]], ref[names[i]], marker='*', s=120, color='red', edgecolor='white', zorder=11)
            # Labels only on the outer edges
            if i == n - 1:
                ax.set_xlabel(names[j])
            else:
                ax.tick_params(labelbottom=False)
            if j == 0 and i > 0:
                ax.set_ylabel(names[i])
            elif j > 0:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=7)

    # Row of trace plots (one per parameter, all chains overlaid) in the lower-right slot
    # (the empty top row of the subgrid leaves room for the x labels of the scatter matrix above)
    sub_tr = gs[1, 1].subgridspec(2, n, height_ratios=[0.5, 1], wspace=0.45, hspace=0)
    ax_tr = []
    for k, name in enumerate(names):
        ax = fig.add_subplot(sub_tr[1, k])
        chains = trace.posterior[name].values  # (chain, draw)
        ax.plot(np.arange(chains.shape[1]), chains.T, lw=0.3, alpha=0.25, color='teal', rasterized=True)
        ax.axhspan(pct['p16'][name], pct['p84'][name], color='gray', alpha=0.3)
        ax.axhline(pct['p50'][name], color='black', lw=1.5)
        if p0 is not None:
            ax.axhline(p0[name], color='blue', ls='--', lw=1.5)
        if ref is not None:
            ax.axhline(ref[name], color='red', ls='-', lw=1.5)
        ax.set_xlabel(name, fontsize=9)
        ax.tick_params(labelsize=7)
        ax_tr.append(ax)

    fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.09, wspace=0.15, hspace=0.1)

    if savefile:
        fig.savefig(savefile, bbox_inches='tight')
    else:
        plt.show()

    return fig, (ax0, ax1, mat, ax_tr)


def summary_ssp(trace, truth_dict=None, var_names=None):

    # Orignal traces to get from the
    df = az.summary(trace, var_names=var_names)

    # Insert p50, p16, p84 as the 3rd, 4th and 5th columns
    post = trace.posterior
    for i, (col, q) in enumerate(zip(('p50', 'p16', 'p84'), (50, 16, 84)), start=2):
        df.insert(i, col, [np.percentile(post[name].values, q) for name in df.index])

    # Insert the true values from a simulation: the argument takes precedence over those stored in the trace
    truth_vals = None
    if truth_dict is not None:
        truth_vals = truth_to_params(truth_dict)
    elif hasattr(trace, 'fit_results') and 'truth' in trace.fit_results:
        truth_vals = dict(zip(PARAM_NAMES, trace.fit_results['truth'].values))

    if truth_vals is not None:
        df.insert(5, 'true_values', [truth_vals[name] for name in df.index])

    return df


def single_pop_model(bounds, lum, lum_err, mask, ssp_sampler):

    (t_lo, t_hi), (z_lo, z_hi), (e_lo, e_hi), (a_lo, a_hi) = bounds

    @as_op(itypes=[pt.dscalar] * 4, otypes=[pt.dscalar])
    def _loglike(log_age, log_met, ebv, log_amp):
        ll = log_likelihood([float(log_age), float(log_met), float(ebv), float(log_amp)],
                            lum, lum_err, mask, ssp_sampler)
        return np.asarray(ll, dtype='float64')

    with pm.Model() as model:

        # Flat priors in the native parameter space (no interval transform): hard walls as in emcee
        log_age = pm.Uniform('log_age', t_lo, t_hi, default_transform=None)
        log_z = pm.Uniform('log_Z', z_lo, z_hi, default_transform=None)
        ebv = pm.Uniform('ebv', e_lo, e_hi, default_transform=None)
        log_amp = pm.Uniform('log_A', a_lo, a_hi, default_transform=None)

        pm.Potential('loglike', _loglike(log_age, log_z, ebv, log_amp))

    return model


def chain_chi2(idata, d=None):
    """chi2 of every posterior sample (nearest-node model), shape (chain, draw)."""

    d = ssp_inputs(idata) if d is None else d
    theta = np.stack([idata.posterior[name].values for name in PARAM_NAMES], axis=-1)    # (chain, draw, param)
    flat = theta.reshape(-1, len(PARAM_NAMES))
    rows = d.rows(flat[:, 0], flat[:, 1])

    good = d.mask & np.isfinite(d.lum) & (d.lum_err > 0)
    y, err, flux, red = d.lum[good], d.lum_err[good], d.grid_flux[:, good], d.red_corr[good]

    chi2 = np.full(len(flat), np.inf)
    for i, (row, (_, _, ebv, log_amp)) in enumerate(zip(rows, flat)):
        if row >= 0:
            chi2[i] = np.sum(((y - 10.0 ** log_amp * flux[row] * red ** ebv) / err) ** 2)

    return chi2.reshape(theta.shape[:2])


def stuck_chains(idata, k=1.5, dchi2=20.0):
    """
    Boolean mask of the chains to keep. A chain is dropped only if its median chi2 is both more than `dchi2` above
    the best chain and above the outlier fence Q3 + k * IQR of the chain medians.
    """

    med = np.median(chain_chi2(idata), axis=1)
    q1, q3 = np.percentile(med, (25, 75))
    keep = (med - med.min() <= dchi2) | (med <= q3 + k * (q3 - q1))

    if not keep.all():
        _logger.warning(f'{int((~keep).sum())} of the {keep.size} chains dropped as stuck: median chi2 above the '
                        f'best chain by {np.round(np.sort(med[~keep] - med.min()), 1).tolist()}.')

    return keep


def update_fit_results(idata, band_draws=1000, seed=None):
    """Recomputes the fit_results group (best fit, node table, predictive band, marginal stats) from the posterior."""

    d = ssp_inputs(idata)
    old = idata.fit_results

    theta_best, best_model, tab = best_fit(idata, d=d)
    p16, p50, p84 = predictive_band(idata, n_draws=band_draws, seed=seed, d=d)
    stats = marginal_stats(idata)

    idata['fit_results'] = xr.Dataset(
        {'theta': ('param', np.asarray(theta_best, dtype=float)),
         **{key: (old[key].dims, old[key].values) for key in ('p0', 'p0_arr', 'truth') if key in old},
         'best_model': ('pixel', best_model),
         'band_p16': ('pixel', p16), 'band_p50': ('pixel', p50), 'band_p84': ('pixel', p84),
         'node_log_age': ('rank', tab['log_age'].to_numpy()), 'node_Z': ('rank', tab['Z'].to_numpy()),
         'node_ebv_p50': ('rank', tab['ebv_p50'].to_numpy()), 'node_log_A_p50': ('rank', tab['log_A_p50'].to_numpy()),
         'node_prob': ('rank', tab['prob'].to_numpy()),
         'marginal_p16': ('param', stats['p16'].to_numpy()), 'marginal_p50': ('param', stats['p50'].to_numpy()),
         'marginal_p84': ('param', stats['p84'].to_numpy()),
         **({'marginal_ess': ('param', stats['ess'].to_numpy()),
             'marginal_tau': ('param', stats['tau'].to_numpy())} if 'ess' in stats else {})},
        coords={'param': list(PARAM_NAMES)})

    return idata




class SSP_sampler:

    def __init__(self, spec, stellar_binaries, red_law='CCM89', r_v=3.1, **params):

        # Inputs
        self.spec = spec
        self.stellar_ssp = stellar_binaries

        # Extinction
        self.red_law = red_law
        self.r_v = DEFAULT_RV.get(red_law, 3.1) if r_v is None else r_v
        self.red_corr = None

        # Model variables
        self.z_log = None
        self.age_log = None
        self.grid_theo = None
        self.node_idx = None
        self.grid = None

        # Nebular continuum
        self.add_nebular = False
        self.grid_stellar = None
        self.grid_neb = None

        # Observation
        self.wave_arr = None
        self.lum = None
        self.lum_err = None
        self.mask = None
        self.mask = None

        self.age_node = None
        self.z_node = None
        self.bounds = None
        self.p0 = None
        self.p0_arr = None

        # Sampling
        self.pm_model = None
        self.label_model = None
        self.idata = None

        return

    def security_checks(self, wave_arr):

        # Check SSP wavelength uniformity
        if not self.stellar_ssp.uniform_dispersion:
            raise SpecSyError(f'The input SSP data set does not have a unique wavelength array')

        # Check SSP wavelength matches observation
        if not np.all(np.isclose(self.stellar_ssp.wave_series.at[0], wave_arr)):
            raise SpecSyError(f'The input SSP wavelength array does not match the observation wavelength array')

        return


    def prepare_inputs(self, wave_arr, flux_arr, err_arr, mask_arr, model=None, prior_bounds=None, initial_values=None,
                       add_nebular=False, ion_frame=None, gamma_table=None, neb_grid=None, p0_bootstrap=False,
                       rng_seed=None, **params):

        # Check the object spectrum and SSP match
        self.security_checks(wave_arr)

        # Check model
        self.label_model = 'single_population_nearest' if model is None else model

        # Prepare the inputs for the model
        match self.label_model:
            case 'single_population_nearest':
                subset =  self.stellar_ssp.frame.loc[self.stellar_ssp._select(**params)]
                pairs = subset[['age', 'metallicity']].drop_duplicates().sort_values(['metallicity', 'age'])

                self.age_node = pairs['age'].to_numpy(dtype=float)
                self.z_node = np.log10(pairs['metallicity'].to_numpy(dtype=float))

                # Grid axes and flux row of each (metallicity, age) node (-1 for combinations missing in the grid)
                self.age_log, i_age = np.unique(self.age_node, return_inverse=True)
                self.z_log, i_z = np.unique(self.z_node, return_inverse=True)

                self.node_idx = np.full((self.z_log.size, self.age_log.size), -1, dtype=int)
                self.node_idx[i_z, i_age] = np.arange(self.age_node.size)

                # Compute the Theoretical fluxes for the model
                self.grid_theo = np.stack(self.stellar_ssp.flux_series.loc[pairs.index].to_numpy())

                # Nebular continuum: fixed per node (set by its Q), so it is added to the grid once
                self.grid_stellar, self.grid_neb, self.add_nebular = self.grid_theo, None, add_nebular
                if add_nebular:
                    if neb_grid is None:
                        if ion_frame is None:
                            raise SpecSyError('add_nebular=True requires ion_table (from load_ionization_table) '
                                              'or a precomputed neb_grid')
                        if isinstance(ion_frame, (str, bytes)) or hasattr(ion_frame, '__fspath__'):
                            ion_frame = load_ionization_frame(ion_frame)
                        log_q = ion_table_log_q(ion_frame, self.age_node, self.z_node)
                        gamma_table = (NEB_GAMMA_WAVE, NEB_GAMMA_COEF) if gamma_table is None else gamma_table
                        neb_grid = nebular_continuum(wave_arr, log_q, *gamma_table)

                    self.grid_neb = np.asarray(neb_grid, dtype=float)
                    if self.grid_neb.shape != self.grid_stellar.shape:
                        raise SpecSyError(f'Nebular grid shape {self.grid_neb.shape} does not match the SSP grid '
                                          f'{self.grid_stellar.shape}')
                    self.grid_theo = self.grid_stellar + self.grid_neb

                # Observation and fitted pixels (unmasked, finite flux and positive uncertainty)
                self.wave_arr, self.lum, self.lum_err, self.mask = wave_arr, flux_arr, err_arr, mask_arr

                # Get the extinciton curve
                self.red_corr = np.power(10.0, -0.4 * extinction_curve(wave_arr, self.red_law, self.r_v))

                # Unpack the model priors
                if prior_bounds is None:
                    prior_bounds = np.array(specsy_cfg['stellar']['ssp'][self.stellar_ssp.source]['bounds'])
                self.bounds = np.asarray(prior_bounds, dtype=float)

                # Full-grid container used by the model evaluation (sampling) and, reduced to the visited nodes, by
                # the stored trace (grid_flux is the total stellar + nebular flux, as in grid_theo)
                self.grid = SSPGrid(wave=self.wave_arr, lum=self.lum, lum_err=self.lum_err, mask=self.mask,
                                    red_corr=self.red_corr, grid_flux=self.grid_theo, age_node=self.age_node,
                                    z_node=self.z_node, age_log=self.age_log, z_log=self.z_log, bounds=self.bounds,
                                    grid_neb=self.grid_neb if self.add_nebular else None,
                                    red_law=self.red_law, r_v=self.r_v)

                # Get the initial coordinate for the fitting: a single fit to the observed flux or, if requested, the
                # nanmean of the fits to p0_bootstrap flux realizations drawn from the uncertainties (True = 500)
                n_boot = 500 if p0_bootstrap is True else int(p0_bootstrap or 0)
                self.p0_arr = None
                if initial_values is not None:
                    if n_boot > 0:
                        _logger.warning(f'Initial values were provided: the p0 bootstrap ({n_boot} realizations) '
                                        f'is skipped.')
                    self.p0 = initial_values
                elif n_boot > 0:
                    flux_in, err_in = np.asarray(flux_arr, dtype=float), np.asarray(err_arr, dtype=float)
                    noise = np.random.default_rng(rng_seed).standard_normal((n_boot, flux_in.size))
                    self.p0_arr = get_initial_values(flux_in + err_in * noise, err_in, self.mask, self, self.bounds)
                    self.p0 = np.nanmean(self.p0_arr, axis=0)
                else:
                    self.p0 = get_initial_values(flux_arr, err_arr, self.mask, self, self.bounds, show_solution=False)

                p0_label = '' if self.p0_arr is None else f' (nanmean of {len(self.p0_arr)} bootstrap fits)'
                print(f'- p0 values{p0_label}: log(age) = {self.p0[0]:0.4f}, '
                                    f'Z = {np.power(10, self.p0[1]):0.5f}, '
                                    f'E(B-V) = {self.p0[2]:0.3f}, '
                                    f'log(A) = {self.p0[3]:0.4f}')

                # Declaring the model
                self.pm_model = single_pop_model(self.bounds, self.lum, self.lum_err, self.mask, self)

            case _ :
                raise SpecSyError(f'SSP synthesis model "{model}" is not recognized')

        return

    def sample(self, draws=5000, tune=1000, n_walkers=64, cores=None, sampler='DEMetropolis', spread=0.10,
               spread_step=True, step_kwargs=None, rng_seed=None, band_draws=1000, truth_dict=None, **kwargs):

        if sampler not in MCMC_SAMPLERS:
            raise SpecSyError(f'MCMC sampler "{sampler}" is not recognized. Available samplers: {MCMC_SAMPLERS}')

        # Initial values
        initvals, p0_walkers = sample_truncnorm(self.p0, self.bounds[:, 0], self.bounds[:, 1], n=n_walkers,
                                                stellar_binaries=self.stellar_ssp, spread=spread, seed=rng_seed)

        # Cores mechanic
        if cores is None:
            cores = 1 if sampler in ('DEMetropolis', 'DEMetropolisZ') else 4

        with self.pm_model:
            # with warnings.catch_warnings():
            #     warnings.filterwarnings('ignore', message='overflow encountered in exp', category=RuntimeWarning)

                #
                # if quiet_overflow:
                #

            # Define the sampler and its step (scale from the shifted p0, floor for parameters at zero)
            if 'step' in kwargs:
                if step_kwargs:
                    _logger.warning(f'A step object was provided: the step_kwargs {list(step_kwargs)} are ignored.')
            else:
                step_args = dict(step_kwargs or {})
                if spread_step and 'S' not in step_args:
                    S = np.abs(p0_walkers) * spread
                    S[S == 0] = spread
                    step_args['S'] = S
                kwargs['step'] = getattr(pm, sampler)(**step_args)

            self.idata = pm.sample(draws=draws, tune=tune, chains=n_walkers, cores=cores, initvals=initvals,
                                   random_seed=rng_seed, **kwargs)

        # Add the inputs to the trace for future reconstruction
        self._parse_results(PARAM_NAMES, band_draws, rng_seed, truth_dict)

        return self.idata

    def summary(self, truth_dict=None, var_names=None, verbose=True):

        # Get the results table
        results_df = summary_ssp(self.idata, truth_dict, var_names=var_names)

        # Print the true results versus output
        if verbose:
            if truth_dict is not None:
                print(f' - True values: log(age) = {results_df.loc["log_age", 'true_values']:0.3f}, '
                                        f'Z = {results_df.loc["log_Z", 'true_values']:0.5f}, '
                                        f'E(B-V) = {results_df.loc["ebv", 'true_values']:0.3f}, '
                                        f'log(A) = {results_df.loc["log_A", 'true_values']:0.5f}')

                print(f' - Initial values: log(age) = {self.p0[0]:0.3f}, '
                                        f'Z = {self.p0[1]:0.5f}, '
                                        f'E(B-V) = {self.p0[2]:0.3f}, '
                                        f'log(A) = {self.p0[3]:0.5f}')

                print(f' - Fitted values: log(age) = {results_df.loc["log_age", 'p50']:0.3f}, '
                                        f'Z = {results_df.loc["log_Z", 'p50']:0.5f}, '
                                        f'E(B-V) = {results_df.loc["ebv", 'p50']:0.3f}, '
                                        f'log(A) = {results_df.loc["log_A", 'p50']:0.5f}')

        return results_df

    def nearest(self, log_age, log_met):
        """Flux row of the grid node nearest to (log_age, log_met), -1 if the combination is missing."""
        return int(self.grid.rows(log_age, log_met)[0])

    def model(self, theta):
        return self.grid.model(theta)

    def _parse_results(self, param_names, band_draws, random_seed, truth_dict=None):

        # Posterior samples as an (n_samples, 4) array with columns in param_names
        post = self.idata.posterior
        theta = np.column_stack([post[name].values.ravel() for name in param_names])

        # Nearest grid node of every sample (ties go to the upper node); -1 marks (age, Z) cells missing in the grid
        nearest = []
        for axis, x in ((self.z_log, theta[:, 1]), (self.age_log, theta[:, 0])):
            if axis.size == 1:
                nearest.append(np.zeros(x.size, dtype=int))
                continue
            i = np.clip(np.searchsorted(axis, x, side='left'), 1, axis.size - 1)
            nearest.append(i - ((x - axis[i - 1]) < (axis[i] - x)).astype(int))
        rows = self.node_idx[nearest[0], nearest[1]]

        # Inputs: observation, extinction, bounds, full axes and only the grid nodes visited by the posterior
        visited = np.unique(rows[rows >= 0])
        inputs = xr.Dataset({'wave': ('pixel', self.wave_arr),
                                       'lum': ('pixel', self.lum),
                                       'lum_err': ('pixel', self.lum_err),
                                       'mask': ('pixel', self.mask),
                                       'red_corr': ('pixel', self.red_corr),
                                       'grid_flux': (('node', 'pixel'), self.grid_theo[visited]),
                                       'age_node': ('node', self.age_node[visited]),
                                       'z_node': ('node', self.z_node[visited]),
                                       'age_axis': ('age_grid', self.age_log),
                                       'z_axis': ('z_grid', self.z_log),
                                       'bounds': (('param', 'limit'), self.bounds),
                                       **({'grid_neb': (('node', 'pixel'), self.grid_neb[visited])} if self.add_nebular else {})},
                            coords={'param': list(PARAM_NAMES), 'limit': ['lower', 'upper']},
                            attrs={'red_law': str(self.red_law),
                                   'r_v': float(self.r_v),
                                   'units_wave': str(self.spec.units_wave),
                                   'units_flux':  str(self.spec.units_flux),
                                   'add_nebular': int(self.add_nebular)})

        self.idata['inputs'] = inputs

        # Best fit, node table, predictive band and marginal-posterior stats
        params_fit, spectrum_fit, node_table = best_fit(self.idata, d=self.grid)
        p16, p50, p84 = predictive_band(self.idata, n_draws=band_draws, seed=random_seed, d=self.grid)
        stats = marginal_stats(self.idata)

        # Truth in PARAM_NAMES order
        if truth_dict is not None:
            truth = np.array([truth_dict['log_age'], np.log10(truth_dict['metallicity']),
                              truth_dict['ebv'],     np.log10(truth_dict['mass'] / 1e6)], dtype=float)

        outputs = xr.Dataset({'theta': ('param', np.asarray(params_fit, dtype=float)),
                              **({'p0': ('param', np.asarray(self.p0, dtype=float))} if self.p0 is not None else {}),
                              **({'p0_arr': (('p0_draw', 'param'), np.asarray(self.p0_arr, dtype=float))}
                                 if self.p0_arr is not None else {}),
                              **({'truth': ('param', truth)} if truth_dict is not None else {}),
                              'best_model': ('pixel', spectrum_fit),
                              'band_p16': ('pixel', p16), 'band_p50': ('pixel', p50), 'band_p84': ('pixel', p84),
                              'node_log_age': ('rank', node_table['log_age'].to_numpy()),
                              'node_Z': ('rank', node_table['Z'].to_numpy()),
                              'node_ebv_p50': ('rank', node_table['ebv_p50'].to_numpy()),
                              'node_log_A_p50': ('rank', node_table['log_A_p50'].to_numpy()),
                              'node_prob': ('rank', node_table['prob'].to_numpy()),
                              'marginal_p16': ('param', stats['p16'].to_numpy()),
                              'marginal_p50': ('param', stats['p50'].to_numpy()),
                              'marginal_p84': ('param', stats['p84'].to_numpy()),
                              **({'marginal_ess': ('param', stats['ess'].to_numpy()),
                                  'marginal_tau': ('param', stats['tau'].to_numpy())} if 'ess' in stats else {})},
                             coords={'param': list(PARAM_NAMES)},
                             attrs={'param_names': list(PARAM_NAMES), 'model': self.label_model})
        self.idata['outputs'] = outputs

        return

    def save_trace(self, fname, trace=None):

        if trace == None:
            trace = self.idata

        trace.to_netcdf(fname)

        return


def log_likelihood(theta, y, yerr, mask, nodes):
    """Gaussian log-likelihood of the node nearest to theta (drop-in for the MCMC stage)."""

    y_model = nodes.model(theta)
    if y_model is None:
        return -np.inf

    resid = (y[mask] - y_model[mask]) / yerr[mask]

    return -0.5 * (resid @ resid + np.sum(np.log(2 * np.pi * yerr[mask] ** 2)))


def get_initial_values(lum, lum_err, mask, nodes, bounds, ebv_step=0.005, show_solution=False):

    """
    Maximum-likelihood (log age, log10 Z, E(B-V), log A) inside `bounds` by exhaustive search: every grid node,
    E(B-V) on a regular grid no coarser than `ebv_step`, and the amplitude solved analytically by weighted linear
    least squares (clipped to its bounds). The likelihood normalization term does not depend on theta, so
    minimizing chi2 is equivalent to minimizing the negative log-likelihood.
    `lum` is either a single spectrum (n_pixel,), returning [log_age, log_Z, ebv, log_A], or a stack of flux
    realizations (n_real, n_pixel), returning an (n_real, 4) array with one solution per row.
    """

    (t_lo, t_hi), (z_lo, z_hi), (e_lo, e_hi), (a_lo, a_hi) = bounds
    mask = np.asarray(mask, dtype=bool)
    lum = np.asarray(lum, dtype=float)
    single = lum.ndim == 1
    lum = np.atleast_2d(lum)

    # Fitted pixels and their weights
    good = mask & np.isfinite(lum).all(axis=0) & np.isfinite(lum_err) & (lum_err > 0)
    n_bad = int(mask.sum() - good.sum())
    if n_bad > 0:
        _logger.warning(f'{n_bad} of the {int(mask.sum())} unmasked pixels have a non-finite flux or a non-positive '
                        f'uncertainty. They are excluded from the fit.')
    y_arr, w = lum[:, good], 1.0 / lum_err[good] ** 2

    # Grid nodes inside the age and metallicity bounds
    in_box = ((nodes.age_node >= t_lo) & (nodes.age_node <= t_hi)
              & (nodes.z_node >= z_lo) & (nodes.z_node <= z_hi))
    if not in_box.any():
        raise SpecSyError(f'No grid nodes inside the bounds log(age)={bounds[0]}, log10(Z)={bounds[1]}. '
                          f'Available: log(age) {nodes.age_log[0]}-{nodes.age_log[-1]}, '
                          f'log10(Z) {nodes.z_log[0]:.3f}-{nodes.z_log[-1]:.3f}.')
    idcs_box = np.flatnonzero(in_box)
    flux = nodes.grid_theo[idcs_box][:, good]

    n_nan = int((~np.isfinite(flux).all(axis=1)).sum())
    if n_nan > 0:
        _logger.warning(f'{n_nan} of the {idcs_box.size} grid nodes have non-finite fluxes in the fitted pixels. '
                        f'They are excluded from the search.')

    # Transmission on the E(B-V) grid and the model-model term, shared by all the flux realizations
    ebv_grid = np.linspace(e_lo, e_hi, int(np.ceil((e_hi - e_lo) / ebv_step)) + 1)
    trans = nodes.red_corr[good] ** ebv_grid[:, None]
    w_trans = w * trans
    s_mm = (w_trans * trans) @ (flux ** 2).T

    # chi2(A) = s_yy - 2 A s_ym + A^2 s_mm for every (E(B-V), node) pair of each realization
    theta, chi2_min = np.full((y_arr.shape[0], 4), np.nan), np.full(y_arr.shape[0], np.nan)
    for r, y in enumerate(y_arr):
        s_ym = (y * w_trans) @ flux.T
        with np.errstate(divide='ignore', invalid='ignore'):
            amp = np.clip(s_ym / s_mm, 10.0 ** a_lo, 10.0 ** a_hi)
        chi2 = np.sum(w * y ** 2) - 2 * amp * s_ym + amp ** 2 * s_mm

        i_ebv, i_node = np.unravel_index(np.nanargmin(chi2), chi2.shape)
        node = idcs_box[i_node]
        theta[r] = nodes.age_node[node], nodes.z_node[node], ebv_grid[i_ebv], np.log10(amp[i_ebv, i_node])
        chi2_min[r] = chi2[i_ebv, i_node]

    # Solutions on the edge of the search space
    ranges = (('log(age)', nodes.age_node[in_box]), ('log10(Z)', nodes.z_node[in_box]),
              ('E(B-V)', ebv_grid), ('log(A)', np.array([a_lo, a_hi])))
    on_edge = np.column_stack([(rng.min() < rng.max()) & (np.isclose(theta[:, k], rng.min())
                                                          | np.isclose(theta[:, k], rng.max()))
                               for k, (_, rng) in enumerate(ranges)])
    if on_edge.any():
        if single:
            at_edge = [f'{name}={theta[0, k]:.3f}' for k, (name, _) in enumerate(ranges) if on_edge[0, k]]
            _logger.warning(f'The initial solution is on the edge of the search space for: {at_edge}. '
                            f'Consider widening the bounds.')
        else:
            at_edge = {name: int(n) for (name, _), n in zip(ranges, on_edge.sum(axis=0)) if n > 0}
            _logger.warning(f'{int(on_edge.any(axis=1).sum())} of the {len(theta)} bootstrap solutions are on the '
                            f'edge of the search space (solutions per parameter: {at_edge}). '
                            f'Consider widening the bounds.')

    if show_solution:
        n_dof = int(good.sum()) - theta.shape[1]
        if single:
            print(f'log(age) = {theta[0, 0]:.2f}, log10(Z) = {theta[0, 1]:.3f} (Z = {10 ** theta[0, 1]:.4f}), '
                  f'E(B-V) = {theta[0, 2]:.3f}, log(A) = {theta[0, 3]:.3f}\n'
                  f'chi2 = {chi2_min[0]:.1f}, reduced chi2 = {chi2_min[0] / n_dof:.2f} ({n_dof} dof), '
                  f'searched {idcs_box.size} nodes x {ebv_grid.size} E(B-V) values')
        else:
            mu, sd = np.nanmean(theta, axis=0), np.nanstd(theta, axis=0)
            print(f'{len(theta)} realizations: log(age) = {mu[0]:.2f} +/- {sd[0]:.2f}, '
                  f'log10(Z) = {mu[1]:.3f} +/- {sd[1]:.3f}, E(B-V) = {mu[2]:.3f} +/- {sd[2]:.3f}, '
                  f'log(A) = {mu[3]:.3f} +/- {sd[3]:.3f}\n'
                  f'median reduced chi2 = {np.nanmedian(chi2_min) / n_dof:.2f} ({n_dof} dof), '
                  f'searched {idcs_box.size} nodes x {ebv_grid.size} E(B-V) values')

    return [float(v) for v in theta[0]] if single else theta