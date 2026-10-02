from pathlib import Path
from time import perf_counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pymc as pm
import arviz as az
import lime
from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import (SSP_sampler, plot_best_fit, load_ionization_frame, summary_ssp,
                                             truth_to_params)

# Stellar grid properties
bpass_pname = Path('/home/vital/PycharmProjects/specsy/examples/stellar_synthesis/SESAMME_BAPSS-AP_v3.fits')
bpass_grid = StellarBinaries.from_fits(bpass_pname)

# Nebular continuum
ion_pname = Path('/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt')
ion_frame = load_ionization_frame(ion_pname)

# Configuration and results folder (mock spectra are read from here, each sampler configuration gets a subfolder)
synth_file = Path('./synth_sample.toml')
cfg = lime.load_cfg(synth_file)
results_folder = Path('/home/vital/Astrodata/SpecSy_tests/SSP_sampling')
tests_folder = results_folder / 'sampler_tests_2'

plot_fit = True     # best-fit figure per observation and configuration
overwrite = False   # False: traces already on disk are reloaded, only the missing configurations are sampled
rng_seed = 1234     # same seed for every configuration

# Second sweep, as overrides of the previous winner (DEM_w128). Budget = n_walkers x (draws + tune) evaluations:
# 64k unless marked x2 (128k). Step arguments are those of pm.DEMetropolis: lamb (jump factor, default 2.38/sqrt(2*4)
# = 0.84), scaling (size of the small random noise added to each jump, default 0.001), tune ('lambda' adapts lamb,
# 'scaling' adapts the noise, None = no adaptation) and tune_interval (steps between adaptations, default 100)
baseline = dict(sampler='DEMetropolis', n_walkers=128, draws=125, tune=375, spread=0.10, step={})
pymc_cfgs = {
    # Reference and budget control
    'W128_ref':       {},                                                   # previous winner
    'W128_x2':        dict(draws=250, tune=750),                            # x2: does more length help at all?

    # More walkers (the one change that helped)
    'W256':           dict(n_walkers=256, draws=62, tune=188),
    'W256_x2':        dict(n_walkers=256, draws=125, tune=375),             # x2
    'W512_x2':        dict(n_walkers=512, draws=62, tune=188),              # x2

    # Noise added to the DE jumps: lets chains leave a node when the difference vectors collapse onto it
    'W128_scal01':    dict(step=dict(scaling=0.01)),
    'W128_scal10':    dict(step=dict(scaling=0.1)),

    # Jump factor
    'W128_lamb06':    dict(step=dict(lamb=0.6)),
    'W128_lamb12':    dict(step=dict(lamb=1.2)),

    # Over-dispersed starting cloud (a more honest r_hat, and more nodes visited at the start)
    'W128_spread30':  dict(spread=0.30),

    # Lambda tuning (the DEM_tunelamb idea) with enough walkers and tuning steps to adapt
    'TL_w128':        dict(step=dict(tune='lambda')),
    'TL_w128_t875':   dict(draws=125, tune=875, step=dict(tune='lambda')),  # x2: longer tuning phase
    'TL_w128_int50':  dict(step=dict(tune='lambda', tune_interval=50)),     # faster adaptation
    'TL_w256_x2':     dict(n_walkers=256, draws=125, tune=375, step=dict(tune='lambda')),  # x2
}


def make_step(sampler, method, spread, **step_kwargs):
    """Step method on the sampler's model, with the S scale of SSP_sampler.sample (|p0| x spread, floored)."""
    S = np.abs(np.asarray(sampler.p0, dtype=float)) * spread
    S[S == 0] = spread
    with sampler.pm_model:
        return getattr(pm, method)(S=S, **step_kwargs)


def fit_metrics(trace, synth_dict):
    """Per-parameter accuracy and convergence table, plus the run-level scores used to rank the configurations."""
    df = summary_ssp(trace, truth_dict=synth_dict)
    half_width = 0.5 * (df['p84'] - df['p16']).replace(0, np.nan)
    df['z_err'] = (df['p50'] - df['true_values']) / half_width
    df['covered'] = df['true_values'].between(df['p16'], df['p84'])

    # Posterior mass on the grid node nearest to the truth (the mock truths are not exactly on grid nodes)
    truth, fr, cd = truth_to_params(synth_dict), trace.fit_results, trace.constant_data
    age_ax, z_ax = cd['age_axis'].values, cd['z_axis'].values
    age_true, z_true = age_ax[np.abs(age_ax - truth['log_age']).argmin()], z_ax[np.abs(z_ax - truth['log_Z']).argmin()]
    on_truth = (np.isclose(fr['node_log_age'].values, age_true, atol=1e-3)
                & np.isclose(np.log10(fr['node_Z'].values), z_true, atol=1e-3))

    stats = trace.sample_stats
    run = {'wall_s': trace.posterior.attrs.get('wall_time', np.nan),
           'accept_rate': float(stats['accepted'].mean()) if 'accepted' in stats else np.nan,
           'max_r_hat': df['r_hat'].max(), 'min_ess_bulk': df['ess_bulk'].min(),
           'mean_abs_z': df['z_err'].abs().mean(), 'coverage': df['covered'].mean(),
           'top_node_true': bool(on_truth[0]), 'p_true_node': float(fr['node_prob'].values[on_truth].sum())}
    run['min_ess_per_s'] = run['min_ess_bulk'] / run['wall_s']

    return df, run


# Fit all the mock observations with every configuration
run_rows, param_rows = [], []
for synth_spec, synth_dict in cfg.items():

    print(f'\n{synth_spec}')
    add_nebular = synth_dict.get('add_nebular', False)

    # Load the mock spectrum
    spec = lime.Spectrum.from_file(results_folder / f'{synth_spec}.txt', instrument='text')
    spec.unit_conversion(distance=synth_dict['distance'])
    wave_arr, flux_arr, err_arr, mask_arr = spec.retrieve.spectrum(return_arrays=True)
    mask_arr = ~mask_arr

    # Grid nodes, p0 and model declaration: computed once and shared by all the configurations
    obj_ssp = bpass_grid.to_stellar_binaries(wmin=wave_arr[0], wmax=wave_arr[-1])
    sampler = SSP_sampler(spec, obj_ssp, red_law='CCM', r_v=3.1)
    sampler.prepare_inputs(wave_arr, flux_arr, err_arr, mask_arr, add_nebular=add_nebular, ion_frame=ion_frame)

    for label, overrides in pymc_cfgs.items():

        conf = {**baseline, **overrides}
        step_kwargs = conf.pop('step')
        cfg_folder = tests_folder / label
        cfg_folder.mkdir(parents=True, exist_ok=True)
        trace_pname = cfg_folder / f'{synth_spec}_trace'

        # Sampling (the default step from SSP_sampler.sample unless step arguments are given)
        if overwrite or not trace_pname.exists():
            print(f'- {label}: {conf} {step_kwargs}')
            extra = {'step': make_step(sampler, conf['sampler'], conf['spread'], **step_kwargs)} if step_kwargs else {}
            t0 = perf_counter()
            sampler.sample(**conf, rng_seed=rng_seed, truth_dict=synth_dict, **extra)
            sampler.idata.posterior.attrs['wall_time'] = perf_counter() - t0
            sampler.save_trace(trace_pname)
        else:
            print(f'- {label}: loaded from {trace_pname}')

        # Reload the results and score the configuration
        trace_data = az.from_netcdf(trace_pname)
        results_df, run = fit_metrics(trace_data, synth_dict)
        results_df.to_csv(cfg_folder / f'{synth_spec}_summary.csv')

        print(f'  {run["wall_s"] / 60:.2f} min')

        n_eval = conf['n_walkers'] * (conf['draws'] + conf['tune'])
        run_rows.append({'spec': synth_spec, 'config': label, 'n_eval': n_eval, **run})
        param_rows.append(results_df.assign(spec=synth_spec, config=label))

        if plot_fit:
            plot_best_fit(trace_data, title=f'{synth_spec} ({label})',
                          savefile=cfg_folder / f'{synth_spec}_diagnostic.png')
            plt.close('all')

# Comparison tables: every run, every parameter and the configurations ranked over all the mocks
runs_df = pd.DataFrame(run_rows)
params_df = pd.concat(param_rows).rename_axis('param').reset_index()
runs_df.to_csv(tests_folder / 'runs_comparison.csv', index=False)
params_df.to_csv(tests_folder / 'params_comparison.csv', index=False)

ranking = (runs_df.groupby('config', sort=False)
           .agg(n_eval=('n_eval', 'first'), total_min=('wall_s', lambda t: t.sum() / 60),
                median_min=('wall_s', lambda t: t.median() / 60), accept_rate=('accept_rate', 'mean'),
                max_r_hat=('max_r_hat', 'max'), min_ess_per_s=('min_ess_per_s', 'min'),
                mean_abs_z=('mean_abs_z', 'mean'), coverage=('coverage', 'mean'),
                top_node_true=('top_node_true', 'mean'), p_true_node=('p_true_node', 'mean'))
           .sort_values(['top_node_true', 'p_true_node', 'min_ess_per_s'], ascending=False))
ranking.to_csv(tests_folder / 'config_ranking.csv')

with pd.option_context('display.width', 200, 'display.max_columns', None):
    print(f'\n{ranking.round(3)}')