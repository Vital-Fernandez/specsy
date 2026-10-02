from pathlib import Path
import numpy as np
import pandas as pd

from sesamme import models, vis
import sesamme.mcmc as stats
from sesamme.models import load_ionization_table, load_ssp_cube, get_model

import lime
import emcee


from astropy.table import Table
from astropy import units as u
from specsy.models.ssp import StellarBinaries
from scipy.optimize import minimize
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


def unpack_spectrum(fname, distance_mpc, luminosity_sun=3.828e33, MPC_CM = 3.0857e24):

    spec = lime.Spectrum.from_file(fname, instrument='text', norm_flux=1)
    # spec.plot.spectrum()
    wl = spec.wave.data
    flux = spec.flux.data
    flux_err = spec.err_flux.data
    mask = ~spec.flux.mask

    distance_cm = (distance_mpc * u.Mpc).to(u.cm).value
    conversion = 4 * np.pi * distance_cm**2 / luminosity_sun
    lum = flux * conversion
    lum_err = flux_err * conversion

    return wl, lum, lum_err, mask


def log_likelihood_model(theta, wl, y, yerr, mask, model_cube, ion_table, add_nebular):
    y_model = models.get_model(theta, wl, model_cube, ion_table, add_nebular)

    masked_spec = np.array((y[mask] - y_model[mask]) / yerr[mask])
    masked_err = np.array(np.sqrt(2) * np.sqrt(np.pi) * yerr[mask])

    resid = -0.5 * (np.dot(masked_spec, masked_spec) + np.log(np.dot(masked_err, masked_err)))

    return resid


def get_initial_values(initial, wl, lum, lum_err, mask, modelcube, bounds, iontable, add_nebular, show_solution=True):

    # -- this is the likelihood function that is being minimized
    nll = lambda *args: -log_likelihood_model(*args)

    solution = minimize(nll, initial, method='Nelder-Mead', bounds=bounds,
                        args=(wl, lum, lum_err, mask, modelcube, iontable, add_nebular),
                        options={'maxiter': 8000})

    if show_solution:
        print(solution)

    initial_optimized = solution.x.tolist()

    return initial_optimized


def plot_sesamme_samples(x, y, windowlist, flat_samples, add_nebular=True, plot_median=True, median_params=None,
                         plot_random_draws=True, n_draws=50, model_cube=None, ion_table=None, title=None, savefile_name=None):
    """
    A plotting function for examining the goodness-of-fit for models in the final sampler object after the MCMC run.

    Parameters
    ----------
    x : array-like
        Wavelength array
    y : array-like
        Flux array
    windowlist : list or array
        List of (low, high) regions to mask
    flat_samples : np.ndarray
        Flattened MCMC chain
    add_nebular : Boolean
        Determines whether to include nebular continuum emission in plotted models; optional.
    plot_median : Boolean
        Determines whether to plot a single model as a "best fit"; optional.
    median_params : list or array
        Model parameters for plot_median. The median of flat_samples is used if None; optional.
    plot_random_draws : Boolean
        Generates and plots random curves sampled from the flattened MCMC chain; optional.
    n_draws : int
        Number of random draws to plot; optional.
    model_cube : dict
        SSP cube. Falls back on the module level "modelcube" if None; optional.
    ion_table : object
        Ionization table. Falls back on the module level "iontable" if None; optional.
    title : str
        Set figure title; optional
    savefile_name : str
        Output file name if savefile set to True; optional
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    flat_samples = np.atleast_2d(np.asarray(flat_samples, dtype=float))

    # Fall back on the module level objects if the caller does not provide them
    model_cube = globals().get('modelcube') if model_cube is None else model_cube
    ion_table = globals().get('iontable') if ion_table is None else ion_table

    # Posterior median unless the caller specifies the parameters
    if median_params is None:
        median_params = np.median(flat_samples, axis=0)

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10, 6), gridspec_kw={'height_ratios': [3, 1]})

    ### Plot the data
    ax[0].step(x, y, color='black', lw=1, label="Data", zorder=1)

    ### Mark intervals that were masked during fitting (axvspan spans the axis, so it does not depend on the flux units)
    if windowlist is not None:
        for low, high in np.atleast_2d(windowlist):
            ax[0].axvspan(low, high, alpha=0.5, color='lightgrey', lw=0)
            ax[1].axvspan(low, high, alpha=0.5, color='lightgrey', lw=0)

    ### Optionally plot random draws from the final walker ensemble
    if plot_random_draws and (flat_samples.shape[0] > 1):
        rng = np.random.default_rng(99)
        idcs = rng.choice(flat_samples.shape[0], size=min(n_draws, flat_samples.shape[0]), replace=False)

        for idx in idcs:
            draw_model = get_model(flat_samples[idx], x, model_cube, ion_table, add_nebular)
            ax[0].step(x, draw_model, alpha=0.15, lw=1, ls=':', color='teal', zorder=0)
            ax[1].step(x, (y - draw_model) / y, ls=':', alpha=0.05, lw=1, color='teal', zorder=1)

    ### Optionally plot an individual model (typically a "best fit")
    total_model = None
    if plot_median:
        total_model = get_model(median_params, x, model_cube, ion_table, add_nebular)
        ax[0].step(x, total_model, color='royalblue', lw=1, label='Optimal Model', zorder=500)
        ax[1].step(x, (y - total_model) / y, color='royalblue', lw=1, zorder=2)

    ### Formatting the upper panel and setting the plot title
    ax[1].set_xlabel(r"Wavelength ($\AA$)")

    y_ref = y[np.isfinite(y) & (y != 0.)]
    ylim_high = 2 * np.median(y_ref) if y_ref.size > 0 else 1.
    ax[0].set_ylim(0, ylim_high)
    ax[0].set_ylabel(r"L$_{\odot}$ $\AA^{-1}$")

    ### Warn if the model is off the scale of the data, otherwise the curve is drawn outside the axes
    if total_model is not None:
        model_ref = total_model[np.isfinite(total_model) & (total_model != 0.)]
        if (model_ref.size > 0) and (y_ref.size > 0):
            ratio = np.median(model_ref) / np.median(y_ref)
            if not (0.01 < ratio < 100):
                print(f'The model/data median ratio is {ratio:.3e}: the model curve falls outside the plot '
                      f'limits. Check that the observation and the SSP cube share the same flux units.')

    if title is not None:
        ax[0].set_title(title, weight='semibold')

    ### Legend, adding an entry for the random draws
    handles, labels = ax[0].get_legend_handles_labels()
    if plot_random_draws and (flat_samples.shape[0] > 1):
        handles.append(Line2D([0], [0], label='Random Draw from PDF', ls=":", alpha=0.5, lw=1, color='teal'))
    ax[0].legend(loc='best', handles=handles)

    ### Lower panel formatting
    ax[1].axhline(0, color='black')
    ax[1].set_ylabel("Residuals")
    ax[1].set_ylim(-0.7, 0.7)

    plt.tight_layout()

    ### Save the file?
    if savefile_name:
        plt.savefig(savefile_name, bbox_inches='tight')
    else:
        plt.show()

    return fig, ax


def stats_line(fname, prefix='- Fit values'):
    """Read the table written by save_stats and return one line with median (+err/-err)."""
    df = pd.read_csv(fname, sep='\t', index_col='parameter')
    lo, med, hi = df['16th'].to_numpy(), df['50th'].to_numpy(), df['84th'].to_numpy()

    def fmt(i, nd=3):
        return f'{med[i]:0.{nd}f} (+{hi[i] - med[i]:0.{nd}f}/-{med[i] - lo[i]:0.{nd}f})'

    # Column 1 holds log10(Z), so exponentiate the percentiles, not the errors
    z, z_lo, z_hi = 10**med[1], 10**lo[1], 10**hi[1]

    return (f'{prefix}: log(age) = {fmt(0)}, Z = {z:0.5f} (+{z_hi - z:0.5f}/-{z - z_lo:0.5f}), '
            f'E(B-V) = {fmt(2)}, log(A) = {fmt(3)}')


bpass_pname = Path('/home/vital/PycharmProjects/specsy/examples/stellar_synthesis/SESAMME_BAPSS-AP_v3.fits')
bpass_grid = StellarBinaries.from_fits(bpass_pname)
bpass_bounds = np.array([(6, 7.5), np.log10([0.00001, 0.04]), (0.0, 3), (-20, 2)])
# bpass_grid.plot_age_metallicity()

ion_pname =  Path("/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt")
ion_table = load_ionization_table(ion_pname)

results_folder = Path('/home/vital/Astrodata/SpecSy_tests/SSP_sampling')

synth_file = './synth_sample.toml'
synth_spec = f'mock_observation_1'
synth_dict = lime.load_cfg(synth_file)[synth_spec]
print(f'- True values: log(age) = {synth_dict["log_age"]:0.3f}, Z = {synth_dict["metallicity"]:0.5f},'
                    f' E(B-V) = {synth_dict["ebv"]:0.3f}, log(A) = {np.log10(synth_dict["mass"]/1e6):0.5f}')
add_nebular = synth_dict.get('add_nebular', False)

# Load example spectrum
wave_arr, lum_arr, err_arr, mask = unpack_spectrum(results_folder/f'{synth_spec}.txt', synth_dict['distance'])

# Get the object SSP
obj_SSP_pname = results_folder/'obj_SSPs.fits'
bpass_grid.to_fits(obj_SSP_pname, wmin=wave_arr[0], wmax=wave_arr[-1])
obj_SSPs = load_ssp_cube(obj_SSP_pname)

# Extinttion law
models.set_ext_law('CCM')

# Get initial values for the fit
p0_opt = [6.5, np.log10(0.001), 0.05, 1]
p0_opt = get_initial_values(p0_opt, wave_arr, lum_arr, err_arr, mask, obj_SSPs, bpass_bounds, ion_table,
                            add_nebular=add_nebular, show_solution=False)
print(f'- p0 values: log(age) = {p0_opt[0]:0.3f}, Z = {10**p0_opt[1]:0.5f},'
                    f' E(B-V) = {p0_opt[2]:0.3f}, log(A) = {p0_opt[3]}')

# Prepare the sampler
stats.set_walker_size(128)
stats.set_chain_size(2000)
stats.set_initial_positions(p0_opt)

# Set the priors
prior_lowbounds = [6.0, np.log10(0.00001), 0.0, -20.]
prior_highbounds = [7.5, np.log10(0.04), 3.0, 2.0]
stats.set_prior_bounds(stats.prior_dict, prior_lowbounds, prior_highbounds)

# New fitting
runname = f'{synth_spec}_v1'
filename = f"SynthTesting_test_{runname}"
stats.run_sesamme(f'{filename}.h5', runname, wave_arr, lum_arr, err_arr, obj_SSPs, ion_table, mask,
                  add_nebular=add_nebular)

# Load the traces
print(f'Loading {filename} for run {runname}')
reader = emcee.backends.HDFBackend(f'{filename}.h5', name=runname)
tau = reader.get_autocorr_time(tol=0)
flat_samples = reader.get_chain(discard=500, thin=int(np.nanmean(tau)), flat=True)

# Save the measurements
vis.save_stats(flat_samples, output_path=results_folder, run_name=runname)
fname = f'{filename}_results'
plot_sesamme_samples(wave_arr, lum_arr, None, flat_samples, add_nebular=add_nebular,
                     model_cube=obj_SSPs, savefile_name=f'{fname}.pdf', ion_table=ion_table)
print(stats_line(f'{runname}_stats.txt'))