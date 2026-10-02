from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from sesamme import models, vis
import sesamme.mcmc as stats
from sesamme.models import load_ionization_table, load_ssp_cube
import emcee
import corner
from astropy.table import Table
from astropy import units as u
from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import marginal_log_likelihood_profiled

# New (PyMC-marginalised) build/search -- adjust the import path to wherever sesamme_pymc.py lives
# from sesamme_pymc import build_grid, get_initial_values as get_initial_values_pymc
from specsy.models.stellar_synthesis import build_grid, get_initial_values as get_initial_values_pymc

def load_spectrum(fname, distance_mpc, col_wl='WL', col_flux='FLUX', col_err='ERROR', mask_arr=None):

    """
    Load a spectrum from an ASCII file, convert flux density to luminosity density,
    and compute a wavelength mask for spectral windows.

    Parameters
    ----------
    fname : str
        Path to the ASCII spectral file.
    distance_mpc : float
        Distance to the source in Mpc, used for flux-to-luminosity conversion.
    col_wl : str, optional
        Column name for wavelength. Default is 'WL'.
    col_flux : str, optional
        Column name for flux. Default is 'FLUX'.
    col_err : str, optional
        Column name for flux error. Default is 'ERROR'.
    windowlist : array-like or None, optional
        2D array of [wl_min, wl_max] pairs defining the clean spectral windows.
        If None, defaults to the M83-8 window list.

    Returns
    -------
    wl : array
        Wavelength array.
    lum : array
        Luminosity density array (L_Sun / A).
    lum_err : array
        Luminosity density uncertainty array (L_Sun / A).
    mask : array
        Boolean mask array for the spectral windows.
    """

    specfile = Table.read(fname, format='ascii')
    wl = specfile[col_wl]
    flux = specfile[col_flux]
    flux_err = specfile[col_err]

    distance_cm = (distance_mpc * u.Mpc).to(u.cm).value
    conversion = 4 * np.pi * distance_cm**2 / 3.83e33
    lum = flux * conversion
    lum_err = flux_err * conversion
    mask = models.get_mask(mask_arr, wl)

    return wl, lum, lum_err, mask


def log_likelihood_model(theta, wl, y, yerr, mask, model_cube, ion_table, add_nebular):
    y_model = models.get_model(theta, wl, model_cube, ion_table, add_nebular)

    masked_spec = np.array((y[mask] - y_model[mask]) / yerr[mask])
    masked_err = np.array(np.sqrt(2) * np.sqrt(np.pi) * yerr[mask])

    resid = -0.5 * (np.dot(masked_spec, masked_spec) + np.log(np.dot(masked_err, masked_err)))

    return resid


def get_initial_values(initial, wl, lum, lum_err, mask, modelcube, bounds, show_solution=True):

    # -- this is the likelihood function that is being minimized
    nll = lambda *args: -log_likelihood_model(*args)

    solution = minimize(nll, initial, method='Nelder-Mead', bounds=bounds,
                        args=(wl, lum, lum_err, mask, modelcube, None, False),
                        options={'maxiter': 8000})

    if show_solution:
        print(solution)

    initial_optimized = solution.x.tolist()

    return initial_optimized


# Load example spectrum
fname = '/home/vital/PycharmProjects/SESAMME/docs/example_data/M83-8_preprocess.txt'
windowlist = np.array([[912, 1133], [1172, 1176.5],  [1188, 1202], [1203.8, 1222],
                       [1257, 1262], [1299, 1305], [1331, 1336.2], [1276, 1286], [1454, 1456], [1465, 1469],
                       [1390.5, 1393.5], [1399.5, 1403.3],
                       [1523.5, 1529], [1543, 1550.], [1608.5, 1621], [1656, 1659], [1666, 1674], [1795, np.max(3000)] ])
wl, lum, lum_err, mask = load_spectrum(fname, distance_mpc=4.8, mask_arr=windowlist)


# # Load the stellar binaries and prepare them for the object
bpass_ssp_pname = Path(f'../SESAMME_BAPSS-AP_v3.fits')
bpass_SSPs = load_ssp_cube(bpass_ssp_pname)
ion_table = load_ionization_table("/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt")
specsy_ssp = StellarBinaries.from_fits(bpass_ssp_pname, source='BPASS')

# Bounds
bpass_bounds = np.array([(6, 7.5), np.log10([1e-05, 0.04]), (0.0, 0.2),  (-3, 2)])

p0_init = [6.5, np.log10(0.001), 0.05, -0.7]

t0 = time.perf_counter()
p0_opt = get_initial_values(p0_init, wl, lum, lum_err, mask, bpass_SSPs, bpass_bounds, show_solution=True)
t_emcee = time.perf_counter() - t0
print(p0_opt)


# =============================================================================
# New: marginalised (PyMC-style) initial value search, same cube and mask
# =============================================================================

age_bounds = tuple(bpass_bounds[0])
met_bounds = tuple(bpass_bounds[1])
ebv_bounds = tuple(bpass_bounds[2])
centre = np.log10(np.median(np.asarray(lum)[mask]))
ampl_bounds = (centre - 3, centre + 3)

# add_nebular=False here to match the original get_initial_values() call above,
# which hardcodes add_nebular=False (and passes ion_table=None, unused as a result)
t0 = time.perf_counter()
flat, klam, ages, mets = build_grid(wl, specsy_ssp, ion_table, add_nebular=False,
                                    age_bounds=age_bounds, met_bounds=met_bounds)

p0_opt_pymc = get_initial_values_pymc(lum, lum_err, mask, flat, klam, ages, mets,
                                      initial=(p0_init[2], p0_init[3]),
                                      bounds=(ebv_bounds, ampl_bounds),
                                      show_solution=True)
t_pymc = time.perf_counter() - t0
print(p0_opt_pymc)


# =============================================================================
# Comparison
# =============================================================================

labels = ['log(age)', 'log(Z)', 'E(B-V)', 'log(A)']
print(f"\n{'':<10}{'emcee-style':>14}{'marginalised':>14}{'diff':>10}")
for label, a, b in zip(labels, p0_opt, p0_opt_pymc):
    print(f"{label:<10}{a:>14.4f}{b:>14.4f}{a - b:>10.4f}")

print(f"\nemcee-style search : {t_emcee:6.2f} s "
     f"({bpass_bounds.shape[0]}D Nelder-Mead, one get_model() call per iteration)")
print(f"marginalised search: {t_pymc:6.2f} s "
     f"(grid build once + 2D Nelder-Mead over the whole grid per iteration)")


print(specsy_ssp.wave_arr.shape, np.asarray(wl).shape)
assert np.allclose(specsy_ssp.wave_arr, wl), "model and observed wavelength grids differ"



print(np.allclose(specsy_ssp.wave_arr, np.asarray(wl)))
print(specsy_ssp.wave_arr[:3], np.asarray(wl)[:3])
print(specsy_ssp.wave_arr[-3:], np.asarray(wl)[-3:])

print(flat.shape, ages.size, mets.size, ages.size * mets.size)

i, j = 15, 0                      # the cell the argmax picked
w, f = specsy_ssp.get_spectrum(ages[i], np.power(10, mets[j]))
print("row matches get_spectrum:", np.allclose(flat[i*mets.size + j], f))
print(f"model median = {np.median(f):.3e}, data median = {np.median(lum):.3e}")

print(f"masked wl range: {np.asarray(wl)[mask].min()} - {np.asarray(wl)[mask].max()}")
print(f"zeros inside mask: {np.count_nonzero(np.asarray(lum)[mask] == 0)} / {np.sum(mask)}")
print(f"nonzero data range: {np.asarray(wl)[np.asarray(lum) != 0].min()} - "
      f"{np.asarray(wl)[np.asarray(lum) != 0].max()}")

norm = np.median(flat[:, mask], axis=1)
print(f"model medians span: {norm.min():.3e} to {norm.max():.3e} "
      f"({np.log10(norm.max()/norm.min()):.2f} dex)")

ebv0 = p0_opt_pymc[2]
_, ll, A_k = marginal_log_likelihood_profiled(ebv0, np.asarray(lum), np.asarray(lum_err),
                                              mask, flat, klam, np.log(np.full(flat.shape[0], 1/flat.shape[0])))
p_k = np.exp(ll - logsumexp(ll))
joint = p_k.reshape(ages.size, mets.size)

print(f"best -2logL = {-2*ll.max():.4e}, chi2/dof = {-2*ll.max()/np.sum(mask):.2f}")
print(f"ll spread = {ll.max() - ll.min():.4e}")
print("met marginal:", np.array2string(joint.sum(axis=0), precision=3))
print("age marginal:", np.array2string(joint.sum(axis=1), precision=3))

best = int(np.argmax(ll))
mod_best = A_k[best] * flat[best, mask] * 10**(-0.4 * ebv0 * klam[mask])
resid = (np.asarray(lum)[mask] - mod_best) / np.asarray(lum_err)[mask]
print(f"residual mean = {resid.mean():.2f}, std = {resid.std():.2f}")

_, ll, A_k = marginal_log_likelihood_profiled(ebv0, np.asarray(lum), np.asarray(lum_err),
                                              mask, flat, klam,
                                              np.log(np.full(flat.shape[0], 1/flat.shape[0])))
k_map = int(np.argmax(ll))
i, j = np.unravel_index(k_map, (ages.size, mets.size))
print(f"k_map = {k_map}, i = {i}, j = {j}, age = {ages[i]}, met = {mets[j]}")
print(f"A_k[k_map] = {A_k[k_map]:.4e}, row_norm[k_map] = {np.median(flat[k_map, mask]):.4e}")
print(f"ampl = {np.log10(A_k[k_map] * np.median(flat[k_map, mask])):.4f}")

print(f"centre = {centre:.4f}, ampl returned = {p0_opt_pymc[3]:.4f}")

import matplotlib.pyplot as plt
wl_m = np.asarray(wl)[mask]
fig, ax = plt.subplots(figsize=(11, 4))
ax.plot(wl_m, resid, lw=0.7)
ax.axhline(0, color='k', lw=0.8)
ax.set_xlabel('Wavelength (Å)'); ax.set_ylabel('residual / σ')
plt.show()