from pathlib import Path
import numpy as np
import pymc as pm
import arviz as az
from matplotlib import pyplot as plt

from astropy.table import Table
from astropy import units as u

from sesamme import models
from sesamme.models import load_ionization_table

from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import (build_grid, build_model, set_prior_bounds,
                                             get_initial_values, posterior_grid, credible_cells)


def load_spectrum(fname, distance_mpc, col_wl='WL', col_flux='FLUX', col_err='ERROR', mask_arr=None):

    """
    Load a spectrum from an ASCII file, convert flux density to luminosity density,
    and compute a wavelength mask for spectral windows.
    """

    specfile = Table.read(fname, format='ascii')
    wl = specfile[col_wl]
    flux = specfile[col_flux]
    flux_err = specfile[col_err]

    distance_cm = (distance_mpc * u.Mpc).to(u.cm).value
    conversion = 4 * np.pi * distance_cm ** 2 / 3.83e33
    lum = flux * conversion
    lum_err = flux_err * conversion
    mask = models.get_mask(mask_arr, wl)

    return np.asarray(wl), np.asarray(lum), np.asarray(lum_err), mask


# =============================================================================
# Inputs
# =============================================================================

fname = '/home/vital/PycharmProjects/SESAMME/docs/example_data/M83-8_preprocess.txt'
windowlist = np.array([[912, 1133], [1172, 1176.5],  [1188, 1202], [1203.8, 1222],
                       [1257, 1262], [1299, 1305], [1331, 1336.2], [1276, 1286], [1454, 1456], [1465, 1469],
                       [1390.5, 1393.5], [1399.5, 1403.3],
                       [1523.5, 1529], [1543, 1550.], [1608.5, 1621], [1656, 1659], [1666, 1674], [1795, np.max(3000)]])
windowlist = np.array([w for w in windowlist if not (1390 < w[0] < 1405 or 1520 < w[0] < 1555)])

wl, lum, lum_err, mask = load_spectrum(fname, distance_mpc=4.8, mask_arr=windowlist)

bpass_ssp_pname = Path('../SESAMME_BAPSS-AP_v3.fits')
ion_table = load_ionization_table("/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt")
specsy_ssp = StellarBinaries.from_fits(bpass_ssp_pname, source='BPASS')


# =============================================================================
# Priors and grid
# =============================================================================

prior_dict = set_prior_bounds(prior_dict=dict.fromkeys(['age', 'met', 'ebv', 'amp']),
                              prior_lowbounds=[6.0, np.log10(1e-5), 0.0, 0.0],
                              prior_highbounds=[7.5, np.log10(0.04), 0.2, 1.0])

# ampl multiplies unit-median rows, so its bounds come from the data scale,
# not from SESAMME's mass-scale log(A) convention
centre = np.log10(np.median(lum[mask]))
prior_dict['amp'] = [centre - 3, centre + 3]

flat, klam, ages, mets = build_grid(wl, specsy_ssp, ion_table, ext_law='CCM', add_nebular=False, prior_bounds=prior_dict)

print(f"grid: {flat.shape[0]} spectra ({ages.size} ages x {mets.size} metallicities), "
      f"{np.sum(mask)} fitted pixels")
print(f"ampl prior: [{prior_dict['amp'][0]:.2f}, {prior_dict['amp'][1]:.2f}]")


# =============================================================================
# Starting point and sampling
# =============================================================================

age0, met0, ebv0, ampl0 = get_initial_values(lum, lum_err, mask, flat, klam, ages, mets,
                                             initial=(0.1, centre),
                                             bounds=(tuple(prior_dict['ebv']), tuple(prior_dict['amp'])),
                                             show_solution=False)
print(f"start: ebv = {ebv0:.4f}, ampl = {ampl0:.4f}  (best cell: age {age0}, log Z {met0})")

model = build_model(wl, lum, lum_err, flat, klam, mask, prior_dict)


with model:
    idata = pm.sample(draws=200, tune=1000, chains=10, cores=10, target_accept=0.9,
                      initvals={"ebv": ebv0, "ampl": ampl0, "lnf": 1.3},
                      nuts_sampler='numpyro', progressbar='combined',
                      random_seed=42)

print(az.summary(idata, var_names=["ebv", "ampl", "lnf"]))
print(f"error inflation factor = {np.exp(float(idata.posterior['lnf'].mean())):.2f}")


# =============================================================================
# Diagnostics
# =============================================================================



print(az.summary(idata, var_names=["ebv", "ampl"]))

n_div = int(idata.sample_stats["diverging"].values.sum())
if n_div > 0:
    print(f"WARNING: {n_div} divergences -- raise target_accept or check the ebv/ampl correlation")

print(np.corrcoef(idata.posterior["ebv"].values.ravel(), idata.posterior["ampl"].values.ravel())[0, 1])

az.plot_trace(idata, var_names=["ebv", "ampl"])
plt.tight_layout()

az.plot_pair(idata, var_names=["ebv", "ampl"])
plt.tight_layout()




# =============================================================================
# Posterior over the grid
# =============================================================================

res = posterior_grid(idata, {"ages": ages, "mets": mets, "flat": flat, "klam": klam})

print(f"\nlog(age) MAP = {res['age_map']:.2f}, weighted mean = {res['age_mean']:.2f}")
print(f"log(Z)   MAP = {res['met_map']:.2f}, weighted mean = {res['met_mean']:.2f}")
print(f"log(A)       = {res['logA']:.2f}")

top = np.argsort(res['joint'].ravel())[::-1][:6]
for k in top:
    i, j = np.unravel_index(k, res['joint'].shape)
    print(f"log(age) = {ages[i]:.1f}, log(Z) = {mets[j]:.2f}, p = {res['joint'][i, j]:.4f}")
print(f"90% set: {len(credible_cells(res['joint'], 0.90))} cells")

print(f"\nlog(age) MAP = {res['age_map']:.2f}, weighted mean = {res['age_mean']:.2f}")
print(f"log(Z)   MAP = {res['met_map']:.2f}, weighted mean = {res['met_mean']:.2f}")
print(f"log(A)       = {res['logA']:.2f}")

cells = credible_cells(res['joint'], level=0.68)
print(f"68% credible set: {len(cells)} of {ages.size * mets.size} cells")
for i, j in cells:
    print(f"   log(age) = {ages[i]:.1f}, log(Z) = {mets[j]:.2f}, p = {res['joint'][i, j]:.3f}")

# Joint posterior over the grid
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(res['joint'].T, origin='lower', aspect='auto', cmap='viridis',
               extent=[ages.min(), ages.max(), mets.min(), mets.max()])
ax.set_xlabel('log(age / yr)')
ax.set_ylabel('log(Z)')
ax.set_title('Joint posterior over the model grid')
fig.colorbar(im, ax=ax, label='P(age, Z | data)')
plt.tight_layout()

# Marginals -- bars, not curves: the distribution is discrete
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
ax1.bar(ages, res['age_pdf'], width=0.06)
ax1.set_xlabel('log(age / yr)'); ax1.set_ylabel('P')
ax2.bar(mets, res['met_pdf'], width=0.12)
ax2.set_xlabel('log(Z)')
plt.tight_layout()

# =============================================================================
# Best-fit spectrum and residuals
# =============================================================================

ampl_med = float(idata.posterior["ampl"].median())
ebv_med = float(idata.posterior["ebv"].median())
row_norm = idata.constant_data["row_norm"].values
klam_mean = float(idata.constant_data["klam_mean"].values)

k_map = int(np.argmax(res['joint']))
mod = 10 ** ampl_med * (flat[k_map] / row_norm[k_map]) * 10 ** (-0.4 * ebv_med * (klam - klam_mean))

for name, lo, hi in [('Si IV', 1390, 1406), ('C IV', 1520, 1555)]:
    sel = mask & (wl > lo) & (wl < hi)
    if sel.sum():
        r = (lum[sel] - mod[sel]) / lum_err[sel]
        print(f"{name}: {sel.sum()} px, resid mean {r.mean():+.2f}, std {r.std():.2f}")

print(np.array2string(specsy_ssp.metallicities, precision=5))
print(f"mets used: {np.array2string(mets, precision=2)}")


# Residuals over the fitted pixels only; NaN elsewhere so the masked stretches
# appear as gaps rather than as lines drawn across them
resid_full = np.full(wl.size, np.nan)
resid_full[mask] = (lum[mask] - mod[mask]) / lum_err[mask]

fig, (ax, axr) = plt.subplots(2, 1, sharex=True, figsize=(13, 7),
                              gridspec_kw={'height_ratios': [3, 1]})

# Masked regions, shaded on both panels
for w0, w1 in windowlist:
    ax.axvspan(w0, w1, color='0.88', zorder=0)
    axr.axvspan(w0, w1, color='0.88', zorder=0)

ax.step(wl, lum, where='mid', lw=0.8, color='k', label='data', zorder=2)
ax.step(wl, mod, where='mid', lw=1.0, color='C0', label='median model', zorder=3)
ax.set_xlim(wl[mask].min() - 20, wl[mask].max() + 20)
ax.set_ylim(0, np.nanpercentile(lum[mask], 99.5) * 1.3)
ax.set_ylabel(r'$L_\odot\ \AA^{-1}$')
ax.set_title(f"log(age) = {res['age_map']:.1f}, log(Z) = {res['met_map']:.2f}, "
             f"E(B-V) = {ebv_med:.3f}, log(A) = {res['logA']:.2f}")
ax.legend(loc='upper right')

axr.step(wl, resid_full, where='mid', lw=0.7, color='C3', zorder=2)
axr.axhline(0, color='k', lw=0.8, zorder=1)
axr.set_xlabel(r'Wavelength ($\AA$)')
axr.set_ylabel(r'residual / $\sigma$')

# Label the features that carry the age and metallicity signal
lines = {'N V': 1240., 'Si IV': 1397., 'C IV': 1549., 'He II': 1640.}
for name, w in lines.items():
    if ax.get_xlim()[0] < w < ax.get_xlim()[1]:
        ax.axvline(w, color='C1', lw=0.8, ls=':', zorder=1)
        ax.text(w, ax.get_ylim()[1] * 0.95, name, rotation=90, fontsize=8,
                ha='right', va='top', color='C1')

resid = resid_full[mask]
print(f"\nresidual mean = {resid.mean():.2f}, std = {resid.std():.2f}")

plt.tight_layout()
plt.show()

