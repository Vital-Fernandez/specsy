import numpy as np
from matplotlib import pyplot as plt

from sesamme import models
from sesamme.models import load_ionization_table, load_ssp_cube
from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import grid_axes, stellar_grid, nebular_shape, nebular_scale, extinction_curve
from pyneb import RedCorr


def report(label, a, b):

    """Print the worst absolute and fractional disagreement between two arrays."""

    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    diff = np.abs(a - b)
    scale = np.maximum(np.abs(a), np.abs(b))
    frac = np.where(scale > 0, diff / scale, 0.0)
    msg = f"{label:<22}\nmax abs = {diff.max():.3e} max frac = {frac.max():.3e} \n allclose = {np.allclose(a, b, rtol=1e-10, atol=0)}"

    return msg


# Input data
CUBE = '../SESAMME_BAPSS-AP_v3.fits'
Q_TABLE = '/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt'
nebx = np.array([912., 913., 1300., 1500., 1800., 2200., 2855., 3331., 3421., 3422.,
                 3642., 3648., 5700., 7000., 8207., 8209., 14583., 14585., 22787., 22789.,
                 32813., 32815., 44680., 44682.])

ion_table = load_ionization_table(Q_TABLE)
sesamme_SSPs = load_ssp_cube(CUBE)
specsy_ssp = StellarBinaries.from_fits(CUBE, source='BPASS')

# Test point
age = 7.0 #6.0
metallicity = 0.02# 0.004
ebv = 0.3
ampl = -1.5
extinction_law = 'CCM'
add_nebular = True

logZ = np.log10(metallicity)
theta = [age, logZ, ebv, ampl]

# ---------------- Spectra extraction
models.set_ext_law(extinction_law)

# Age and metallicity coordinates
met_key, met_val = models._nearest_metallicity(logZ)
age_key, age_val = models._nearest_age(age)

print(f"SESAMME selection:")
print(f' - age         {age_key}: {age_val}')
print(f' - metallicity {met_key}: {met_val}')
if not (np.isclose(age_val, age) and np.isclose(met_val, metallicity)):
    print(" - WARNING: SESAMME coordinate is far from user input")
x = np.asarray(sesamme_SSPs[met_key].data['WL'], dtype=np.float64)
y = np.asarray(sesamme_SSPs[met_key].data[age_key], dtype=np.float64)

# SpecSy interpolation
wave, flux = specsy_ssp.get_spectrum(age, metallicity)

# Plot both
fig, ax = plt.subplots()
ax.step(x, y, where='mid', label='SESAMME')
ax.step(wave, flux, where='mid', label='SpecSy', linestyle=':')
ax.set_title(report('Extinction continuum', y, flux))
ax.legend()
plt.show()

# ---------------- Extinction calculation
ones = np.ones_like(x)
atten_sesamme = models.apply_ext_law(x, ebv, ones)

klam = extinction_curve(x, extinction_law)
atten_new = 10 ** (-0.4 * ebv * klam)

# Plot both
fig, ax = plt.subplots()
ax.step(x, atten_sesamme, where='mid', label='SESAMME')
ax.step(wave, atten_new, where='mid', label='SpecSy', linestyle=':')
ax.set_title(report('Extinction continuum', atten_sesamme, atten_new))
ax.legend()
plt.show()




# ---------------- Nebular calculation
neb_sesamme = models.nebular_continuum(x, [age, logZ, ebv, ampl], ion_table)

q = nebular_scale(ion_table, [age_key], [met_key])
neb_new = 10 ** ampl * q[0] * nebular_shape(x)

# Plot both
fig, ax = plt.subplots()
ax.step(x, neb_sesamme, where='mid', label='SESAMME')
ax.step(wave, neb_new, where='mid', label='SpecSy', linestyle=':')
ax.set_title(report('Nebular continuum', neb_sesamme, neb_new))
ax.legend()
plt.show()


# Spectrum calculation
model_sesamme = models.get_model(theta, x, sesamme_SSPs, ion_table, add_nebular)

ages, mets, met_vals = grid_axes(specsy_ssp)
flat = stellar_grid(specsy_ssp, ages, met_vals)

if add_nebular:
    age_keys = [models._nearest_age(a)[0] for a in ages]
    met_keys = [models._nearest_metallicity(np.log10(m))[0] for m in met_vals]
    q_all = nebular_scale(ion_table, age_keys, met_keys)
    flat = flat + q_all[:, None] * nebular_shape(x)[None, :]

for i, age_i in enumerate(ages):
    for j, met_j in enumerate(met_vals):
        theta_ij = [age_i, np.log10(met_j), ebv, ampl]
        ref = models.get_model(theta_ij, x, sesamme_SSPs, ion_table, add_nebular)
        new = 10 ** ampl * flat[i * met_vals.size + j] * 10 ** (-0.4 * ebv * klam)
        assert np.allclose(ref, new, rtol=1e-10), f"mismatch at age {age_i}, Z {met_j}"
print(f"all {ages.size * met_vals.size} grid cells match")
