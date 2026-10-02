from pathlib import Path

import lime
from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import load_ionization_frame, ion_table_grid, nebular_continuum
from lime import Spectrum

import numpy as np
import astropy.units as u
from dust_extinction.parameter_averages import CCM89
from scipy.interpolate import RegularGridInterpolator


def interp_log_q(ion_table, log_age, Z):
    """log10 Q at (log_age, Z) by bilinear interpolation in (log Z, log age); exact on the table nodes."""
    tab_z, tab_age, q_mat = ion_table_grid(ion_table)
    return float(RegularGridInterpolator((tab_z, tab_age), q_mat)([[np.log10(Z), log_age]])[0])


def mock_observation(grid, log_age, Z, ebv, mass, dist_mpc, masks, rv=3.1, snr=5.0, wave_range=(1150, 1800),
                     m_ref=1e6, seed=42, luminosity_sun=3.828e33, MPC_CM=3.0857e24, ion_table=None,
                     plot_interp=False):
    """Mock cluster observation from a BPASS grid.

    mass        : cluster mass [Msun]
    dist_mpc    : distance [Mpc]
    ebv         : E(B-V) [mag], Cardelli+89 law
    snr         : S/N per wavelength bin (scalar or array matching the cut wavelength grid)
    masks       : (wmin, wmax) intervals in Angstrom flagged as bad (mask=True -> excluded)
    ion_table   : SESAMME ionization table (DataFrame). If given, the nebular continuum is added
    plot_interp : plot the grid spectra used in the age-metallicity interpolation
    """

    rng = np.random.default_rng(seed)

    wave, lum = grid.get_spectrum(log_age, Z, interpolate=True, plot_results=plot_interp)  # Lsun/AA for m_ref
    sel = (wave >= wave_range[0]) & (wave <= wave_range[1])
    wave, lum = wave[sel], lum[sel]

    # Nebular continuum (for m_ref, like the SSP), added before rescaling and reddening as in SESAMME
    if ion_table is not None:
        log_q = interp_log_q(ion_table, log_age, Z)
        neb = nebular_continuum(wave, [log_q])[0]
        print(f'- Nebular continuum: log(Q) = {log_q:0.3f}, nebular fraction at 1500 A = '
              f'{np.interp(1500, wave, neb / (lum + neb)):0.3f}')
        lum = lum + neb

    lum = lum * (mass / m_ref)

    # Reddening
    lum_red = lum * CCM89(Rv=rv).extinguish(wave * u.AA, Ebv=ebv)

    # Noise (Gaussian, sigma = flux / SNR)
    sigma = lum_red / snr
    lum_obs = lum_red + rng.normal(0, sigma)

    # Observed flux [erg/s/cm2/AA] -> what a telescope would record
    to_flux = luminosity_sun / (4 * np.pi * (dist_mpc * MPC_CM) ** 2)

    spec = Spectrum(wave, lum_obs * to_flux, sigma * to_flux, redshift=0)
    spec = spec.retrieve.spectrum(mask_intvls=masks)

    return spec


# Inputs
bpass_pname = Path('/home/vital/PycharmProjects/specsy/examples/stellar_synthesis/SESAMME_BAPSS-AP_v3.fits')
bpass_grid = StellarBinaries.from_fits(bpass_pname)

# Ionizing photon rates for the nebular continuum
ion_pname = Path('/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt')
ion_table = load_ionization_frame(ion_pname)

synth_file = Path('./synth_sample.toml')
cfg = lime.load_cfg(synth_file)
results_folder = Path('/home/vital/Astrodata/SpecSy_tests/SSP_sampling')

# Plot switches
plot_interp = False   # grid spectra used in the interpolation (one figure per observation)
plot_spec = False     # final mock spectrum

# Generate all the mock observations in the configuration file
for synth_spec, synth_dict in cfg.items():

    if not synth_spec.startswith('mock_observation'):
        continue

    print(f'\n{synth_spec}')
    add_nebular = synth_dict.get('add_nebular', False)

    obj_spec = mock_observation(bpass_grid,
                                synth_dict['log_age'], synth_dict['metallicity'], synth_dict['ebv'],
                                synth_dict['mass'], synth_dict['distance'], synth_dict.get('masks'),
                                ion_table=ion_table if add_nebular else None,
                                plot_interp=plot_interp)

    if plot_spec:
        obj_spec.plot.spectrum()

    obj_spec.save_spectrum(results_folder / f'{synth_spec}.txt')