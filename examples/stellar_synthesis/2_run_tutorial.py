from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

from specsy.models.stellar import StellarBinaries
from sesamme import models, vis
import sesamme.mcmc as stats
from sesamme.models import load_ionization_table, load_ssp_cube
import emcee
import corner
from astropy.table import Table
from astropy import units as u


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


# Load example spectrum
fname = '/home/vital/PycharmProjects/SESAMME/docs/example_data/M83-8_preprocess.txt'
windowlist = np.array([[912, 1133], [1172, 1176.5],  [1188, 1202], [1203.8, 1222],
                       [1257, 1262], [1299, 1305], [1331, 1336.2], [1276, 1286], [1454, 1456], [1465, 1469],
                       [1390.5, 1393.5], [1399.5, 1403.3],
                       [1523.5, 1529], [1543, 1550.], [1608.5, 1621], [1656, 1659], [1666, 1674], [1795, np.max(3000)] ])
wl, lum, lum_err, mask = load_spectrum(fname, distance_mpc=4.8, mask_arr=windowlist)


# # Load the stellar binaries and prepare them for the object
bpass_SSPs = load_ssp_cube(f'./SESAMME_BAPSS-AP_v3.fits')
ion_table = load_ionization_table("/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt")

# Extinttion law
models.set_ext_law('CCM')

# Prepare the sampler
stats.set_walker_size(128)
stats.set_chain_size(5000)
stats.set_initial_positions([7., -2.1, 0.3, -2.])

# Set the priors
prior_lowbounds = [6.0, np.log10(0.008), 0.01, -20.]
prior_highbounds = [7.5, np.log10(0.03), 1.0, 1.0]
stats.set_prior_bounds(stats.prior_dict, prior_lowbounds, prior_highbounds)

# New fitting
runname = 'specsy_script_v4'
filename = f"M83_test_{runname}.h5"
stats.run_sesamme(filename, runname, wl, lum, lum_err, bpass_SSPs, ion_table, mask, True)

# Load the traces
print(f'Loading {filename} for run {runname}')
reader = emcee.backends.HDFBackend(filename, name=runname)
flat_samples = reader.get_chain(discard=40, thin=10, flat=True)

# Save the measurements
fname = f'{filename}_results'
vis.save_stats(flat_samples, output_path='./', run_name=runname)