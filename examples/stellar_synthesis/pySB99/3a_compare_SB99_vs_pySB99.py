from pathlib import Path
import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt, rc_context
from specsy.models.stellar import StellarBinaries
from astropy.io import fits
import lime

# Starburst99
logan_pname = Path('/home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits')
logan_cube = StellarBinaries.from_fits(logan_pname, 'starburst99')

#pyStarburst99
vital_pname = Path('/home/vital/Astrodata/SESAMME/pySB99_SSP_Grid_v5.fits')
vital_cube = StellarBinaries.from_fits(vital_pname, 'starburst99')
vital_cube.plot_age_metallicity()

# Original files
ages_pname = Path('./pysb_testmod/timesteps.txt')
wave_pname = Path('./pysb_testmod/spectrum_wavelength.txt')
flux_pname = Path('./pysb_testmod/pySB_hires_spectrum.npy')
wave_arr = np.loadtxt(wave_pname)
ages_arr = np.loadtxt(ages_pname)
flux_grid = np.load(flux_pname)

target_ages = np.array([1000000, 10000000, 48900000])
metallicity = 0.02


fig_cfg = lime.theme.fig_defaults(user_fig={"figure.figsize" : (8, 8), "figure.dpi" : 100, "axes.titlesize": 20,
                                            "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16,
                                            "legend.fontsize" : 16})
with rc_context(fig_cfg):

    fig, ax1 = plt.subplots()
    colors = cm.viridis(np.linspace(0, 1, len(target_ages)))
    for color, age in zip(colors, target_ages):

        idx = np.searchsorted(np.power(10, logan_cube.age_log), age)
        age_logan = logan_cube.age_log[idx]
        wave_logan, flux_logan = logan_cube.get_spectrum(age=age_logan, metallicity=metallicity)

        idx = np.searchsorted(np.power(10, vital_cube.age_log), age)
        age_vital = vital_cube.age_log[idx]
        wave_vital, flux_vital = vital_cube.get_spectrum(age=age_vital, metallicity=metallicity)

        ax1.step(wave_logan, flux_logan, where='mid', color=color, label=f'log(age) = {age_logan:0.4f} (Starburst)')
        ax1.step(wave_vital, flux_vital, where='mid', color=color, linestyle='--', label=f'log(age) = {age_vital:0.4f} (pySB99)')

    ax1.set_title('Starburst99 versus pyStarburst99 (Z = 0.02)')
    ax1.set_ylabel('Flux $L_{sun}$')
    ax1.set_xlabel('Wavelength (Å)')
    ax1.set_yscale('log')

    ax1.legend(loc='upper right')
    plt.show()

