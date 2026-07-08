from pathlib import Path

from matplotlib import pyplot as plt, rc_context

import lime
from specsy.models.stellar import StellarBinaries
from astropy.io import fits
import numpy as np
import specsy as sy


root_folder = Path('/home/vital/Astrodata/pystarburst99/compiled_spectra_v0')


# pysb99_ssps = StellarBinaries.from_source('pystarburst99', source_folder=root_folder, load_spectra=True)

fname = '/home/vital/Astrodata/pystarburst99/pySTARBURST99_SESAMME_cube_v0/pySB99_SESAMMEcube_v0.fits'
pysb99_ssps = StellarBinaries.from_fits(fname, source='pystarburst99')

print(pysb99_ssps)
pysb99_ssps.plot_age_metallicity(fig_cfg=lime.theme.fig_defaults())

# fname = '/home/vital/Astrodata/pystarburst99/pySTARBURST99_SESAMME_cube_v0/pySB99_SESAMMEcube_v0.fits'
# pysb99_ssps.to_fits(fname, libraries='WM')


# Compiled file

# pyStarburst file
fname = '/home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits'
starburst99 = StellarBinaries.from_fits(fname, source='starburst99')
print(starburst99)
starburst99.plot_age_metallicity(fig_cfg = lime.theme.fig_defaults())

wave_orig, flux_orig = starburst99.get_spectrum(age=6.0, metallicity=0.020)
wave_new, flux_new = pysb99_ssps.get_spectrum(age=6.0, metallicity=0.020)

plt_cfg = {'font.size': 14, 'axes.labelsize': 16, 'xtick.labelsize': 14,
           'ytick.labelsize': 14, 'axes.titlesize': 20, 'legend.fontsize': 20}

with rc_context(plt_cfg):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.step(wave_orig, flux_orig, where='mid', label='Starburst99', color='C0')
    ax2.step(wave_new, flux_new, where='mid', label='pyStarburst99', color='C1')

    ax1.set_ylabel('Starburst99 flux', color='C0')
    ax2.set_ylabel('pyStarburst99 flux', color='C1')
    ax1.tick_params(axis='y', labelcolor='C0')
    ax2.tick_params(axis='y', labelcolor='C1')

    ax1.set_xlabel('Wavelength (Å)')
    ax1.set_title('Starburst99 versus pyStarburst99, log(age)=6.0, Z=0.020')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2)

    plt.tight_layout()
    plt.show()



# for metallicity in ['XMP']:
#     data_set = root_folder / f'pySB99_{metallicity}_v0'
#
#     fname = data_set/'SED_wavelength.txt'
#     wave_SED_arr = np.loadtxt(fname)
#
#     fname = data_set/'timesteps.txt'
#     time_arr = np.loadtxt(fname)
#
#     fname = data_set/'pySB_SED_stellar_and_nebular.npy'
#     flux_matrix = np.load(fname)
#
#     fname = data_set/'spectrum_wavelength.txt'
#     wave_arr = np.loadtxt(fname)
#
#     fname = data_set/'pySB_hires_spectrum.npy'
#     flux_matrix_hires = np.load(fname)
#
#     fig, ax = plt.subplots()
#     ax.step(wave_arr, flux_matrix_hires[0, :], label='Hi-Res', where='mid')
#     plt.show()