from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt, rc_context
from astropy.io import fits
import lime
import specsy as sy
from specsy.models.stellar import StellarBinaries


# Generate from source
root_folder = Path('/home/vital/Astrodata/pystarburst99/compiled_spectra_v4')
pysb99_ssps = StellarBinaries.from_source('pystarburst99', source_folder=root_folder, load_spectra=True)

fname = '/home/vital/Astrodata/pystarburst99/compiled_spectra_v4/pySB99_SSP_Grid_v4.fits'
pysb99_ssps.to_fits(fname, libraries='WM')

# Reload fits
fname = '/home/vital/Astrodata/pystarburst99/compiled_spectra_v4/pySB99_SSP_Grid_v4.fits'
pysb99_ssps = StellarBinaries.from_fits(fname, source='pystarburst99')

print(pysb99_ssps)
fname = f'/home/vital/Astrodata/pystarburst99/compiled_spectra_v4/pySB99_SESAMME_age_metallicity_grid_v4.png'
pysb99_ssps.plot_age_metallicity(fig_cfg=lime.theme.fig_defaults(), fname=fname)

fname = '/home/vital/Astrodata/pystarburst99/pySTARBURST99_SESAMME_cube_v0/pySB99_SESAMMEcube_v0.fits'
pysb99_ssps_v0 = StellarBinaries.from_fits(fname, source='pystarburst99')

wave_orig, flux_orig = pysb99_ssps_v0.get_spectrum(age=6.0, metallicity=0.00200)
with rc_context():

    fig, ax = plt.subplots()
    for met in [0.00000, 0.00001, 0.00040, 0.00200, 0.00600, 0.01400, 0.02000]:
        wave_new, flux_new = pysb99_ssps.get_spectrum(age=5.991458368939762, metallicity=met)
        ax.step(wave_new, flux_new, where='mid', label=f'z = {met}, log(age) = 5.99')

    ax.legend()
    ax.set_xlabel('Wavelength $(\AA)$')
    # plt.show()
    plt.savefig(f'/home/vital/Astrodata/pystarburst99/compiled_spectra_v4/pySB99_SESAMME_spectra_comparison_v4.png')




# # fname = '/home/vital/Astrodata/pystarburst99/pySTARBURST99_SESAMME_cube_v0/pySB99_SESAMMEcube_v0.fits'
# # pysb99_ssps.to_fits(fname, libraries='WM')
#
#
# # Compiled file
#
# pyStarburst file
fname = '/home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits'
starburst99 = StellarBinaries.from_fits(fname, source='starburst99')
print(starburst99)
starburst99.plot_age_metallicity(fig_cfg = lime.theme.fig_defaults())

# wave_orig, flux_orig = starburst99.get_spectrum(age=6.0, metallicity=0.020)
# wave_new, flux_new = pysb99_ssps.get_spectrum(age=6.0, metallicity=0.020)
#
# plt_cfg = {'font.size': 14, 'axes.labelsize': 16, 'xtick.labelsize': 14,
#            'ytick.labelsize': 14, 'axes.titlesize': 20, 'legend.fontsize': 20}
#
# with rc_context(plt_cfg):
#
#     fig, ax1 = plt.subplots()
#     ax2 = ax1.twinx()
#
#     ax1.step(wave_orig, flux_orig, where='mid', label='Starburst99', color='C0')
#     ax2.step(wave_new, flux_new, where='mid', label='pyStarburst99', color='C1')
#
#     ax1.set_ylabel('Starburst99 flux', color='C0')
#     ax2.set_ylabel('pyStarburst99 flux', color='C1')
#     ax1.tick_params(axis='y', labelcolor='C0')
#     ax2.tick_params(axis='y', labelcolor='C1')
#
#     ax1.set_xlabel('Wavelength (Å)')
#     ax1.set_title('Starburst99 versus pyStarburst99, log(age)=6.0, Z=0.020')
#
#     lines1, labels1 = ax1.get_legend_handles_labels()
#     lines2, labels2 = ax2.get_legend_handles_labels()
#     ax1.legend(lines1 + lines2, labels1 + labels2)
#
#     plt.tight_layout()
#     plt.show()
#
#
#
# for metallicity in ['XMP']:
#     # data_set = root_folder / f'pySB99_{metallicity}_v0'
#     data_set = root_folder / f'{metallicity}'
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
#     idcs = (wave_SED_arr >= 900) & (wave_SED_arr <= 3000)
#
#
#     with rc_context(lime.theme.fig_defaults()):
#
#         fig, ax = plt.subplots()
#         ax2 = ax.twinx()
#
#         ax.step(wave_SED_arr[idcs], flux_matrix[0, idcs], label='Stellar + Nebular', where='mid')
#         ax2.step(wave_arr, flux_matrix_hires[0, :], label='Hi-Res', where='mid', color='C1')
#
#         # Combine legends from both axes
#         handles1, labels1 = ax.get_legend_handles_labels()
#         handles2, labels2 = ax2.get_legend_handles_labels()
#         ax.legend(handles1 + handles2, labels1 + labels2)
#         ax.set_title(f'XMP metallicity with log(age) = {np.log10(time_arr[0])}')
#
#         plt.show()