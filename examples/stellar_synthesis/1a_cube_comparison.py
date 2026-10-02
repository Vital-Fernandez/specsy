from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt, rc_context
from astropy.io import fits
import lime
import specsy as sy
from specsy.models.stellar import StellarBinaries

name: /home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       4   ()
  1  Z001          1 BinTableHDU    109   4199R x 50C
  3  Z008          1 BinTableHDU    109   4199R x 50C
  4  Z020          1 BinTableHDU    109   4199R x 50C
  5  Z040          1 BinTableHDU    109   4199R x 50C

Filename: /home/vital/Astrodata/SESAMME/pySB99_SSP_Grid_v5.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      12   ()
  1  Z000          1 BinTableHDU   1010   4200R x 500C
  2  Z00001        1 BinTableHDU   1010   4200R x 500C
  3  Z0004         1 BinTableHDU   1010   4200R x 500C
  4  Z002          1 BinTableHDU   1010   4200R x 500C
  5  Z006          1 BinTableHDU   1010   4200R x 500C
  6  Z014          1 BinTableHDU   1010   4200R x 500C
  7  Z020          1 BinTableHDU   1010   4200R x 500C



# StarBurst99 file
fname = '/home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits'
print(fits.info(fname))


# pyStarBurst99 file
fname = '/home/vital/Astrodata/SESAMME/pySB99_SSP_Grid_v5.fits'
print(fits.info(fname))


# # StarBurst99 file
# fname = '/home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits'
# starburst99 = StellarBinaries.from_fits(fname, source='starburst99')
# print(starburst99)
#
# # pyStarBurst99 file
# fname = '/home/vital/Downloads/pySB99_SESAMME_Cube_v4/pySB99_SSP_Grid_v4.fits'
# pySB99 = StellarBinaries.from_fits(fname, source='starburst99')
# print(pySB99)


wave_sb99, flux_sb99 = starburst99.get_spectrum(age=7.69, metallicity=0.02)
wave_py99, flux_py99 = pySB99.get_spectrum(age=7.69897, metallicity=0.02)

starburst99.plot_age_metallicity(fig_cfg = lime.theme.fig_defaults())
pySB99.plot_age_metallicity(fig_cfg = lime.theme.fig_defaults())

plt_cfg = {'font.size': 14, 'axes.labelsize': 16, 'xtick.labelsize': 14,
           'ytick.labelsize': 14, 'axes.titlesize': 20, 'legend.fontsize': 20}

with rc_context(plt_cfg):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.step(wave_sb99, flux_sb99, where='mid', label='Starburst99', color='C0')
    ax2.step(wave_py99, flux_py99, where='mid', label='pyStarburst99', color='C1')

    ax1.set_ylabel('Starburst99 flux', color='C0')
    ax2.set_ylabel('pyStarburst99 flux', color='C1')
    ax1.tick_params(axis='y', labelcolor='C0')
    ax2.tick_params(axis='y', labelcolor='C1')

    ax1.set_xlabel('Wavelength (Å)')
    ax1.set_title('Starburst99 versus pyStarburst99, log(age)=7.69, Z=0.020')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2)

    plt.tight_layout()
    plt.show()