import numpy as np

from specsy.models.stellar import StellarBinaries
from astropy.io import fits
from pathlib import Path
import lime

# lime.theme.set_style('dark')

# Get all files and subfolder files recursively
starburst_folder = '/home/vital/Astrodata/pystarburst99'
files = [p for p in Path(starburst_folder).rglob("*") if p.is_file()]
print(files)

# Compiled file
fname = '/home/vital/Astrodata/SESAMME/SB99_Metallicity_Table.fits'

# Starburst file
sesamme_starburst99 = StellarBinaries.from_fits(fname, source='starburst99')
print(sesamme_starburst99)
sesamme_starburst99.plot_age_metallicity(fig_cfg = lime.theme.fig_defaults())

# # Dataframe with the binary properties
# binaries_df = sesamme_starburst99.to_dataframe()
# print(binaries_df)
# np.logspace()
#
# start_age = np.log10(1e6)
# end_age = np.log10(50e6)
# n_steps = 50
# block = np.round(np.log10(np.logspace(start_age, end_age, n_steps)), 3)