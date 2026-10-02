from pathlib import Path
from specsy.models.ssp import StellarBinaries
import numpy as np


ssp_source = 'VMS_BPASS'
bin_folder = Path('/home/vital/Astrodata/VMS_BPASS/bursts_MP25')

names = [f.name for f in bin_folder.iterdir() if f.is_file()]

# Load SSPs binaries
VMS_SSPs = StellarBinaries.from_source(ssp_source, bin_folder, file_list=names, load_spectra=True)
VMS_SSPs.plot_age_metallicity()
print(VMS_SSPs)