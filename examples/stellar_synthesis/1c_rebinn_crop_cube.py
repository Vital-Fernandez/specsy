from pathlib import Path
import numpy as np
import lime
from specsy.models.stellar import StellarBinaries


# # Load the cube
ssp_path = Path('/home/vital/Downloads/pySB99_SESAMME_Cube_v4/pySB99_SSP_Grid_v4.fits')
pysb99_ssps = StellarBinaries.from_fits(ssp_path, source='pystarburst99')
print(pysb99_ssps)
pysb99_ssps.plot_age_metallicity()
np.savetxt(f'/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/pyStarburst99_wave_arr.txt', pysb99_ssps.wave_arr)


# Load the object
fname = '/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/SBS0335052_pyStarburst99_spectrum.txt'
spec = lime.Spectrum.from_file(fname, instrument='text')
print(spec)
spec.plot.spectrum()

# fname = '/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/SBS0335052_pyStarburst99_cube.fits'
# pysb99_ssps.to_fits(fname, disp_intvl=spec.wave.data)


fname = '/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/SBS0335052_pyStarburst99_cube.fits'
SSP_obj = StellarBinaries.from_fits(fname)
print(SSP_obj)
SSP_obj.plot_age_metallicity()