from pathlib import Path
import lime
from specsy.models.stellar import StellarBinaries
from astropy.io import fits
import numpy as np
import specsy as sy


# Locate data
fits_pname = Path(f'./sdss_dr18_0358-51818-0504.fits')
spec = lime.Spectrum.from_file(fits_pname, instrument='sdss')
print(spec)

bands_pname = Path(f'./SHOC579_dr18_bands.txt')
bands = lime.load_frame(bands_pname)


# Rebin the spectrum
pixel_width = 2
spec = lime.Spectrum.from_file(fits_pname, instrument='sdss')
spec_r = spec.retrieve.rebinned(pixel_width=pixel_width, return_spectrum=True, rest_frame=True)

print(spec_r)
spec.plot.spectrum(log_scale=True, in_fig=None, rest_frame=True)
spec.plot.ax.step(spec_r.wave, spec_r.flux, label='Resampled')
spec.plot.ax.legend()
spec.plot.show()

# Mask the spectrum features
spec_cont = spec_r.retrieve.spectrum(mask_intvls=bands, obj_redshift=True)
print(spec_cont)
spec_cont.plot.spectrum(log_scale=True)

# Save the rebined masked spectrum
fname = f'./SHOC579_dr18_continuum.txt'
spec_cont.save_spectrum(fname)
spec_cont = lime.Spectrum.from_file(fname, instrument='text')
spec_cont.plot.spectrum(log_scale=True)
