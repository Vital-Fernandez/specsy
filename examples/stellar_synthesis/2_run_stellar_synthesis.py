from pathlib import Path
import lime
from specsy.models.stellar import StellarBinaries


# Load the stellar binaries
ssps_fname = f'./BPASS-AP_Binaries_v0.fits'
bpass_SSPs = StellarBinaries.from_fits(ssps_fname)
print(bpass_SSPs)

# Save the rebined masked spectrum
fname = f'./SHOC579_dr18_continuum.txt'
spec_cont = lime.Spectrum.from_file(fname, instrument='text')
print(spec_cont)
spec_cont.plot.spectrum(log_scale=True)

# Prepare the binaries for the object
fname = f'./SHOC579_dr18_bpass_ssp.fits'
bpass_SSPs.to_fits(fname, disp_intvl=spec_cont.wave_rest.data) #, wmin=wmin, wmax=wmax, pixel_width=pixel_width)
sdss_bpass = StellarBinaries.from_fits(fname)
print(sdss_bpass)