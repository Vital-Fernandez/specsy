from pathlib import Path

from matplotlib import pyplot as plt

from specsy.models.stellar import StellarBinaries
from astropy.io import fits
import numpy as np
import specsy as sy


ssp_source = 'BPASS'
bin_folder = Path('/home/vital/Astrodata/BPASS_v2.3/binaries')

nion_filenames = ['spectra-bin-imf135_300.AP.z001.a+00.dat', 'spectra-bin-imf135_300.AP.z002.a+00.dat',
                  'spectra-bin-imf135_300.AP.z003.a+00.dat', 'spectra-bin-imf135_300.AP.z004.a+00.dat',
                  'spectra-bin-imf135_300.AP.z006.a+00.dat', 'spectra-bin-imf135_300.AP.z008.a+00.dat',
                  'spectra-bin-imf135_300.AP.z010.a+00.dat', 'spectra-bin-imf135_300.AP.z014.a+00.dat',
                  'spectra-bin-imf135_300.AP.z020.a+00.dat', 'spectra-bin-imf135_300.AP.z030.a+00.dat',
                  'spectra-bin-imf135_300.AP.z040.a+00.dat', 'spectra-bin-imf135_300.AP.zem4.a+00.dat',
                  'spectra-bin-imf135_300.AP.zem5.a+00.dat']
matched_files = [f.name for f in bin_folder.iterdir() if any(term in f.name for term in nion_filenames)]

# Load SSPs binaries
bpass_SSPs = StellarBinaries.from_source(ssp_source, bin_folder, file_list=matched_files, load_spectra=True)
print(bpass_SSPs)

# Show the SSPs metallicity age range
bpass_SSPs.plot_age_metallicity()

# Save library to a certain format
output_cube = f'./BPASS-AP_Binaries_v0.fits'
bpass_SSPs.to_fits(output_cube)

output_cube = f'./SESAMME_BAPSS-AP_v3.fits'
bpass_SSPs.to_fits(output_cube, wmin=912, wmax=3000)

# Load the indexed binaries
sesamme_tutorial_ssps = StellarBinaries.from_fits(output_cube)
print(sesamme_tutorial_ssps)

# Plot binaries spectra range
sesamme_tutorial_ssps.plot_spectra(metallicity_range=[0.00001, 0.04000], age=6.5)

# Get individual spectrum
metal = 0.04
age = 6.5
wave, flux = sesamme_tutorial_ssps.get_spectrum(age=age, metallicity=metal)

fig, ax = plt.subplots()
ax.step(wave, flux, label=f'Spectrum BPASS Z={age}, log(age) = {metal}')
ax.legend()
plt.show()



# # Save formated option to external file
# BPASS_LOG_AGES = sy.cfg['stellar']['ssp']['BPASS']['z_keys']
# BPASS_LOG_AGES = {value: key for key, value in BPASS_LOG_AGES.items()}


# file_orig = "/home/vital/Astrodata/BPASS_v2.3/Demo_Met_Table.fits"
# file_new = "./SESAMME_BAPSS-AP_v1.fits"
# print(fits.info(file_orig))
# print(fits.info(file_new))
#
# with fits.open(file_orig) as hdul_orig, fits.open(file_new) as hdul_new:
#
#     # Map extname -> HDU for each file
#     orig_exts = {hdu.name: hdu for hdu in hdul_orig[1:]}
#     new_exts  = {hdu.name: hdu for hdu in hdul_new[1:]}
#
#     # Check for missing extensions
#     only_orig = set(orig_exts) - set(new_exts)
#     only_new  = set(new_exts)  - set(orig_exts)
#     if only_orig:
#         print(f'Extensions only in original: {only_orig}')
#     if only_new:
#         print(f'Extensions only in new file: {only_new}')
#
#     # Compare common extensions
#     common = set(orig_exts) & set(new_exts)
#     print(f'\nComparing {len(common)} common extensions...\n')
#
#     for extname in sorted(common):
#         orig_data = orig_exts[extname].data
#         new_data  = new_exts[extname].data
#
#         # Compare column names
#         orig_cols = set(orig_data.dtype.names)
#         new_cols  = set(new_data.dtype.names)
#         if orig_cols != new_cols:
#             print(f'[{extname}] Column mismatch: orig={orig_cols} new={new_cols}')
#             continue
#
#         # Compare data column by column
#         max_diff = {}
#         for col in orig_data.dtype.names:
#             diff = np.abs(orig_data[col].astype(float) - new_data[col].astype(float))
#             max_diff[col] = diff.max()
#
#         all_match = all(np.isclose(v, 0) for v in max_diff.values())
#         if all_match:
#             print(f'[{extname}] OK — all columns match')
#         else:
#             print(f'[{extname}] DIFFERENCES FOUND:')
#             for col, diff in max_diff.items():
#                 if not np.isclose(diff, 0):
#                     print(f'    column {col}: max absolute diff = {diff:.6e}')