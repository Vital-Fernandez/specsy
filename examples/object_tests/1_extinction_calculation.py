import lime
import specsy as sy
from pathlib import Path


# Load the data
log_address = 'hlsp_ceers_jwst_nirspec_nirspec7-006649_comb-mgrat_v0.7_x1d-masked_log.txt'
log = sy.load_log(log_address)
lime.normalize_fluxes(log, norm_list='H1_4862A')

# Calculate the extinction coefficients
cHbeta, cHbeta_err = sy.cHbeta_from_log(log, ref_line='H1_4862A', flux_entry='gauss', show_plot=True)

# Get line intensities
sy.reddening_correction(cHbeta, cHbeta_err, log, norm_wavelength=4862, flux_column='line_flux')


print('Intensity ratios')
int_series = log['line_int']
print(f'O3_5008A/O3_4969A', int_series['O3_5008A']/int_series['O3_4960A'])
print(f'N2_6584A/N2_6550A', int_series['N2_6585A']/int_series['N2_6550A'])
print(f'S3_9533A/S3_9071A', int_series['S3_9533A']/int_series['S3_9071A'])

