import lime
import specsy as sy
from pathlib import Path
import pyneb as pn

# Parameters
cHbeta = 0.93
S3, O3 = pn.Atom('S', 3), pn.Atom('O', 3)

# Load the fluxes
log_address = Path(f'D:/Pycharm Projects/CEERs_field/reduction_v2/measurements/006649_hlsp_ceers_jwst_nirspec_nirspec7-006649_comb-mgrat_dr0.7_x1d_masked_log.txt')
log = lime.load_log(log_address)
lime.normalize_fluxes(log, norm_list='H1_4862A')

# Extinction calculation
rc = pn.RedCorr(R_V=3.1, law='G03 LMC', cHbeta=cHbeta)
e_corr = rc.getCorr(log.wavelength.to_numpy(), rel_wave=4862)
log['line_int'] = log['line_flux'] * e_corr
log['line_int_err'] = log['line_flux_err'] * e_corr

RS3 = log.loc['S3_9533A', 'line_int']/log.loc['S3_9071A', 'line_int']
RSIII_a = log.loc['S3_6314A', 'line_int']/log.loc['S3_9071A', 'line_int']
RSIII_b = log.loc['S3_6314A', 'line_int']/log.loc['S3_9533A', 'line_int']

log.loc['S3_9533A', 'intg_flux']/log.loc['S3_9071A', 'intg_flux']

den_low = 100
TSIII_a = S3.getTemDen(RSIII_a, den=den_low, wave1=6312, wave2=9069)
TSIII_b = S3.getTemDen(RSIII_b, den=den_low, wave1=6312, wave2=9531)

print(f'TSIII_a = {TSIII_a}')
print(f'TSIII_a = {TSIII_b}')