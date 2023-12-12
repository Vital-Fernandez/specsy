import lime
import specsy as sy
from pathlib import Path
import pyneb as pn
from specsy.astro.chemistry import sufur_diaz_2022

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

S_H = sufur_diaz_2022(log, S2_lines=("S2_6718A", "S2_6732A"), S3_lines=("S3_9071A", "S3_9533A"),
                S2_norm="H1_6565A", S3_norm="H1_6565A", flux_column='line_int')

print(f'12 + log(S/H) = {S_H[0]:0.2f} + {S_H[1]:0.2f}')