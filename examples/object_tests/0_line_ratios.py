import lime
# import specsy as sy
from pathlib import Path
import pyneb as pn

cHbeta = 0.93

log_address = Path(f'D:/Pycharm Projects/CEERs_field/reduction_v2/measurements/006649_hlsp_ceers_jwst_nirspec_nirspec7-006649_comb-mgrat_dr0.7_x1d_masked_log.txt')
log = lime.load_log(log_address)
lime.normalize_fluxes(log, norm_list='H1_4862A')

rc = pn.RedCorr(R_V=3.1, law='G03 LMC', cHbeta=cHbeta)
e_corr = rc.getCorr(log.wavelength.to_numpy(), rel_wave=4862)
log['line_int'] = log['line_flux'] * e_corr

gauss_fluxes = log.gauss_flux
intg_fluxes = log.intg_flux

print('Gaussian fluxes')
print(f'O3_5008A/O3_4969A', gauss_fluxes['O3_5008A']/gauss_fluxes['O3_4960A'])
print(f'N2_6584A/N2_6550A', gauss_fluxes['N2_6585A']/gauss_fluxes['N2_6550A'])
print(f'S3_9533A/S3_9071A', gauss_fluxes['S3_9533A']/gauss_fluxes['S3_9071A'])

print('Integrated fluxes')
print(f'O3_5008A/O3_4969A', intg_fluxes['O3_5008A']/intg_fluxes['O3_4960A'])
print(f'N2_6584A/N2_6550A', intg_fluxes['N2_6585A']/intg_fluxes['N2_6550A'])
print(f'S3_9533A/S3_9071A', intg_fluxes['S3_9533A']/intg_fluxes['S3_9071A'])
