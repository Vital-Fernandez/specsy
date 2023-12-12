import lime
import specsy as sy
from pathlib import Path
import pyneb as pn

cHbeta = 0.93
S2, O2 = pn.Atom('S', 2), pn.Atom('O', 2)

log_address = Path(f'D:/Pycharm Projects/CEERs_field/reduction_v2/measurements/006649_hlsp_ceers_jwst_nirspec_nirspec7-006649_comb-mgrat_dr0.7_x1d_masked_log.txt')
log = lime.load_log(log_address)
lime.normalize_fluxes(log, norm_list='H1_4862A')

rc = pn.RedCorr(R_V=3.1, law='G03 LMC', cHbeta=cHbeta)
e_corr = rc.getCorr(log.wavelength.to_numpy(), rel_wave=4862)
log['line_int'] = log['line_flux'] * e_corr

g_int = log['line_int']

R_OII = g_int['O2_3727A']/g_int['O2_3730A']
R_SII = g_int['S2_6718A']/g_int['S2_6732A']

temp_low = 20000
neOII = O2.getTemDen(R_OII, tem=temp_low, wave1=3726, wave2=3729)
neSII = S2.getTemDen(R_SII, tem=temp_low, wave1=6717, wave2=6731)

print(f'ne_OII = {neOII} cm^-3')
print(f'ne_SII = {neSII} cm^-3')
