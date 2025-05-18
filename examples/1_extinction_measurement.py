import pandas as pd

import lime
import specsy as sy
import numpy as np
import pyneb as pn

# Load the line measurements
log_file = './sample_data/manga_lines_log.txt'
log_df = lime.load_frame(log_file)
R_V = 3.1
law = 'G03 LMC'

print(log_df.index.to_numpy())

sy.theme.set_style('dark', library='bokeh')


cHbeta, cHbeta_err, log_extinction = sy.extinction_coeff_calc(log_df, norm_line='H1_4861A', R_V=R_V, law=law,
                                                              exclude_list=['H1_3889A_m', 'H1_3970A'],
                                                              plot_results=True)

print(cHbeta, cHbeta_err)

# idcs_line =
# log_theo = pd.DataFrame
#
# tem, den = 10000, 100
# H1 = pn.RecAtom('H', 1)
#
# # Compare with pyneb
# rc = pn.RedCorr(R_V=R_V, law=law)
# obs_ratio = log_df.loc['H1_6563A', 'profile_flux']/log_df.loc['H1_4861A', 'profile_flux']
# theo_ratio = H1.getEmissivity(tem, den, wave=6563)/H1.getEmissivity(tem, den, wave=4861)
# rc.setCorr(obs_over_theo=obs_ratio/theo_ratio, wave1=6563., wave2=4861.)
# print(rc.cHbeta)
# int_Halpha = obs_ratio * rc.getCorrHb(6563)
# print(int_Halpha)
#


# sy.extinction.reddening_correction(cHbeta, cHbeta_err, log_df, norm_wavelength=4861.0, flux_column='line_flux')
# lime.save_frame('./sample_data/gp121903_log_norm_int.txt', log_df)
#
# # Compare with Hbeta
# rc = pn.RedCorr(cHbeta=cHbeta, R_V=3.1, law='G03 LMC')
# sy_corr = log_df.loc["H1_6563A", "line_int"]/log_df.loc["H1_6563A", "line_flux"]
#
# print(f'\nPyneb correction: {rc.getCorr(6563)}')
# print(f'LiMe correction: {1/sy_corr}')


# rc.setCorr(obs_over_theo=5.34/2.86, wave1=6563., wave2=4861.)
# ratio_dis = np.random.normal(obsRatio_uarray[-1].nominal_value, obsRatio_uarray[-1].std_dev, size=1000)/2.86
# rc.setCorr(obs_over_theo=ratio_dis, wave1=6563., wave2=4861.)
# cHb, cHb_err = rc.cHbeta.mean(), rc.cHbeta.std() #(0.8395045358309076, 0.15112567990954212)
# eBV, eBV_err = rc.EbvFromCHbeta(cHb), rc.EbvFromCHbeta(cHb_err)
# print(f'E(B-V) = {eBV:0.2f} +/- {eBV_err:0.2f} || c(Hbeta) = {cHb:0.2f} +/- {cHb_err:0.2f}')
