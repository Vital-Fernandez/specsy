import pandas as pd

import lime
import specsy
import specsy as sy
import numpy as np
import pyneb as pn

def pyneb_cHbeta_calc(measurements_df, labelA, labelB, r_v, red_law, tem, den, ):

    red_corr = pn.RedCorr(R_V=r_v, law=red_law)
    lineA, lineB = lime.Line(labelA, band=measurements_df), lime.Line(labelB, band=measurements_df)

    H1 = pn.RecAtom('H', 1)

    obs_ratio = measurements_df.loc[labelA, 'profile_flux']/measurements_df.loc[labelB, 'profile_flux']
    theo_ratio = H1.getEmissivity(tem=tem, den=den, wave=lineA.wavelength[0]) /H1.getEmissivity(tem=tem, den=den, wave=lineB.wavelength[0])
    red_corr.setCorr(obs_ratio/theo_ratio, wave1=lineA.wavelength[0], wave2=lineB.wavelength[0])
    print(f'cHbeta = {red_corr.cHbeta:0.3f} ({labelA}/{labelB})')
    print(f'cHbeta = {red_corr.cHbeta:0.3f} ({labelA}/{labelB})')
    red_corr.cHbetaFromEbv()

    return


specsy.theme.set_style('dark', library='bokeh')

# Load the line measurements
log_file = './sample_data/manga_lines_log.txt'
log_df = lime.load_frame(log_file)
R_V = 3.1
law = 'G03 LMC'
tem, den = 10000, 100
cHbeta = 1.235
rc = pn.RedCorr(R_V=R_V, law=law, cHbeta=cHbeta)

H1 = pn.RecAtom('H', 1)
line_labels = ['H1_4102A', 'H1_4340A', 'H1_4861A', 'H1_6563A', 'H1_9229A']
wavelength = np.empty(len(line_labels))
profile_flux = np.empty(len(line_labels))
profile_flux_err = np.empty(len(line_labels))
particles = ['H1'] * len(line_labels)

for i, label in enumerate(line_labels):
    line = lime.Line(label)
    wavelength[i] = line.wavelength[0]
    ext_i = rc.getCorr(line.wavelength[0], rel_wave=4861.25)
    profile_flux[i] = H1.getEmissivity(tem=tem, den=den, wave=line.wavelength[0]) / ext_i
    profile_flux_err[i] = profile_flux[i] * 0.05


df = pd.DataFrame(data={"wavelength": wavelength, 'profile_flux': profile_flux, 'profile_flux_err': profile_flux_err,
                        'particle': particles}, index=line_labels)


norm_line = 'H1_4861A'
cHbeta, cHbeta_err, log_ext = sy.extinction_coeff_calc(df, norm_line=norm_line, R_V=R_V, law=law, exclude_list=['H1_3889A_m', 'H1_3970A'], plot_results=False)
print(f'My result ({norm_line}): cHbeta = {cHbeta}')
sy.plotting.plots.extinction_gradient(cHbeta, cHbeta_err, log_ext)

# norm_line = 'H1_6563A'
# cHbeta, cHbeta_err, log_ext = sy.extinction_coeff_calc(df, norm_line=norm_line, R_V=R_V, law=law, exclude_list=['H1_3889A_m', 'H1_3970A'], plot_results=True)
# print(f'My result ({norm_line}): cHbeta = {cHbeta}')
#
# norm_line = 'H1_6563A'
# cHbeta, cHbeta_err, log_ext = sy.extinction_coeff_calc(df, norm_line=norm_line, R_V=R_V, law=law, exclude_list=['H1_3889A_m', 'H1_3970A'], plot_results=True, rel_Hbeta=False)
# print(f'My result ({norm_line}): cHbeta = {cHbeta}')
#
#
# pyneb_cHbeta_calc(df, 'H1_6563A', 'H1_4861A', R_V, law, tem, den)
# pyneb_cHbeta_calc(df, 'H1_9229A', 'H1_4861A', R_V, law, tem, den)
# pyneb_cHbeta_calc(df, 'H1_9229A', 'H1_6563A', R_V, law, tem, den)