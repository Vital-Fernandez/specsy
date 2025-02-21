import lime
import specsy as sy
import pyneb as pn

# State the data files
obsFitsFile = '/home/vital/PycharmProjects/lime/examples/sample_data/spectra/gp121903_osiris.fits'
lineBandsFile = '/home/vital/PycharmProjects/lime/examples/sample_data/osiris_bands.txt'
cfgFile = '/home/vital/PycharmProjects/lime/examples/sample_data/osiris.toml'

# Load configuration
obs_cfg = lime.load_cfg(cfgFile)
z_obj = obs_cfg['sample_data']['z_array'][2]
norm_flux = obs_cfg['sample_data']['norm_flux']

# Declare LiMe spectrum
gp_spec = lime.Spectrum.from_file(obsFitsFile, instrument='osiris', redshift=z_obj, norm_flux=norm_flux)

# Fit the continuum
gp_spec.fit.continuum(degree_list=[3, 6, 6], emis_threshold=[3, 2, 1.5], plot_steps=False)

# Find lines
match_bands = gp_spec.line_detection(lineBandsFile, sigma_threshold=3, plot_steps=False)

# Measure the emission lines
gp_spec.fit.frame(match_bands, obs_cfg, id_conf_prefix='gp121903', line_detection=True)
# gp_spec.plot.spectrum()

# Save the results
log_file = './sample_data/gp121903_log.txt'
gp_spec.save_frame(log_file)

# Reload and normalized
log_df = lime.load_frame(log_file)
lime.normalize_fluxes(log_df, norm_list='H1_4861A')
lime.save_frame('./sample_data/gp121903_log_norm.txt', log_df)

cHbeta, cHbeta_err = sy.cHbeta_from_log(log_df, ref_line='H1_4861A', lines_ignore=['H1_3889A_m', 'H1_3970A'], show_plot=False)
# cHbeta, cHbeta_err = sy.cHbeta_from_log(log_df, ref_line='H1_4861A', lines_ignore=['H1_3889A_m', 'H1_3970A'], flux_entry='line',
#                                         show_plot=True)

sy.extinction.reddening_correction(cHbeta, cHbeta_err, log_df, norm_wavelength=4861.0, flux_column='line_flux')
lime.save_frame('./sample_data/gp121903_log_norm_int.txt', log_df)

# Compare with Hbeta
rc = pn.RedCorr(cHbeta=cHbeta, R_V=3.1, law='G03 LMC')
sy_corr = log_df.loc["H1_6563A", "line_int"]/log_df.loc["H1_6563A", "line_flux"]

print(f'\nPyneb correction: {rc.getCorr(6563)}')
print(f'LiMe correction: {1/sy_corr}')