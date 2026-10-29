from pathlib import Path
import lime
from specsy.io import load_trace
from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import SSP_sampler, plot_best_fit, load_ionization_frame, summary_ssp
from specsy.plotting.plots import ssp_fitting_infograph, ssp_fitting_infograph_bokeh

# Stellar grid properties
bpass_pname = Path('/home/vital/PycharmProjects/specsy/examples/stellar_synthesis/SESAMME_BAPSS-AP_v3.fits')
bpass_grid = StellarBinaries.from_fits(bpass_pname)


# 0.00001, 0.00010, 0.00100, 0.00200, 0.00300, 0.00400, 0.00600
# 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9

bpass_grid.plot_spectra(age=[6.5], metallicity=[0.00001, 0.00010, 0.00100, 0.00200, 0.00300], fig_cfg={'legend.fontsize':16})

# Nebular continuum
ion_pname = Path('/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt')
ion_frame = load_ionization_frame(ion_pname)

# Configuration and results folder (mock spectra are read from and results saved to the same place)
synth_file = Path('./synth_sample.toml')
mock_obs_cfg = lime.load_cfg(synth_file)

# Sampling configuration
sampling_fname = Path('./fitting_configuration_entries.toml')
sampler_cfg = lime.load_cfg(sampling_fname)

# Outputs configuration
results_folder = Path('/home/vital/Astrodata/SpecSy_tests/SSP_sampling')

# Fit all the mock observations in the configuration file, one after another
for i, (synth_spec, synth_dict) in enumerate(mock_obs_cfg.items()):

    if i == 2:

        print(f'\n{synth_spec}')
        spec_pname = results_folder / f'{synth_spec}.txt'
        trace_pname = results_folder / f'{synth_spec}_trace'
        add_nebular = synth_dict.get('add_nebular', False)

        # Load the mock spectrum
        spec = lime.Spectrum.from_file(spec_pname, instrument='text')
        spec.unit_conversion(distance=synth_dict['distance'])
        wave_arr, flux_arr, err_arr, mask_arr = spec.retrieve.spectrum(return_arrays=True)
        mask_arr = ~mask_arr

        # # Grid nodes on the observed wavelengths
        # obj_ssp = bpass_grid.to_stellar_binaries(wmin=wave_arr[0], wmax=wave_arr[-1])
        #
        # # Sampling
        # sampler = SSP_sampler(spec, obj_ssp, red_law='CCM89', r_v=3.1)
        # sampler.prepare_inputs(wave_arr, flux_arr, err_arr, mask_arr, p0_bootstrap=True, rng_seed=42, add_nebular=add_nebular, ion_frame=ion_frame)
        # trace = sampler.sample(**sampler_cfg['W512_x2'], truth_dict=synth_dict)
        # sampler.summary()
        # sampler.save_trace(trace_pname)
        #
        # trace_data = load_trace(trace_pname)
        # results_df = summary_ssp(trace_data, truth_dict=synth_dict)
        # ssp_fitting_infograph(trace_data, title=synth_spec)

        # trace_data = load_trace(trace_pname, drop_stuck=True, k=1.5, dchi2=20.0, seed=42)
        # results_df = summary_ssp(trace_data, truth_dict=synth_dict)
        # plot_best_fit(trace_data, title=synth_spec)


        trace_data = load_trace(trace_pname)
        results_df = summary_ssp(trace_data, truth_dict=synth_dict)

        residual_df = lime.lines_frame(line_list=['C3_1176A', 'C4_1548A', 'N5_1239A', 'He2_1640A'])
        residual_df['w3'] -= 20
        residual_df['w4'] += 20

        ssp_fitting_infograph(trace_data, title=synth_spec, solar_Z=False, show_grid=True, residual_df=residual_df)
        # ssp_fitting_infograph_bokeh(trace_data, title=synth_spec, solar_Z=False)
        # ssp_fitting_infograph(trace_data, title=synth_spec, residual_df=residual_df, solar_Z=False)
        # ssp_fitting_infograph_bokeh(trace_data, title=synth_spec, residual_df=residual_df, solar_Z=True)