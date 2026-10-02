from pathlib import Path
import lime
from specsy.io import load_trace
from specsy.models.ssp import StellarBinaries
from specsy.models.stellar_synthesis import SSP_sampler, plot_best_fit, load_ionization_frame, summary_ssp
from specsy.plotting.plots import ssp_fitting_infograph

# Stellar grid properties
bpass_pname = Path('/home/vital/PycharmProjects/specsy/examples/stellar_synthesis/SESAMME_BAPSS-AP_v3.fits')
bpass_grid = StellarBinaries.from_fits(bpass_pname)

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

        # Grid nodes on the observed wavelengths
        obj_ssp = bpass_grid.to_stellar_binaries(wmin=wave_arr[0], wmax=wave_arr[-1])

        # Sampling
        sampler = SSP_sampler(spec, obj_ssp, red_law='CCM89', r_v=3.1)
        sampler.prepare_inputs(wave_arr, flux_arr, err_arr, mask_arr, p0_bootstrap=True, rng_seed=42, add_nebular=add_nebular, ion_frame=ion_frame)
        trace = sampler.sample(**sampler_cfg['W512_x2'], truth_dict=synth_dict)
        sampler.summary()
        sampler.save_trace(trace_pname)

        # plot_best_fit(trace_data, title=synth_spec, savefile=results_folder/f'{synth_spec}_diagnostic.png')
        trace_data = load_trace(trace_pname)
        results_df = summary_ssp(trace_data, truth_dict=synth_dict)
        ssp_fitting_infograph(trace_data, title=synth_spec)

        # trace_data = load_trace(trace_pname, drop_stuck=True, k=1.5, dchi2=20.0, seed=42)
        # results_df = summary_ssp(trace_data, truth_dict=synth_dict)
        # plot_best_fit(trace_data, title=synth_spec)