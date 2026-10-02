import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from specsy.io import specsy_cfg, SpecSyError

@dataclass
class PySB99Config:
    file_path:          str
    mass_grid:          list
    evo_tracks:         np.ndarray
    minimum_wr_mass:    float
    spectra_grid_file:  str
    hires_flux:         np.ndarray
    hires_cont_flux:    np.ndarray
    hires_params:       np.ndarray
    lowmass_params:     np.ndarray
    lowmass_flux:       np.ndarray
    WN_spec_params:     np.ndarray
    WN_spectra:         np.ndarray
    WC_spec_params:     np.ndarray
    WC_spectra:         np.ndarray
    WN_spec_params_powr: np.ndarray
    WN_spectra_powr:    np.ndarray
    WC_spec_params_powr: np.ndarray
    WC_spectra_powr:    np.ndarray

    def __repr__(self):
        lines = ['PySB99Config(']
        lines.append(f'  file_path         = {self.file_path}')
        lines.append(f'  minimum_wr_mass   = {self.minimum_wr_mass}')
        lines.append(f'  spectra_grid_file = {self.spectra_grid_file}')
        lines.append(f'  mass_grid ({len(self.mass_grid)} masses) = [{self.mass_grid[0]} ... {self.mass_grid[-1]}]')
        for attr in ('evo_tracks', 'hires_flux', 'hires_cont_flux', 'hires_params',
                     'lowmass_params', 'lowmass_flux', 'WN_spec_params', 'WN_spectra',
                     'WC_spec_params', 'WC_spectra', 'WN_spec_params_powr', 'WN_spectra_powr',
                     'WC_spec_params_powr', 'WC_spectra_powr'):
            arr = getattr(self, attr)
            lines.append(f'  {attr:<22} shape={arr.shape}')
        lines.append(')')
        return '\n'.join(lines)


def _metallicity_to_key(met_val, decimals=3):

    """
    Convert a numerical metallicity into the FITS extension name notation read by
    ``_interpret_metallicity_keys`` (e.g. 0.014 -> 'Z014', 1e-05 -> 'Zem5').
    """

    # Zero has no exponent, but it round-trips through the decimal notation
    if met_val == 0:
        return 'Z' + '0' * decimals

    exponent = np.log10(met_val)

    # # Exact powers of ten below the decimal resolution use the exponent notation
    # if (met_val < 10 ** -decimals) and np.isclose(exponent, np.round(exponent)):
    #     return f'Zem{-int(np.round(exponent))}'

    # Otherwise the digits after the decimal point, keeping at least ``decimals`` of them
    digits = max(decimals, -int(np.floor(exponent)))

    return 'Z' + f'{met_val:.{digits}f}'.split('.')[1]


def pystarburst99_file_manager(cfg, M_total, IMF_exponents, IMF_mass_limits, run_speed_mode, times_out_yr,
                               SED_library='WM', POWR=False):

    # ── Unpack config into the bare variable names the engine expects ──────────
    mass_grid         = np.array(cfg.mass_grid)
    evo_tracks        = cfg.evo_tracks
    minimum_wr_mass   = cfg.minimum_wr_mass
    spectra_grid_file = cfg.spectra_grid_file
    hires_wave_grid   = np.load(cfg.file_path + 'hires_wave_grid.npy')
    empty_hires_flux  = np.full_like(hires_wave_grid, 0.0)
    hires_flux        = cfg.hires_flux
    hires_cont_flux   = cfg.hires_cont_flux
    hires_params      = cfg.hires_params
    lowmass_params    = cfg.lowmass_params
    lowmass_flux      = cfg.lowmass_flux
    WN_spec_params    = cfg.WN_spec_params
    WN_spectra        = cfg.WN_spectra
    WC_spec_params    = cfg.WC_spec_params
    WC_spectra        = cfg.WC_spectra
    WN_spec_params_powr = cfg.WN_spec_params_powr
    WN_spectra_powr     = cfg.WN_spectra_powr
    WC_spec_params_powr = cfg.WC_spec_params_powr
    WC_spectra_powr     = cfg.WC_spectra_powr

    # Time grid
    times_spectra_end = max(times_out_yr) * 1.01
    time_step         = 0.1e6                            # 0.1 Myr, matches notebook
    times_steps       = np.arange(0.0, times_spectra_end, time_step)
    times_spectra     = times_steps                      # run all steps

    # Track unpacking
    ids                 = evo_tracks[:, 0]
    ages                = evo_tracks[:, 1]
    masses              = evo_tracks[:, 2]
    luminosities        = evo_tracks[:, 3]
    temperatures        = evo_tracks[:, 4]
    H_abundances        = evo_tracks[:, 5]
    He_abundances       = evo_tracks[:, 6]
    abundance_12C       = evo_tracks[:, 7]
    abundance_14N       = evo_tracks[:, 8]
    abundance_16O       = evo_tracks[:, 9]
    core_temperatures   = evo_tracks[:, 10]
    mass_loss_rates     = evo_tracks[:, 11]

    split_factor         = len(mass_grid)
    track_ids            = np.array_split(ids, split_factor)
    track_masses         = np.array_split(masses, split_factor)
    track_ages           = np.array_split(ages, split_factor)
    track_lums           = np.array_split(luminosities, split_factor)
    track_temps          = np.array_split(temperatures, split_factor)
    track_H_abundances   = np.array_split(H_abundances, split_factor)
    track_He_abundances  = np.array_split(He_abundances, split_factor)
    track_12C_abundances = np.array_split(abundance_12C, split_factor)
    track_14N_abundances = np.array_split(abundance_14N, split_factor)
    track_16O_abundances = np.array_split(abundance_16O, split_factor)
    track_core_temps     = np.array_split(core_temperatures, split_factor)
    track_mass_loss_rates= np.array_split(mass_loss_rates, split_factor)

    for i in track_ages:
        i[0] = 1.0e-3

    for i in range(len(track_temps)):
        track_ages[i]            = np.append(track_ages[i],            [track_ages[i][-1] + 100.,            (track_ages[i][-1] + 100.) * 10**6])
        track_lums[i]            = np.append(track_lums[i],            [-20., -20.])
        track_temps[i]           = np.append(track_temps[i],           [track_temps[i][-1],           track_temps[i][-1]])
        track_masses[i]          = np.append(track_masses[i],          [track_masses[i][-1],          track_masses[i][-1]])
        track_H_abundances[i]    = np.append(track_H_abundances[i],    [track_H_abundances[i][-1],    track_H_abundances[i][-1]])
        track_He_abundances[i]   = np.append(track_He_abundances[i],   [track_He_abundances[i][-1],   track_He_abundances[i][-1]])
        track_12C_abundances[i]  = np.append(track_12C_abundances[i],  [track_12C_abundances[i][-1],  track_12C_abundances[i][-1]])
        track_14N_abundances[i]  = np.append(track_14N_abundances[i],  [track_14N_abundances[i][-1],  track_14N_abundances[i][-1]])
        track_16O_abundances[i]  = np.append(track_16O_abundances[i],  [track_16O_abundances[i][-1],  track_16O_abundances[i][-1]])
        track_core_temps[i]      = np.append(track_core_temps[i],      [track_core_temps[i][-1],      track_core_temps[i][-1]])
        track_mass_loss_rates[i] = np.append(track_mass_loss_rates[i], [track_mass_loss_rates[i][-1], track_mass_loss_rates[i][-1]])

    for i in range(len(track_mass_loss_rates)):
        for j in range(len(track_mass_loss_rates[0])):
            if track_mass_loss_rates[i][j] > -1:
                track_mass_loss_rates[i][j] = -1000

    track_masses_increasing          = np.transpose(track_masses)
    track_ages_massincind            = np.transpose(track_ages)
    track_lums_mass_incind           = np.transpose(track_lums)
    track_temps_mass_incind          = np.transpose(track_temps)
    track_H_abundances_mass_incind   = np.transpose(track_H_abundances)
    track_He_abundances_mass_incind  = np.transpose(track_He_abundances)
    track_12C_abundances_mass_incind = np.transpose(track_12C_abundances)
    track_14N_abundances_mass_incind = np.transpose(track_14N_abundances)
    track_16O_abundances_mass_incind = np.transpose(track_16O_abundances)
    track_core_temps_mass_incind     = np.transpose(track_core_temps)
    track_mass_loss_rates_mass_incind= np.transpose(track_mass_loss_rates)

    # ── Read and reform spectral grids ─────────────────────────────────────────
    spec_params, spectra = pySB99.read_spectra_grid(spectra_grid_file)
    reformed_spec_grid, spec_params_reform, spec_params_teff, spec_params_logl, spec_params_logg = pySB99.reform_spec_grid(spectra, spec_params)

    spectrum    = reformed_spec_grid[0]
    wave_grid   = spectrum[:, 0]
    empty_flux  = np.full_like(wave_grid, 0.0)

    WN_reformed_spec_grid, WN_spec_params_reform, WN_spec_params_teff = pySB99.reform_spec_grid_WR(WN_spectra, WN_spec_params)
    WC_reformed_spec_grid, WC_spec_params_reform, WC_spec_params_teff = pySB99.reform_spec_grid_WR(WC_spectra, WC_spec_params)

    WC_reformed_spec_grid_powr, WC_spec_params_reform_powr, WC_spec_params_teff_powr, WC_spec_params_radius_powr, WC_spec_params_length_powr = pySB99.reform_spec_grid_powr(WC_spectra_powr, WC_spec_params_powr)
    WN_reformed_spec_grid_powr, WN_spec_params_reform_powr, WN_spec_params_teff_powr, WN_spec_params_radius_powr, WN_spec_params_length_powr = pySB99.reform_spec_grid_powr(WN_spectra_powr, WN_spec_params_powr)

    integrated_spectra    = pySB99.integrate_spec_grid(reformed_spec_grid)
    WN_integrated_spectra = pySB99.integrate_spec_grid(WN_reformed_spec_grid)
    WC_integrated_spectra = pySB99.integrate_spec_grid(WC_reformed_spec_grid)
    WN_integrated_spectra_powr = pySB99.integrate_spec_grid(WN_reformed_spec_grid_powr)
    WC_integrated_spectra_powr = pySB99.integrate_spec_grid(WC_reformed_spec_grid_powr)

    # ── Interpolate tracks ─────────────────────────────────────────────────────
    grid_masses_adjinterp_total, grid_ages_adjinterp_total          = pySB99.interpolate_param(track_ages,            track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_lums_adjinterp_total          = pySB99.interpolate_param(track_lums,            track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_temps_adjinterp_total         = pySB99.interpolate_param(track_temps,           track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_H_abundances_adjinterp_total  = pySB99.interpolate_param(track_H_abundances,    track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_He_abundances_adjinterp_total = pySB99.interpolate_param(track_He_abundances,   track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_12C_abundances_adjinterp_total= pySB99.interpolate_param(track_12C_abundances,  track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_14N_abundances_adjinterp_total= pySB99.interpolate_param(track_14N_abundances,  track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_16O_abundances_adjinterp_total= pySB99.interpolate_param(track_16O_abundances,  track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_core_temps_adjinterp_total    = pySB99.interpolate_param(track_core_temps,      track_masses, run_speed_mode)
    grid_masses_adjinterp_total, grid_mass_loss_rates_adjinterp_total = pySB99.interpolate_param(track_mass_loss_rates,track_masses, run_speed_mode)

    grid_masses          = pySB99.rearrange_grid_array(grid_masses_adjinterp_total)
    grid_ages            = pySB99.rearrange_grid_array(grid_ages_adjinterp_total)
    grid_lums            = pySB99.rearrange_grid_array(grid_lums_adjinterp_total)
    grid_temps           = pySB99.rearrange_grid_array(grid_temps_adjinterp_total)
    grid_H_abundances    = pySB99.rearrange_grid_array(grid_H_abundances_adjinterp_total)
    grid_He_abundances   = pySB99.rearrange_grid_array(grid_He_abundances_adjinterp_total)
    grid_12C_abundances  = pySB99.rearrange_grid_array(grid_12C_abundances_adjinterp_total)
    grid_14N_abundances  = pySB99.rearrange_grid_array(grid_14N_abundances_adjinterp_total)
    grid_16O_abundances  = pySB99.rearrange_grid_array(grid_16O_abundances_adjinterp_total)
    grid_core_temps      = pySB99.rearrange_grid_array(grid_core_temps_adjinterp_total)
    grid_mass_loss_rates = pySB99.rearrange_grid_array(grid_mass_loss_rates_adjinterp_total)

    # ── Low-mass and high-res spectral setup ───────────────────────────────────
    lowmass_teffs = lowmass_params[:, 0]
    lowmass_loggs = lowmass_params[:, 1]
    lowmass_wave_grid = wave_grid
    lowmass_spec = [np.column_stack((lowmass_wave_grid, lowmass_flux[i])) for i in range(len(lowmass_flux))]
    lowmass_int_spec = pySB99.integrate_spec_grid(lowmass_spec)

    hires_teffs = hires_params[:, 0]
    hires_loggs = hires_params[:, 1]
    hires_spec = []
    hires_cont = []
    for i in range(len(hires_flux)):
        hires_spec.append(np.column_stack((hires_wave_grid, hires_flux[i])))
        hires_cont.append(np.column_stack((hires_wave_grid, hires_cont_flux[i])))
    hires_int_spec = pySB99.integrate_spec_grid(hires_spec)

    WN_spec_params_logl = np.full_like(WN_spec_params_teff, 1.)
    WC_spec_params_logl = np.full_like(WC_spec_params_teff, 1.)
    lm_spec_params_logl = np.full_like(lowmass_teffs, 1.)

    # ── Initial timestep mass index ────────────────────────────────────────────
    timestep_mass_ind = pySB99.get_timestep_0_ind(timestep=times_steps[0], grid_ages=grid_ages,
                                                  grid_masses=grid_masses, times_steps=times_steps,
                                                  IMF_mass_limits=IMF_mass_limits)

    # ── Main time loop ─────────────────────────────────────────────────────────
    population_flux_total_iterations = []

    for i in range(len(times_steps)):
        timestep_ex = times_steps[i]
        (timestep_ages_final, timestep_temps_final, timestep_lums_final, timestep_masses_final,
         timestep_H_abundances_final, timestep_loggs_final, timestep_mass_loss_rates_final,
         timestep_12C_abundances_final, timestep_14N_abundances_final, timestep_16O_abundances_final,
         timestep_mass_test, timestep_cnr, timestep_coher) = pySB99.get_timestep_params(timestep=timestep_ex,
                                                                                        timestep_mass_ind=timestep_mass_ind,
                                                                                        grid_ages=grid_ages,
                                                                                        grid_masses=grid_masses,
                                                                                        grid_temps=grid_temps,
                                                                                        grid_lums=grid_lums,
                                                                                        grid_H_abundances=grid_H_abundances,
                                                                                        grid_He_abundances=grid_He_abundances,
                                                                                        grid_12C_abundances=grid_12C_abundances,
                                                                                        grid_14N_abundances=grid_14N_abundances,
                                                                                        grid_16O_abundances=grid_16O_abundances,
                                                                                        grid_core_temps=grid_core_temps,
                                                                                        grid_mass_loss_rates=grid_mass_loss_rates)

        if i == 0:
            initial_masses = timestep_masses_final
            No_stars, c_masses, dens, xmhigh, xmlow = pySB99.calc_Nostars(IMF_masses=timestep_masses_final,
                                                                          IMF_exponents=IMF_exponents,
                                                                          IMF_mass_limits=IMF_mass_limits,
                                                                          M_total=M_total)

        timestep_teffs_final = 10 ** timestep_temps_final
        timestep_radii_final = pySB99.compute_radii(timestep_temps_final, timestep_lums_final)
        specsyn_teffs, specsyn_loggs, specsyn_radii, specsyn_cotests, specsyn_bbfluxes = pySB99.get_specsyn_params(
            timestep_temps_final, timestep_masses_final, timestep_lums_final)

        assigned_integrated_spectra, assigned_spec_teff, assigned_spectra, population_choice, assigned_spec_logl, assigned_spec_logg = pySB99.assign_spectra_to_grid_WR(
            wave_grid, empty_flux, minimum_wr_mass, WC_spec_params_reform, WN_spec_params_reform,
            WN_spec_params_teff, WN_spec_params_logl, WN_reformed_spec_grid, WN_integrated_spectra,
            WC_spec_params_teff, WC_spec_params_logl, WC_reformed_spec_grid, WC_integrated_spectra,
            No_stars, reformed_spec_grid, integrated_spectra, lowmass_params, timestep_loggs_final,
            lowmass_teffs, lowmass_loggs, lm_spec_params_logl, lowmass_spec, lowmass_int_spec,
            timestep_temps_final, spec_params_reform, spec_params_teff, spec_params_logl,
            timestep_lums_final, timestep_H_abundances_final, spec_params_logg,
            timestep_masses_final, initial_masses, specsyn_cotests, specsyn_loggs,
            timestep_cnr, timestep_coher)

        population_flux, _ = pySB99.specsyn(assigned_integrated_spectra, specsyn_bbfluxes, assigned_spectra, specsyn_radii, No_stars)

        population_ion_flux = pySB99.ionise(wave_grid, population_flux, 2)
        population_continuum = pySB99.continuum(population_ion_flux[1][0])
        continuum_resampled = np.interp(wave_grid, population_continuum[:, 0], population_continuum[:, 1])
        population_flux_total = population_flux + continuum_resampled

        population_flux_total_iterations.append(population_flux_total)

    # ── Extract at requested ages ──────────────────────────────────────────────
    fluxes = []
    for t_yr in times_out_yr:
        idx = int(np.argmin(np.abs(times_steps - t_yr)))
        fluxes.append(np.array(population_flux_total_iterations[idx]))

    return wave_grid, fluxes


def load_pySB99_config(Z: str, base_path: Path, rot: bool = False) -> PySB99Config:

    z_cfg = specsy_cfg['stellar']['ssp']['pystarburst99'].get(Z)
    if z_cfg is None:
        valid = list(specsy_cfg['stellar']['ssp']['pystarburst99'].keys())
        raise ValueError(f"Unknown metallicity '{Z}'. Valid options: {valid}")

    fp = base_path / z_cfg['folder']
    tracks_key = 'tracks_v40' if rot else 'tracks_v00'
    s = z_cfg['wr_suffix']

    ld = lambda f: np.load(fp / f)
    lp = lambda f: np.load(fp / f, allow_pickle=True)

    return PySB99Config(file_path           = str(fp) + '/',
                        mass_grid           = list(z_cfg['mass_grid']),
                        evo_tracks          = ld(z_cfg[tracks_key]),
                        minimum_wr_mass     = float(z_cfg['minimum_wr_mass']),
                        spectra_grid_file   = str(fp / z_cfg['wm_ob']),
                        hires_flux          = ld(z_cfg['ifa_line']),
                        hires_cont_flux     = ld(z_cfg['ifa_cont']),
                        hires_params        = ld(z_cfg['ifa_params']),
                        lowmass_params      = ld(z_cfg['lm_params']),
                        lowmass_flux        = ld(z_cfg['lm_flux']),
                        WN_spec_params      = ld(f'WN_spec_params_cmfgen_{s}.npy'),
                        WN_spectra          = lp(f'WN_spectra_cmfgen_{s}.npy'),
                        WC_spec_params      = ld(f'WC_spec_params_cmfgen_{s}.npy'),
                        WC_spectra          = lp(f'WC_spectra_cmfgen_{s}.npy'),
                        WN_spec_params_powr = ld(f'WN_spec_params_powr_{s}.npy'),
                        WN_spectra_powr     = lp(f'WN_spectra_powr_{s}.npy'),
                        WC_spec_params_powr = ld(f'WC_spec_params_powr_{s}.npy'),
                        WC_spectra_powr     = lp(f'WC_spectra_powr_{s}.npy'),)





def _parse_pysb99_folder(folder_list, load_spectrum, sun_luminosity=3.83e33):
    """Parse one or more pySB99 output folders into a list of Binary objects."""

    binary_list = []
    for folder in folder_list:
        folder = Path(folder)

        if folder.is_dir():

            # Read metadata
            input_path = folder / 'input.txt'
            if not input_path.exists():
                raise SpecSyError(f'Missing input.txt in {folder}')
            meta = parse_pysb99_input_file(input_path)

            # Read wavelength grid
            wave_path = folder / 'spectrum_wavelength.txt'
            if not wave_path.exists():
                raise SpecSyError(f'Missing SED_wavelength.txt in {folder}')
            wave_arr = np.loadtxt(wave_path)

            # Read time axis (in years, linear)
            times_path = folder / 'timesteps.txt'
            if not times_path.exists():
                raise SpecSyError(f'Missing timesteps.txt in {folder}')

            # Check if the interpolation first value is 0
            age_arr = np.loadtxt(times_path)
            idx_0 = 0 if age_arr[0] != 0 else 1
            log_ages = np.log10(age_arr[idx_0:])

            # Read flux array — shape (n_wave, n_times), one column per age
            if load_spectrum:
                flux_path = folder / 'pySB_hires_spectrum.npy'
                if not flux_path.exists():
                    raise SpecSyError(f'Missing pySB_hires_spectrum.npy in {folder}')
                flux_matrix = np.load(flux_path)  # shape (n_wave, n_times), no conversion needed
                flux_matrix = np.power(10, flux_matrix[idx_0:, :])/sun_luminosity

            for i, log_age in enumerate(log_ages):
                binary_list.append(Binary(
                    source='pystarburst99',
                    library=meta['library'],
                    metallicity=meta['metallicity'],
                    alpha=0.0,
                    imf='kroupa',
                    age=log_age,
                    fpath=str(folder),
                    wavelength=wave_arr if load_spectrum else None,
                    flux=flux_matrix[i, :] if load_spectrum else None,
                    low_imf_exp=meta.get('low_imf_exp'),
                    high_imf_exp=meta.get('high_imf_exp'),
                    low_imf_mass=meta.get('low_imf_mass'),
                    high_imf_mass=meta.get('high_imf_mass'),
                ))

    return binary_list

