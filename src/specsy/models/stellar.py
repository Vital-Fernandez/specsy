from logging import getLogger
from pathlib import Path
from dataclasses import dataclass, field

from astropy.io import fits
from astropy.table import Table

import re
import numpy as np
import pandas as pd

from specsy.io import specsy_cfg, SpecSyError
from lime.workflow import spectrum_resampling
from matplotlib import pyplot as plt, rc_context
from lime.plotting.plots import _NO_FIG, save_close_fig_swicth
import specsy.models.tools_pySB99 as pySB99


# Logger
_logger = getLogger('SpecSy')

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


# def _parse_pysb99_input_file(input_path):
#     """Parse pySB99 input.txt into a metadata dict."""
#
#     params = {}
#     with open(input_path) as f:
#         for line in f:
#             line = line.strip()
#             if line.startswith('Metallicity input'):
#                 z_label = line.split('=')[1].strip()
#                 params['metallicity'] = SB99_Z_MAP[z_label]
#             elif line.startswith('SED library input'):
#                 params['library'] = line.split('=')[1].strip()
#             elif line.startswith('Rotation input'):
#                 params['rotation'] = line.split('=')[1].strip() == 'True'
#             elif line.startswith('IMF_exponents'):
#                 raw = line.split('=')[1].strip()
#                 exponents = [float(x) for x in raw.strip('[]').split(',')]
#                 params['low_imf_exp'] = exponents[0]
#                 params['high_imf_exp'] = exponents[1] if len(exponents) > 1 else None
#             elif line.startswith('IMF_mass_limits'):
#                 raw = line.split('=')[1].strip()
#                 limits = [float(x) for x in raw.strip('()[]').split(',')]
#                 params['low_imf_mass'] = limits[0]
#                 params['high_imf_mass'] = limits[-1]
#
#     return params

def parse_pysb99_input_file(input_path):
    """Parse pySB99 input.txt into a metadata dict."""

    params = {}
    with open(input_path) as f:
        for line in f:
            line = line.strip()
            if not line or '=' not in line:
                continue

            key, _, value = line.partition('=')
            key = key.strip()
            value = value.strip()

            if key == 'Metallicity input':
                z_label = value
                params['metallicity'] = specsy_cfg['stellar']['ssp']['pystarburst99']['z_keys'][z_label]

            elif key == 'spectra library input':
                params['library'] = value

            elif key == 'Rotation input':
                params['rotation'] = value == 'True'

            elif key == 'IMF_exponents':
                exponents = [float(x) for x in value.strip('[]').split(',')]
                params['low_imf_exp'] = exponents[0]
                params['high_imf_exp'] = exponents[1] if len(exponents) > 1 else None

            elif key == 'IMF_mass_limits':
                limits = [float(x) for x in value.strip('()[]').split(',')]
                params['low_imf_mass'] = limits[0]
                params['high_imf_mass'] = limits[-1]

    return params


def _parse_pysb99_folder(folder_list, load_spectrum):
    """Parse one or more pySB99 output folders into a list of Binary objects."""

    binary_list = []
    for folder in folder_list:
        folder = Path(folder)

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
        log_ages = np.log10(np.loadtxt(times_path))

        # Read flux array — shape (n_wave, n_times), one column per age
        if load_spectrum:
            flux_path = folder / 'pySB_hires_spectrum.npy'
            if not flux_path.exists():
                raise SpecSyError(f'Missing pySB_hires_spectrum.npy in {folder}')
            flux_matrix = np.load(flux_path)  # shape (n_wave, n_times), no conversion needed

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


def bpass_file_manager(fname, load_spectrum, metal_map=specsy_cfg['stellar']['ssp']['BPASS']['z_keys']):

    """Parse: spectra-bin-imf135_300.LIBRARY.zXXX.aYYY.dat"""
    name = fname.stem
    pattern = r'spectra-bin-(imf[\w]+)\.([\w]+)\.(z[\w]+)\.(a[+-]\d+)'
    m = re.match(pattern, name)

    if not m:
        raise ValueError(f"Unrecognised filename format for BPASS files: {fname}")

    imf_str, library, z_str, a_str = m.groups()

    # Solar metallicity
    if z_str in metal_map:
        metallicity = metal_map[z_str]
    else:
        metallicity = metal_map[z_str[:1].upper() + z_str[1:]]

    # Alpha enhancement
    alpha = int(a_str.replace('a', '')) / 100

    # Parse IMF string: imf135_300 -> low_imf_exp=1.35, high_imf_mass=300.0
    imf_pattern = r'imf(\d+)_(\d+)'
    imf_match = re.match(imf_pattern, imf_str)
    if imf_match:
        low_imf_exp = float(imf_match.group(1)) / 100
        high_imf_mass = float(imf_match.group(2))
    else:
        low_imf_exp = None
        high_imf_mass = None

    params = dict(library=library, metallicity=metallicity, alpha=alpha, imf=imf_str,
                  low_imf_exp=low_imf_exp, high_imf_exp=None,
                  low_imf_mass=None, high_imf_mass=high_imf_mass)

    data_arr = np.loadtxt(fname) if load_spectrum else None
    wave_arr, flux_arr = (data_arr[:, 0], data_arr[:, 1:]) if data_arr is not None else (None, None)

    return params, wave_arr, flux_arr


def parse_binaries(source, file_list, load_spectrum=True, sed_type='stellar_and_nebular'):

    match source:
        case 'BPASS':

            BPASS_LOG_AGES = specsy_cfg['stellar']['ssp']['BPASS']['log_age_arr']

            binary_list = []
            for fname in file_list:

                bin_params, wave_arr, flux_matrix = bpass_file_manager(fname, load_spectrum)

                bin_age_list = [None] * len(BPASS_LOG_AGES)
                for i, age in enumerate(BPASS_LOG_AGES):
                    idx_age = np.searchsorted(BPASS_LOG_AGES, age)
                    bin_age_list[i] = Binary(source='bpass',
                                             age=age,
                                             fpath=str(fname),
                                             wavelength=wave_arr,
                                             flux=flux_matrix[:, idx_age] if load_spectrum else None,
                                             **bin_params)
                binary_list += bin_age_list

            return binary_list

        case 'pystarburst99':
            return _parse_pysb99_folder(file_list, load_spectrum)

        case _:
            raise SpecSyError(f'Input binary source "{source}" is not supported')


@dataclass
class Binary:
    source: str
    library: str
    metallicity: float
    alpha: float
    age: float
    fpath: str
    imf: str | None = field(default=None)
    flux: np.ndarray | None = field(default=None, repr=False)
    wavelength: np.ndarray | None = field(default=None, repr=False)
    low_imf_exp: float | None = field(default=None)
    high_imf_exp: float | None = field(default=None)
    low_imf_mass: float | None = field(default=None)
    high_imf_mass: float | None = field(default=None)

    @property
    def is_loaded(self):
        return self.flux is not None

    def clear_spectrum(self):
        self.flux = None
        self.wavelength = None


class StellarBinaries:

    def __init__(self, binary_list):

        # Attributes
        self._binaries = {}
        self.uniform_wmin = None
        self.uniform_wmax = None
        self.uniform_deltalambda = None
        self.source = None
        self.library = None

        # Assign an index to each binary
        for entry in binary_list:
            self._binaries[(entry.source, entry.library, entry.imf, entry.age, entry.metallicity, entry.alpha)] = entry

        # Unique metallicities and ages across all binaries
        self.metallicities = np.unique([b.metallicity for b in self._binaries.values()])
        self.ages = np.unique([b.age for b in self._binaries.values()])

        # Check wavelength uniformity across loaded binaries
        self._check_dispersion_uniformity()

        return

    def _check_dispersion_uniformity(self):

        loaded = [b for b in self._binaries.values() if b.is_loaded]
        if len(loaded) == 0:
            return

        # Wavelength uniformity
        wmins = np.array([b.wavelength.min() for b in loaded])
        wmaxs = np.array([b.wavelength.max() for b in loaded])
        deltas = np.array([np.mean(np.diff(b.wavelength)) for b in loaded])

        self.uniform_wmin = float(wmins[0]) if np.all(np.isclose(wmins, wmins[0])) else None
        self.uniform_wmax = float(wmaxs[0]) if np.all(np.isclose(wmaxs, wmaxs[0])) else None
        self.uniform_deltalambda = float(deltas[0]) if np.all(np.isclose(deltas, deltas[0])) else None

        non_uniform = [k for k, v in [('wmin', self.uniform_wmin), ('wmax', self.uniform_wmax),
                                      ('delta_lambda', self.uniform_deltalambda)] if v is None]
        if non_uniform:
            _logger.warning(f'Wavelength grid is not uniform across loaded binaries for: {non_uniform}')

        # Source and library assignment
        unique_sources = sorted({b.source for b in self._binaries.values()})
        unique_libraries = sorted({b.library for b in self._binaries.values()})

        self.source = unique_sources[0] if len(unique_sources) == 1 else unique_sources
        self.library = unique_libraries[0] if len(unique_libraries) == 1 else unique_libraries

        if len(unique_sources) > 1:
            _logger.warning(f'StellarBinaries contains multiple sources: {unique_sources}')
        if len(unique_libraries) > 1:
            _logger.warning(f'StellarBinaries contains multiple libraries: {unique_libraries}')

        return

    @classmethod
    def from_source(cls, source_name, source_folder, root_folder=None, load_spectra=False,
                    ext=None, file_list=None, sed_type='stellar_and_nebular'):

        if source_name not in specsy_cfg['stellar']['source_list']:
            msg = f'Stellar binaries source is not recognized. Supported: {specsy_cfg["stellar"]["source_list"]}'
            raise SpecSyError(msg)

        if root_folder is None:
            source_folder = Path(source_folder).resolve()
        else:
            source_folder = Path(source_folder)
            if source_folder.is_absolute():
                source_folder = source_folder.relative_to("/")
            source_folder = (Path(root_folder) / source_folder).resolve()

        if not source_folder.exists():
            raise SpecSyError(f'Folder not found: {source_folder}')

        # SB99: source_folder is a root containing one subfolder per run
        if source_name == 'SB99':
            if file_list is not None:
                folder_list = [source_folder / f for f in file_list]
                missing = [f for f in folder_list if not f.exists()]
                if missing:
                    raise SpecSyError(f'Folders not found: {missing}')
            else:
                folder_list = [f for f in source_folder.iterdir() if f.is_dir()]
                if len(folder_list) == 0:
                    raise SpecSyError(f'No subfolders found in {source_folder}')

            binary_list = parse_binaries(source_name, folder_list, load_spectra, sed_type)
            return cls(binary_list)

        # BPASS: original file-based logic
        if file_list is not None:
            list_fname = [source_folder / f for f in file_list]
            missing = [f for f in list_fname if not f.exists()]
            if missing:
                raise SpecSyError(f'Files not found in {source_folder}: {missing}')
        else:
            list_fname = list(source_folder.glob(f'*{"" if ext is None else ext}'))
            if len(list_fname) == 0:
                raise SpecSyError(f'No files found in {source_folder} with extension {ext}')

        binary_list = parse_binaries(source_name, list_fname, load_spectra)
        return cls(binary_list)


    def __repr__(self):

        # Sources and libraries from attributes
        sources = [self.source] if isinstance(self.source, str) else self.source
        libraries = [self.library] if isinstance(self.library, str) else self.library
        sources = sources if sources is not None else []
        libraries = libraries if libraries is not None else []

        n_loaded = sum(1 for b in self._binaries.values() if b.is_loaded)

        # Wavelength uniformity # TODO add units attribute
        wmin_str = f'{self.uniform_wmin:.1f} Å' if self.uniform_wmin is not None else 'non-uniform'
        wmax_str = f'{self.uniform_wmax:.1f} Å' if self.uniform_wmax is not None else 'non-uniform'
        delta_str = f'{self.uniform_deltalambda:.3f} Å' if self.uniform_deltalambda is not None else 'non-uniform'

        # Ages and metallicities
        age_str = ', '.join(f'{a:.1f}' for a in self.ages)
        met_str = ', '.join(f'{m:.5f}' for m in self.metallicities)

        return (
            f'\n{"=" * 60}\n'
            f'  StellarBinaries\n'
            f'{"=" * 60}\n'
            f'  Source(s)      : {", ".join(sources)}\n'
            f'  Libraries      : {", ".join(libraries)}\n'
            f'  Total stellar spectra : {len(self._binaries)} ({n_loaded} loaded)\n'
            f'{"-" * 60}\n'
            f'  Wavelength grid\n'
            f'    wmin         : {wmin_str}\n'
            f'    wmax         : {wmax_str}\n'
            f'    delta lambda : {delta_str}\n'
            f'{"-" * 60}\n'
            f'  Ages (log yr)  [{len(self.ages)}]: {age_str}\n'
            f'{"-" * 60}\n'
            f'  Metallicities  [{len(self.metallicities)}]: {met_str}\n'
            f'{"=" * 60}\n'
        )

    def __getitem__(self, key):
        """
        Flexible indexing by any combination of (source, library, imf, age, metallicity, alpha):
            sb['BPASS']                              -> list by source
            sb['BPASS', 'C3K']                       -> list by source + library
            sb['BPASS', 'C3K', 'imf135_300']         -> list by source + library + imf
            sb['BPASS', 'C3K', 'imf135_300', 8.0]    -> list by source + library + imf + age
        """
        if not isinstance(key, tuple):
            key = (key,)
        return self._filter(*key)

    def _filter(self, source=None, library=None, imf=None, age=None, metallicity=None, alpha=None):
        results = [b for (src, lib, imf_key, age_key, met, alp), b in self._binaries.items()
                   if (source is None or src == source)
                   and (library is None or lib == library)
                   and (imf is None or imf_key == imf)
                   and (age is None or np.isclose(age_key, age))
                   and (metallicity is None or np.isclose(met, metallicity))
                   and (alpha is None or np.isclose(alp, alpha))]

        return results if len(results) > 0 else None

    def get_spectrum(self, age, metallicity, library=None, alpha=None, imf=None):

        # Filter matching binaries
        results = self._filter(library=library, imf=imf, age=age, metallicity=metallicity, alpha=alpha)

        # Security checks
        if results is None:
            raise SpecSyError(f'No binary found for age={age}, metallicity={metallicity}.')

        if len(results) > 1:
            libs = [b.library for b in results]
            alphas = [b.alpha for b in results]
            imfs = [b.imf for b in results]
            raise SpecSyError(f'Multiple binaries ({len(results)}) match age={age}, metallicity={metallicity}. '
                              f'Please refine with: library={libs}, alpha={alphas}, imf={imfs}.')

        b = results[0]
        if not b.is_loaded:
            raise SpecSyError(f'Binary with age={age}, metallicity={metallicity} is not loaded.')

        return b.wavelength, b.flux

    def to_fits(self, fname, libraries=None, wmin=None, wmax=None, ext_map=None, disp_intvl=None, pixel_width=None,
                pixel_number=None, constant_pixel_width=True):

        # Check output folder exists
        output_path = Path(fname)
        if not output_path.parent.exists():
            raise SpecSyError(f'Output folder not found: {output_path.parent}')

        # Filter binaries by requested libraries
        all_libraries = sorted({b.library for b in self._binaries.values()})
        plot_libraries = all_libraries if libraries is None else (
            [libraries] if isinstance(libraries, str) else list(libraries))

        # Check all requested libraries exist
        unrecognised = set(plot_libraries) - set(all_libraries)
        if unrecognised:
            _logger.warning(f'Libraries not found in StellarBinaries: {unrecognised}. '
                            f'Available libraries: {all_libraries}')
        plot_libraries = [lib for lib in plot_libraries if lib in all_libraries]
        if len(plot_libraries) == 0:
            raise SpecSyError(f'No valid libraries to export. Available libraries: {all_libraries}')

        # Check spectra are loaded
        for b in self._binaries.values():
            if b.library in plot_libraries and not b.is_loaded:
                raise SpecSyError(f'Spectra are not loaded. Call load() before exporting to fits.')

        met_groups = {}
        for b in self._binaries.values():
            if b.library in plot_libraries:
                met_groups.setdefault(b.metallicity, []).append(b)

        # Dictionary for mapping the fits header metallicity values according to the dataset
        if ext_map is None:
            if self.source == 'bpass':
                BPASS_Z_KEYS = specsy_cfg['stellar']['ssp']['BPASS']['z_keys']
                ext_map = {value: key for key, value in BPASS_Z_KEYS.items()}
            elif self.source == 'pystarburst99':
                pySB99_Z_KEYS = specsy_cfg['stellar']['ssp']['pystarburst99']['z_keys']
                ext_map = {value: key for key, value in pySB99_Z_KEYS.items()}
            else:
                raise SpecSyError(f'Binaries source "{self.source}" does not have a metallicity header map')

        # Check if rebinning is requested
        rebin_args = {'disp_intvl': disp_intvl, 'pixel_width': pixel_width, 'pixel_number': pixel_number}
        active_rebin = [k for k, v in rebin_args.items() if v is not None]
        if len(active_rebin) > 1:
            raise SpecSyError(f'Rebinning arguments {active_rebin} are mutually exclusive. Please provide only one.')
        do_rebin = len(active_rebin) == 1

        # Check if alpha and imf are uniform across selected binaries
        selected = [b for b in self._binaries.values() if b.library in plot_libraries]
        unique_alphas = {b.alpha for b in selected}

        # IMF uniformity checks
        unique_imfs = {b.imf for b in selected}
        unique_low_exp = {b.low_imf_exp for b in selected}
        unique_high_exp = {b.high_imf_exp for b in selected}
        unique_low_mass = {b.low_imf_mass for b in selected}
        unique_high_mass = {b.high_imf_mass for b in selected}

        # Build FITS primary header
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['SOURCE'] = ', '.join({b.source for b in selected})
        primary_hdu.header['LIBS'] = ', '.join(plot_libraries)

        if len(unique_alphas) == 1:
            primary_hdu.header['ALPHA'] = float(next(iter(unique_alphas)))
        if len(unique_imfs) == 1:
            v = next(iter(unique_imfs))
            if v is not None:
                primary_hdu.header['IMF'] = v
        if len(unique_low_exp) == 1:
            v = next(iter(unique_low_exp))
            if v is not None:
                primary_hdu.header['LOW_IMFE'] = float(v)
        if len(unique_high_exp) == 1:
            v = next(iter(unique_high_exp))
            if v is not None:
                primary_hdu.header['HIGH_IMFE'] = float(v)
        if len(unique_low_mass) == 1:
            v = next(iter(unique_low_mass))
            if v is not None:
                primary_hdu.header['LOW_IMFM'] = float(v)
        if len(unique_high_mass) == 1:
            v = next(iter(unique_high_mass))
            if v is not None:
                primary_hdu.header['HIGH_IMFM'] = float(v)

        stellar_list = fits.HDUList([primary_hdu])
        for metallicity in sorted(met_groups.keys()):
            binaries = sorted(met_groups[metallicity], key=lambda b: b.age)

            # Use wavelength from first binary and apply wmin/wmax mask
            wavelength = binaries[0].wavelength
            wave_mask = np.ones(wavelength.size, dtype=bool)
            if wmin is not None:
                wave_mask &= wavelength >= wmin
            if wmax is not None:
                wave_mask &= wavelength <= wmax

            # Build spectra table (one row per wavelength, one column per age)
            table_data = {}
            for b in binaries:
                wave_arr = b.wavelength[wave_mask]
                flux_arr = b.flux[wave_mask]
                if do_rebin:
                    wave_arr, flux_arr, _ = spectrum_resampling(disp_intvl, pixel_width, pixel_number,
                                                                constant_pixel_width, wave_arr, flux_arr,
                                                                err_arr=None, mask_arr=None)
                table_data['WL'] = wave_arr
                table_data[str(b.age)] = flux_arr

            spec_table = Table(table_data)
            hdu = fits.BinTableHDU(spec_table)
            hdu.header['EXTNAME'] = ext_map[metallicity]
            hdu.header['METAL'] = metallicity
            stellar_list.append(hdu)

        stellar_list.writeto(output_path, overwrite=True)

        return

    @classmethod
    def from_fits(cls, fname, source=None):

        fname = Path(fname)
        if not fname.exists():
            raise SpecSyError(f'File not found: {fname}')

        binary_list = []

        with fits.open(fname) as hdul:

            # Read source and libraries from primary header
            primary_header = hdul[0].header
            file_source = primary_header.get('SOURCE', source)
            file_libs = primary_header.get('LIBS', None)

            # IMF metadata from primary header
            file_imf = primary_header.get('IMF', None) or None
            file_low_imfe = primary_header.get('LOW_IMFE', None)
            file_high_imfe = primary_header.get('HIGH_IMFE', None)
            file_low_imfm = primary_header.get('LOW_IMFM', None)
            file_high_imfm = primary_header.get('HIGH_IMFM', None)

            # Security check
            if file_source is None:
                raise SpecSyError(
                    f'No source information found in {fname}. Please provide it via the "source" argument.')

            # Get the metallicity map to reverse-lookup the extension names
            sources = [s.strip() for s in file_source.split(',')]
            if len(sources) == 1 and sources[0] in specsy_cfg['stellar']['ssp']:
                z_keys = specsy_cfg['stellar']['ssp'][sources[0]]['z_keys']
                ext_to_metal = {v_str: float(v) for v_str, v in z_keys.items()}
            else:
                ext_to_metal = None

            # Loop through extensions (skip primary)
            for hdu in hdul[1:]:

                # Get metallicity from header card first, then fall back to ext_map
                if 'METAL' in hdu.header:
                    metallicity = float(hdu.header['METAL'])
                elif ext_to_metal is not None:
                    if sources[0] == 'bpass':
                        metallicity = ext_to_metal[hdu.name.lower()]
                    elif sources[0] == 'starburst99':
                        metallicity = ext_to_metal[hdu.name]
                    elif sources[0] == 'pystarburst99':
                        metallicity = ext_to_metal[hdu.name]
                    else:
                        raise SpecSyError(
                            f'Cannot match the metallicity source "{sources[0]}" in "{hdu.name}" of file {fname}')
                else:
                    raise SpecSyError(f'Cannot determine metallicity for extension "{hdu.name}" in {fname}')

                # Read wavelength and flux columns
                table = Table(hdu.data)
                wave_arr = np.array(table['WL'])
                age_cols = [c for c in table.colnames if c != 'WL']

                # Determine library from header or filename
                library = file_libs.split(',')[0].strip() if file_libs else None

                for age_str in age_cols:
                    flux_arr = np.array(table[age_str])
                    binary_list.append(Binary(
                        source=file_source.strip(),
                        library=library,
                        metallicity=metallicity,
                        alpha=0.0,
                        imf=file_imf,
                        age=float(age_str),
                        fpath=str(fname),
                        flux=flux_arr,
                        wavelength=wave_arr,
                        low_imf_exp=float(file_low_imfe) if file_low_imfe is not None else None,
                        high_imf_exp=float(file_high_imfe) if file_high_imfe is not None else None,
                        low_imf_mass=float(file_low_imfm) if file_low_imfm is not None else None,
                        high_imf_mass=float(file_high_imfm) if file_high_imfm is not None else None,
                    ))

        return cls(binary_list)


    def to_dataframe(self, libraries=None, age_range=None, metallicity_range=None):

        # Filter binaries by requested libraries
        all_libraries = sorted({b.library for b in self._binaries.values()})
        if libraries is None:
            plot_libraries = all_libraries
        else:
            plot_libraries = [libraries] if isinstance(libraries, str) else list(libraries)
            unrecognised = set(plot_libraries) - set(all_libraries)
            if unrecognised:
                _logger.warning(f'Libraries not found in StellarBinaries: {unrecognised}. '
                                f'Available libraries: {all_libraries}')
            plot_libraries = [lib for lib in plot_libraries if lib in all_libraries]

        if len(plot_libraries) == 0:
            raise SpecSyError(f'No valid libraries to export. Available libraries: {all_libraries}')

        rows = []
        for b in self._binaries.values():
            if b.library not in plot_libraries:
                continue
            if age_range is not None and not (age_range[0] <= b.age <= age_range[1]):
                continue
            if metallicity_range is not None and not (metallicity_range[0] <= b.metallicity <= metallicity_range[1]):
                continue
            rows.append({'source': b.source,
                         'library': b.library,
                         'imf': b.imf,
                         'age': b.age,
                         'metallicity': b.metallicity,
                         'alpha': b.alpha,
                         'fpath': b.fpath,
                         'wmin': float(b.wavelength.min()) if b.is_loaded else np.nan,
                         'wmax': float(b.wavelength.max()) if b.is_loaded else np.nan})

        return pd.DataFrame(rows)


    def plot_age_metallicity(self, libraries=None, fig_cfg=None):

        # Filter binaries by requested libraries
        all_libraries = sorted({b.library for b in self._binaries.values()})
        if libraries is None:
            plot_libraries = all_libraries
        else:
            libraries = [libraries] if isinstance(libraries, str) else list(libraries)
            unrecognised = set(libraries) - set(all_libraries)
            if unrecognised:
                _logger.warning(f'Libraries not found in StellarBinaries: {unrecognised}. '
                                f'Available libraries: {all_libraries}')
            plot_libraries = [lib for lib in libraries if lib in all_libraries]

        if len(plot_libraries) == 0:
            raise SpecSyError(f'No valid libraries to plot. Available libraries: {all_libraries}')

        # Group by library
        library_groups = {}
        for b in self._binaries.values():
            if b.library in plot_libraries:
                library_groups.setdefault(b.library, ([], []))
                library_groups[b.library][0].append(b.age)
                library_groups[b.library][1].append(b.metallicity)


        use_colors = len(library_groups) > 1

        with rc_context(fig_cfg):

            fig, ax = plt.subplots(figsize=(8, 6))
            for i, (lib, (ages, metals)) in enumerate(library_groups.items()):
                ax.scatter(ages, metals, alpha=0.6, s=15,
                           color=f'C{i}' if use_colors else 'C0',
                           label=lib if use_colors else None)

            sources = sorted({b.source for b in self._binaries.values()})
            ax.set_title(f"Age-Metallicity coverage — source: {', '.join(sources)}")
            ax.set_xlabel("log(Age/yr)")
            ax.set_ylabel("Metallicity (Z)")
            ax.set_yscale('log')

            if use_colors:
                ax.legend(title='Library', markerscale=2)

            plt.tight_layout()
            plt.show()

        return

    def plot_spectra(self, fname=None, metallicity=None, age=None, metallicity_range=None, age_range=None,
                     libraries=None, log_scale=False, in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None, maximize=False):

        # Display check for input figures
        display_check = True if in_fig is _NO_FIG else False

        # Filter binaries by requested libraries
        all_libraries = sorted({b.library for b in self._binaries.values()})
        if libraries is None:
            plot_libraries = all_libraries
        else:
            plot_libraries = [libraries] if isinstance(libraries, str) else list(libraries)
            unrecognised = set(plot_libraries) - set(all_libraries)
            if unrecognised:
                _logger.warning(f'Libraries not found in StellarBinaries: {unrecognised}. '
                                f'Available libraries: {all_libraries}')
            plot_libraries = [lib for lib in plot_libraries if lib in all_libraries]

        if len(plot_libraries) == 0:
            raise SpecSyError(f'No valid libraries to plot. Available libraries: {all_libraries}')

        # Normalise metallicity and age selection to lists
        met_values = [metallicity] if metallicity is not None and not np.ndim(metallicity) else metallicity
        age_values = [age] if age is not None and not np.ndim(age) else age

        # Select binaries matching the filters
        selected = []
        for b in self._binaries.values():

            if b.library not in plot_libraries:
                continue
            if not b.is_loaded:
                continue

            # Exact metallicity match
            if met_values is not None:
                if not any(np.isclose(b.metallicity, m) for m in met_values):
                    continue

            # Metallicity range
            if metallicity_range is not None:
                if not (metallicity_range[0] <= b.metallicity <= metallicity_range[1]):
                    continue

            # Exact age match
            if age_values is not None:
                age_values = np.atleast_1d(age_values)
                if not any(np.isclose(b.age, a) for a in age_values):
                    continue

            # Age range
            if age_range is not None:
                age_range = np.atleast_1d(age_range)
                if not (age_range[0] <= b.age <= age_range[1]):
                    continue

            selected.append(b)

        if len(selected) == 0:
            raise SpecSyError('No loaded binaries match the requested filters.')

        # Adjust the default theme
        plt_cfg = {}
        if fig_cfg is not None:
            plt_cfg.update(fig_cfg)

        with rc_context(plt_cfg):

            self.fig = plt.figure() if (in_fig is None) or (in_fig is _NO_FIG) else in_fig
            self.ax = self.fig.add_subplot()

            # Default axis labels
            ax_labels_cfg = {
                'xlabel': 'Wavelength (Å)',
                'ylabel': 'Flux',
                'title': f'Stellar binary spectra — source: {self.source}',
            }
            if ax_cfg is not None:
                ax_labels_cfg.update(ax_cfg)

            # Color cycle across unique metallicities for visual grouping
            unique_mets = sorted({b.metallicity for b in selected})
            met_color = {m: f'C{i}' for i, m in enumerate(unique_mets)}

            legend_entries = {}
            for b in sorted(selected, key=lambda x: (x.metallicity, x.age)):
                color = met_color[b.metallicity]
                label_key = f'Z={b.metallicity:.5f} | {b.library}'
                self.ax.step(b.wavelength, b.flux, where='mid', alpha=0.7,
                             linewidth=0.8, color=color,
                             label=label_key if label_key not in legend_entries else None)
                legend_entries[label_key] = True

            self.ax.set(**ax_labels_cfg)

            if log_scale:
                self.ax.set_yscale('log')

            if len(legend_entries) <= 20:
                self.ax.legend(fontsize='small', title='Metallicity | Library')

            plt.tight_layout()
            save_close_fig_swicth(fname, 'tight', self.fig, maximize, display_check)

        return