from logging import getLogger
from pathlib import Path
from dataclasses import dataclass, field, replace

import re
import numpy as np
import pandas as pd
import xarray as xr

from matplotlib import pyplot as plt, rc_context
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.table import Table

from lime import save_frame, load_frame
from lime.workflow import spectrum_resampling
from lime.plotting.plots import _NO_FIG, save_close_fig_swicth
from specsy.io import specsy_cfg, SpecSyError
from specsy.models import tools_pySB99 as pySB99
from lime.tools import au
from astropy import constants as const

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


def bpass_file_manager(fname, load_spectrum, metal_map=specsy_cfg['stellar']['ssp']['bpass']['z_keys']):

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


# ── VMS_BPASS (Martins+25) ───────────────────────────────────────────────────
# Two filename flavours coexist in the same folder:
#   pop_vms_bpass_ZKEY_AGEmyr_MXXX_(gr|grzsc)[.neb].dat  -> age < 2.5 Myr, VMS up to 300 Msun
#   bpass_(sing|bin)_imfXXX_YYY.ZKEY.AGEmyr[.neb].dat    -> older ages, no VMS
# "gr"/"grzsc": mass loss rates are not/are scaled with metallicity.

_VMS_BPASS_DEFAULTS = {'z_keys': {'zsmc':  0.002, 'z0p1':  0.001, 'z0p01': 0.0001,      # pop_vms_ files
                                  'z002':  0.002, 'z001':  0.001, 'z1em4': 0.0001}}     # bpass_sing_ files

VMS_BPASS_CFG = specsy_cfg['stellar']['ssp'].get('VMS_BPASS', _VMS_BPASS_DEFAULTS)

VMS_ZERO_AGE_MYR = 0.0001      # The "0myr" models have no logarithm: they are floored at this age

_VMS_PATTERN = re.compile(r'^pop_vms_bpass_(z[\w]+?)_([\dp]+)myr_M(\d+)_(gr|grzsc)(\.neb)?$')

_SING_PATTERN = re.compile(r'^bpass_(sing|bin)_(imf\d+_\d+)\.(z[\w]+)\.([\dp]+)myr(\.neb)?$')


def _myr_string_to_log_age(age_str, zero_age_myr=VMS_ZERO_AGE_MYR):

    """Convert the filename age token into log(age/yr): '0p5' -> 5.70, '12p5' -> 7.10."""

    age_myr = float(age_str.replace('p', '.'))
    zero_age = age_myr == 0
    age_myr = zero_age_myr if zero_age else age_myr

    return np.log10(age_myr * 1e6), zero_age


def vms_bpass_file_manager(fname, load_spectrum, met_map=None, sed_type='stellar_and_nebular',
                           vms_imf_exp=1.35, zero_age_myr=VMS_ZERO_AGE_MYR):

    """
    Parse one VMS_BPASS file name into the ``Binary`` metadata plus its spectrum.

    The mass loss prescription (gr/grzsc), the population type (sing/bin) and the presence of the
    nebular continuum are encoded in the ``library`` entry ('VMS_grzsc', 'BPASS_sing',
    'VMS_gr_neb', ...), so files sharing age and metallicity do not overwrite each other and can be
    selected later with the ``libraries`` argument of ``to_stellar_binaries``, ``to_fits`` or ``plot_spectra``.

    ``sed_type`` selects the files to keep: 'stellar' (no .neb), 'nebular' (only .neb) or
    'stellar_and_nebular' (both). Returns ``None`` for a discarded file.
    """

    met_map = VMS_BPASS_CFG['z_keys'] if met_map is None else met_map
    fname = Path(fname)
    name = fname.stem

    m_vms, m_sing = _VMS_PATTERN.match(name), _SING_PATTERN.match(name)

    # Very massive star populations
    if m_vms is not None:
        z_str, age_str, high_mass, mass_loss, neb_str = m_vms.groups()
        library = f'VMS_{mass_loss}'
        imf_str = f'imf{int(round(vms_imf_exp * 100))}_{high_mass}'
        low_imf_exp, high_imf_mass = vms_imf_exp, float(high_mass)

    # Standard BPASS populations for the older ages
    elif m_sing is not None:
        pop_str, imf_str, z_str, age_str, neb_str = m_sing.groups()
        library = f'BPASS_{pop_str}'
        imf_match = re.match(r'imf(\d+)_(\d+)', imf_str)
        low_imf_exp, high_imf_mass = float(imf_match.group(1)) / 100, float(imf_match.group(2))

    else:
        raise ValueError(f'Unrecognised filename format for VMS_BPASS files: {fname}')

    # Nebular files get their own library entry so they do not collide with the stellar ones
    nebular = neb_str is not None
    if (sed_type == 'stellar') and nebular:
        return None
    if (sed_type == 'nebular') and not nebular:
        return None
    library = f'{library}_neb' if nebular else library

    # Metallicity and age
    if z_str not in met_map:
        raise SpecSyError(f'Metallicity key "{z_str}" from {fname.name} is not in the VMS_BPASS map: '
                          f'{list(met_map)}')
    metallicity = met_map[z_str]
    log_age, zero_age = _myr_string_to_log_age(age_str, zero_age_myr)

    params = dict(library=library, metallicity=metallicity, alpha=0.0, imf=imf_str,
                  low_imf_exp=low_imf_exp, high_imf_exp=None,
                  low_imf_mass=None, high_imf_mass=high_imf_mass)

    # One age per file: a single flux column
    wave_arr, flux_arr = None, None
    if load_spectrum:
        data_arr = np.loadtxt(fname)
        wave_arr, flux_arr = data_arr[:, 0], data_arr[:, 1]
        if data_arr.shape[1] > 2:
            _logger.warning(f'File {fname.name} has {data_arr.shape[1]} columns: only the first two '
                            f'(wavelength, flux) are read for the VMS_BPASS source.')

    return params, log_age, zero_age, wave_arr, flux_arr


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

        case 'VMS_BPASS':

            binary_list, zero_age_count = [], 0
            for fname in file_list:

                parsed = vms_bpass_file_manager(fname, load_spectrum, sed_type=sed_type)

                # File discarded by the sed_type selection
                if parsed is None:
                    continue

                bin_params, log_age, zero_age, wave_arr, flux_arr = parsed
                zero_age_count += int(zero_age)

                binary_list.append(Binary(source='vms_bpass',
                                          age=log_age,
                                          fpath=str(fname),
                                          wavelength=wave_arr,
                                          flux=flux_arr,
                                          **bin_params))

            if zero_age_count > 0:
                _logger.warning(f'{zero_age_count} VMS_BPASS files declare a 0 Myr age, which has no logarithm. '
                                f'They were assigned {VMS_ZERO_AGE_MYR} Myr '
                                f'(log(age) = {np.log10(VMS_ZERO_AGE_MYR * 1e6):.2f}).')

            if len(binary_list) == 0:
                raise SpecSyError(f'No VMS_BPASS files left after the sed_type="{sed_type}" selection of '
                                  f'{len(file_list)} files.')

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


def to_dataframe(binary_list, columns=None):

    """
    Table with one row per binary and one column per property. By default, the columns are the
    specsy_cfg['stellar']['ssp_params'] entries.
    """

    columns = specsy_cfg['stellar']['ssp_params'] if columns is None else list(columns)

    return pd.DataFrame([{col: getattr(b, col) for col in columns} for b in binary_list], columns=columns)


def _to_python_scalar(value):
    """Frame cell to Binary field: NaN -> None and numpy scalars -> python scalars."""
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


_NOT_LOADED_MSG = 'Please create the object with StellarBinaries.from_source(..., load_spectra=True).'


# Physical constants
L_SUN_CGS = const.L_sun.cgs.value                                       # erg s^-1
HC_CGS = (const.h * const.c).cgs.value                                  # erg cm
HC_EV_ANGSTROM = (const.h * const.c).to(au.eV * au.AA).value            # eV Å (lambda = hc / E)
Q_LIMITS_LAMBDA = {'HI': 911.753, 'HeI': 504.259, 'HeII': 227.838}      # Å


# Energy (eV) to create the ion of each MIRI line, i.e. the ionization potential of the previous stage
MIRI_ION_LINES = {r'[Ne II] 12.81 $\mu$m': 21.56,
                 r'[S IV] 10.51 $\mu$m': 34.86,
                 r'[Ne III] 15.56 $\mu$m': 40.96,
                 r'[O IV] 25.89 $\mu$m': 54.9,
                 r'[Ar V] 7.90|13.10 $\mu$m': 59.8,
                 r'[Ne V] 14.32|24.32 $\mu$m': 97.1,
                 r'[Ne VI] 7.65 $\mu$m': 126.2}


def log_ionizing_photons(wave_arr, flux_arr, units_wave, units_flux, lambda_limit=Q_LIMITS_LAMBDA['HI']):

    """
    Logarithm of the photon rate blueward of a wavelength limit.

    Computes the number of photons per second emitted below `lambda_limit` by a spectrum given as a luminosity
    per unit wavelength:

    .. math::  Q = \\int_0^{\\lambda_{lim}} L_\\lambda \\, \\frac{\\lambda}{hc} \\, d\\lambda

    The integral is evaluated with the trapezoidal rule and stops exactly at `lambda_limit`, using the flux
    linearly interpolated at that wavelength. The units are not converted: the spectrum must be in Angstroms and
    solar luminosities per Angstrom (the BPASS convention).

    Parameters
    ----------
    wave_arr : array_like
        Wavelength array in ascending order, in Angstroms.
    flux_arr : array_like
        Luminosity per unit wavelength in Lsun / Angstrom. It must have the same length as `wave_arr`.
    units_wave : str or astropy.units.Unit
        Units of `wave_arr`. It must be Angstrom.
    units_flux : str or astropy.units.Unit
        Units of `flux_arr`. It must be Lsun / Angstrom.
    lambda_limit : float, optional
        Wavelength limit in Angstroms. The default is the hydrogen ionization edge (911.753 Å).

    Returns
    -------
    float
        log10 of the photon rate Q in s^-1, or -inf if the spectrum has no flux below `lambda_limit`.

    Raises
    ------
    SpecSyError
        If a unit is not a valid astropy unit or is not Angstrom (wavelength) or Lsun / Angstrom (flux), or if
        `wave_arr` starts redward of `lambda_limit`.

    Notes
    -----
    The flux below the first wavelength of the grid is taken as zero, so the result is a lower limit if the
    spectrum does not reach the ionizing regime of interest. NaN values in `flux_arr` give -inf.

    Examples
    --------
    >>> log_q = log_ionizing_photons(wave, flux, 'Angstrom', 'Lsun / Angstrom')
    >>> log_q_he = log_ionizing_photons(wave, flux, 'Angstrom', 'Lsun / Angstrom', lambda_limit=Q_LIMITS_LAMBDA['HeII'])
    """

    # The units are not converted: they must be Å and Lsun/Å
    for name, unit, expected in (('units_wave', units_wave, au.AA), ('units_flux', units_flux, au.Lsun / au.AA)):
        try:
            unit = au.Unit(unit)
        except (ValueError, TypeError) as e:
            raise SpecSyError(f'The {name} "{unit}" is not a valid astropy unit: {e}')
        if unit != expected:
            raise SpecSyError(f'The {name} "{unit}" is not "{expected}". The units are not converted: '
                              f'please provide the spectrum in {expected}.')

    if wave_arr[0] >= lambda_limit:
        raise SpecSyError(f'The spectrum starts at {wave_arr[0]:.1f} Å, redward of the {lambda_limit:.1f} Å edge. '
                          f'Please compute Q before cropping the wavelength range.')

    # Integrate up to the edge exactly
    mask = wave_arr < lambda_limit
    wave_int = np.append(wave_arr[mask], lambda_limit)
    flux_int = np.append(flux_arr[mask], np.interp(lambda_limit, wave_arr, flux_arr))

    # Number of ionization photons: Lsun -> erg/s, and 1e-8 takes the lambda in hc/lambda from Å to cm
    q_ion = np.trapezoid(flux_int * wave_int, wave_int) * au.Lsun.to(au.erg / au.s) * au.AA.to(au.cm) / HC_CGS

    return np.log10(q_ion) if q_ion > 0 else -np.inf


@dataclass(eq=False)
class SSPGrid:

    """
    Everything needed to evaluate the nearest-node SSP model: observation, extinction transmission, the (visited or
    full) grid fluxes with their (log age, log10 Z) node coordinates, the full age/Z axes and the prior bounds.

    It is built by SSP_sampler.prepare_inputs (full grid), reduced to the posterior-visited nodes with `subset` for
    storage, and rebuilt from a saved trace with `from_idata`. The same lookup/model code is therefore used while
    sampling and when post-processing.

    Fields
    ------
    grid_flux : (n_node, n_pixel) total model flux per node (stellar + nebular if the latter was fitted), Lsun/A.
    grid_neb  : (n_node, n_pixel) nebular part of grid_flux, or None if it was not fitted.
    age_node, z_node : (n_node,) log(age) and log10(Z) of each row of grid_flux.
    age_log, z_log   : full ascending axes (unique node values); node_idx[i_z, i_age] gives the row (-1 = missing).
    bounds : (4, 2) prior limits in PARAM_NAMES order.

    Dict-style access (grid['wave']) is kept so code written for the former dict returned by ssp_inputs still works.
    """

    wave: np.ndarray
    lum: np.ndarray
    lum_err: np.ndarray
    mask: np.ndarray
    red_corr: np.ndarray
    grid_flux: np.ndarray
    age_node: np.ndarray
    z_node: np.ndarray
    age_log: np.ndarray
    z_log: np.ndarray
    bounds: np.ndarray
    grid_neb: np.ndarray = None
    red_law: str = 'CCM'
    r_v: float = 3.1
    node_idx: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):

        for name in ('wave', 'lum', 'lum_err', 'red_corr', 'grid_flux', 'age_node', 'z_node', 'age_log', 'z_log',
                     'bounds'):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        self.mask = np.asarray(self.mask, dtype=bool)
        if self.grid_neb is not None:
            self.grid_neb = np.asarray(self.grid_neb, dtype=float)

        # Flux row of each (Z, age) cell, -1 for combinations missing in the grid (nodes are taken from the axes,
        # so the matches are exact)
        i_z = np.searchsorted(self.z_log, self.z_node)
        i_age = np.searchsorted(self.age_log, self.age_node)
        self.node_idx = np.full((self.z_log.size, self.age_log.size), -1, dtype=int)
        self.node_idx[i_z, i_age] = np.arange(self.age_node.size)

        return

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def add_nebular(self):
        return self.grid_neb is not None

    def rows(self, log_age, log_met):
        """Flux row of the grid node nearest to each (log_age, log_met), -1 if the combination is missing."""
        return self.node_idx[self._nearest(self.z_log, log_met), self._nearest(self.age_log, log_age)]

    def model(self, theta):
        """Rescaled, reddened model of the node nearest to theta (None if the node is missing from the grid)."""
        log_age, log_met, ebv, log_amp = theta
        row = int(self.rows(log_age, log_met)[0])
        return None if row < 0 else np.power(10.0, log_amp) * self.grid_flux[row] * np.power(self.red_corr, ebv)

    def components(self, theta):
        """Stellar and nebular parts of the rescaled, reddened model (nebular is None if it was not fitted)."""
        log_age, log_met, ebv, log_amp = theta
        row = int(self.rows(log_age, log_met)[0])
        if row < 0:
            return None, None

        scale = np.power(10.0, log_amp) * np.power(self.red_corr, ebv)
        total = scale * self.grid_flux[row]
        if self.grid_neb is None:
            return total, None

        neb = scale * self.grid_neb[row]
        return total - neb, neb

    def subset(self, rows):
        """Copy keeping only the given flux rows (e.g. the nodes visited by the posterior); axes are unchanged."""
        return replace(self, grid_flux=self.grid_flux[rows], age_node=self.age_node[rows], z_node=self.z_node[rows],
                       grid_neb=None if self.grid_neb is None else self.grid_neb[rows])

    @classmethod
    def from_idata(cls, idata):
        """Rebuilds the grid from the constant_data group of a trace."""
        cd = idata.constant_data
        return cls(wave=cd['wave'].values, lum=cd['lum'].values, lum_err=cd['lum_err'].values,
                   mask=cd['mask'].values, red_corr=cd['red_corr'].values, grid_flux=cd['grid_flux'].values,
                   age_node=cd['age_node'].values, z_node=cd['z_node'].values,
                   age_log=cd['age_axis'].values, z_log=cd['z_axis'].values, bounds=cd['bounds'].values,
                   grid_neb=cd['grid_neb'].values if 'grid_neb' in cd else None,
                   red_law=cd.attrs['red_law'], r_v=cd.attrs['r_v'])

    # inside SSPGrid
    @staticmethod
    def _nearest(axis, x):
        """Vectorized nearest-node lookup: index of the axis value nearest to each x (ties go to the upper node)."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        if axis.size == 1:
            return np.zeros(x.size, dtype=int)
        i = np.clip(np.searchsorted(axis, x, side='left'), 1, axis.size - 1)
        return i - ((x - axis[i - 1]) < (axis[i] - x)).astype(int)

    def rows(self, log_age, log_met):
        """Flux row of the grid node nearest to each (log_age, log_met), -1 if the combination is missing."""
        return self.node_idx[self._nearest(self.z_log, log_met), self._nearest(self.age_log, log_age)]


class StellarBinaries:


    def __init__(self, binary_list):

        binary_list = list(binary_list)
        if len(binary_list) == 0:
            raise SpecSyError('Cannot create a StellarBinaries object from an empty binary list.')

        # Properties table and spectra series (the Binary objects are not stored)
        self._set_data(frame=to_dataframe(binary_list),
                       wave_list=[b.wavelength if b.is_loaded else None for b in binary_list],
                       flux_list=[b.flux for b in binary_list])

        return


    def __repr__(self):

        # Sources and libraries from attributes
        sources = [self.source] if isinstance(self.source, str) else self.source
        libraries = [self.library] if isinstance(self.library, str) else self.library
        sources = sources if sources is not None else []
        libraries = libraries if libraries is not None else []

        n_loaded = int(self._loaded_mask().sum())

        # Wavelength uniformity # TODO add units attribute
        wmin_str = f'{self.uniform_wmin:.1f} Å' if self.uniform_wmin is not None else 'non-uniform'
        wmax_str = f'{self.uniform_wmax:.1f} Å' if self.uniform_wmax is not None else 'non-uniform'
        delta_str = f'{self.uniform_deltalambda:.3f} Å' if self.uniform_deltalambda is not None else 'non-uniform'

        # Ages and metallicities
        age_str = ', '.join(f'{a:.1f}' for a in self.ages)
        met_str = ', '.join(f'{m:.5f}' for m in self.metallicities)

        return (f'\n{"=" * 60}\n'
                f'  StellarBinaries\n'
                f'{"=" * 60}\n'
                f'  Source(s)      : {", ".join(map(str, sources))}\n'
                f'  Libraries      : {", ".join(map(str, libraries))}\n'
                f'  Total stellar spectra : {len(self.frame)} ({n_loaded} loaded)\n'
                f'{"-" * 60}\n'
                f'  Wavelength grid\n'
                f'    wmin         : {wmin_str}\n'
                f'    wmax         : {wmax_str}\n'
                f'    delta lambda : {delta_str}\n'
                f'{"-" * 60}\n'
                f'  Ages (log yr)  [{len(self.ages)}]: {age_str}\n'
                f'{"-" * 60}\n'
                f'  Metallicities  [{len(self.metallicities)}]: {met_str}\n'
                f'{"=" * 60}\n')


    @property
    def is_loaded(self):
        return bool(self._loaded_mask().all())


    @classmethod
    def from_source(cls, source_name, source_folder, root_folder=None, load_spectra=False, ext=None, file_list=None,
                    sed_type='stellar_and_nebular'):

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

        # pySB99: source_folder is a root containing one subfolder per run
        if source_name == 'pystarburst99':
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

        # BPASS and VMS_BPASS: file-based logic
        if file_list is not None:
            list_fname = [source_folder / f for f in file_list]
            missing = [f for f in list_fname if not f.exists()]
            if missing:
                raise SpecSyError(f'Files not found in {source_folder}: {missing}')
        else:
            list_fname = list(source_folder.glob(f'*{"" if ext is None else ext}'))
            if len(list_fname) == 0:
                raise SpecSyError(f'No files found in {source_folder} with extension {ext}')

        binary_list = parse_binaries(source_name, list_fname, load_spectra, sed_type)
        return cls(binary_list)

    @classmethod
    def from_frame(cls, frame, wave_list=None, flux_list=None, quiet=False, **kwargs):

        """
        New StellarBinaries from a properties table. `frame` is either:
          - a pandas DataFrame already selected in memory, paired with `wave_list`/`flux_list` holding the
            matching spectra in the same row order (internal use, e.g. to_stellar_binaries).
          - a file path, read with lime.load_frame(frame, **kwargs). A saved frame has no spectra, so the
            resulting object is unloaded; use from_source(..., load_spectra=True) to get the spectra too.

        `quiet=True` silences the column/index warnings below — used internally, where the frame is
        already known to follow the expected structure.
        """

        # Read the frame from disk if a path (or anything else lime.load_frame accepts) was given
        if not isinstance(frame, pd.DataFrame):
            frame = load_frame(frame, **kwargs)
            # wave_list = [None] * len(frame)
            # flux_list = [None] * len(frame)

        # The frame must have every expected column; anything else is dropped
        expected = specsy_cfg['stellar']['ssp_params']
        missing = [col for col in expected if col not in frame.columns]
        if missing:
            raise SpecSyError(f'The frame is missing the columns: {missing}. Expected: {expected}')

        unknown = [col for col in frame.columns if col not in expected]
        if unknown:
            if not quiet:
                _logger.warning(f'The frame has columns not recognised by StellarBinaries: {unknown}. '
                                f'They will be ignored.')
            frame = frame[expected]

        # The frame index must line up positionally with wave_list/flux_list; _set_data resets it regardless
        if not quiet and not pd.api.types.is_numeric_dtype(frame.index):
            _logger.warning('The frame index is not numeric. It will be reset to a default integer index.')

        new_obj = cls.__new__(cls)
        new_obj._set_data(frame, wave_list, flux_list)

        return new_obj


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

            # Alpha and IMF metadata from primary header
            file_alpha = primary_header.get('ALPHA', 0.0)
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
            if len(sources) == 1 and sources[0] == 'vms_bpass':
                ext_to_metal = {_metallicity_to_key(v): float(v) for v in set(VMS_BPASS_CFG['z_keys'].values())}
            elif len(sources) == 1 and sources[0] in specsy_cfg['stellar']['ssp']:
                z_keys = specsy_cfg['stellar']['ssp'][sources[0]]['z_keys']
                ext_to_metal = {v_str: float(v) for v_str, v in z_keys.items()}
            else:
                ext_to_metal = None

            # Library for files without a LIBRARY card in the extensions
            default_library = file_libs.split(',')[0].strip() if file_libs else None

            # Loop through extensions (skip primary)
            for hdu in hdul[1:]:

                # Get metallicity from header card first, then fall back to ext_map
                if 'METAL' in hdu.header:
                    metallicity = float(hdu.header['METAL'])
                elif ext_to_metal is not None:
                    if sources[0] == 'bpass':
                        metallicity = ext_to_metal[hdu.name.lower()]
                    elif sources[0] in ('pystarburst99', 'vms_bpass'):
                        metallicity = ext_to_metal[hdu.name]
                    else:
                        raise SpecSyError(
                            f'Cannot match the metallicity source "{sources[0]}" in "{hdu.name}" of file {fname}')
                else:
                    raise SpecSyError(f'Cannot determine metallicity for extension "{hdu.name}" in {fname}')

                # Library from the extension header, then from the primary header
                library = hdu.header.get('LIBRARY', default_library)

                # Read wavelength and flux columns
                table = Table(hdu.data)
                wave_arr = np.array(table['WL'])
                age_cols = [c for c in table.colnames if c != 'WL']

                for age_str in age_cols:
                    flux_arr = np.array(table[age_str])
                    binary_list.append(Binary(
                        source=file_source.strip(),
                        library=library,
                        metallicity=metallicity,
                        alpha=float(file_alpha),
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


    @staticmethod
    def _write_fits(output_path, selection, wave_list, flux_list, ext_map):

        # Frame labels matching the positions in the spectra lists
        selection = selection.reset_index(drop=True)
        sel_sources = sorted(selection['source'].unique(), key=str)
        sel_source = sel_sources[0] if len(sel_sources) == 1 else sel_sources
        sel_mets = np.unique(selection['metallicity'])

        # Dictionary for mapping the fits header metallicity values according to the dataset
        if ext_map is None:
            if sel_source == 'bpass':
                BPASS_Z_KEYS = specsy_cfg['stellar']['ssp']['bpass']['z_keys']
                ext_map = {value: key for key, value in BPASS_Z_KEYS.items()}
            elif sel_source == 'pystarburst99':
                pySB99_Z_KEYS = specsy_cfg['stellar']['ssp']['pystarburst99']['z_keys']
                ext_map = {met_val: _metallicity_to_key(met_val) for name, met_val in pySB99_Z_KEYS.items()}
            elif sel_source == 'vms_bpass':
                ext_map = {met_val: _metallicity_to_key(met_val) for met_val in sel_mets}
                if len(set(ext_map.values())) != len(ext_map):
                    raise SpecSyError(f'The VMS_BPASS metallicities {list(sel_mets)} do not produce unique '
                                      f'extension names ({sorted(ext_map.values())}). Please provide an "ext_map".')
            else:
                raise SpecSyError(f'Binaries source "{sel_source}" does not have a metallicity header map')

        missing_mets = [met for met in sel_mets if met not in ext_map]
        if missing_mets:
            raise SpecSyError(f'The metallicities {missing_mets} are not in the extension map. '
                              f'Please provide an "ext_map".')

        # Build FITS primary header
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['SOURCE'] = ', '.join(map(str, sel_sources))
        primary_hdu.header['LIBS'] = ', '.join(sorted(selection['library'].dropna().unique()))

        # Coverage actually written to the file
        primary_hdu.header['MET_MIN'] = float(selection['metallicity'].min())
        primary_hdu.header['MET_MAX'] = float(selection['metallicity'].max())
        primary_hdu.header['AGE_MIN'] = float(selection['age'].min())
        primary_hdu.header['AGE_MAX'] = float(selection['age'].max())

        # Alpha and IMF cards if uniform across the selection
        for card, col in (('ALPHA', 'alpha'), ('IMF', 'imf'), ('LOW_IMFE', 'low_imf_exp'),
                          ('HIGH_IMFE', 'high_imf_exp'), ('LOW_IMFM', 'low_imf_mass'), ('HIGH_IMFM', 'high_imf_mass')):
            values = selection[col].unique()
            if (len(values) == 1) and not pd.isna(values[0]):
                primary_hdu.header[card] = values[0] if isinstance(values[0], str) else float(values[0])

        # Metallicities shared by several libraries get the library in the extension name
        libs_per_met = selection.groupby('metallicity')['library'].nunique(dropna=False)

        stellar_list = fits.HDUList([primary_hdu])
        for (metallicity, library), group in selection.groupby(['metallicity', 'library'], sort=True, dropna=False):

            library = None if pd.isna(library) else library

            # Spectra with the same age would overwrite each other's column
            if group['age'].duplicated().any():
                varying = [col for col in group.columns
                           if (col not in ('age', 'fpath')) and (group[col].nunique(dropna=False) > 1)]
                raise SpecSyError(f'Several spectra share the same age for metallicity={metallicity}, '
                                  f'library={library}. They differ in: {varying if varying else "none"}. '
                                  f'Please export them in separate files.')

            group = group.sort_values('age')

            # Build spectra table (one row per wavelength, one column per age)
            wave_ref = wave_list[group.index[0]]
            table_data = {'WL': wave_ref}
            for idx, age in zip(group.index, group['age']):
                wave_arr = wave_list[idx]
                if (wave_arr is not wave_ref) and not np.array_equal(wave_arr, wave_ref):
                    raise SpecSyError(f'The spectra with metallicity={metallicity}, library={library} do not share '
                                      f'the same wavelength array. Please resample them with the "disp_intvl", '
                                      f'"pixel_width" or "pixel_number" arguments.')
                table_data[str(float(age))] = flux_list[idx]

            extname = ext_map[metallicity]
            if libs_per_met.loc[metallicity] > 1:
                extname = f'{extname}_{library}'

            hdu = fits.BinTableHDU(Table(table_data))
            hdu.header['EXTNAME'] = extname
            hdu.header['METAL'] = float(metallicity)
            if library is not None:
                hdu.header['LIBRARY'] = library
            stellar_list.append(hdu)

        stellar_list.writeto(output_path, overwrite=True)

        return


    @staticmethod
    def _crop_spectrum(wave_arr, flux_arr, wmin, wmax):

        if (wmin is None) and (wmax is None):
            return wave_arr, flux_arr

        wave_mask = np.ones(wave_arr.size, dtype=bool)
        if wmin is not None:
            wave_mask &= wave_arr >= wmin
        if wmax is not None:
            wave_mask &= wave_arr <= wmax

        return wave_arr[wave_mask], flux_arr[wave_mask]


    def get_binaries(self, **params):

        """
        List of Binary objects matching the requested frame values:
            sb.get_binaries(library='AP')
            sb.get_binaries(age=[7.0, 8.0], metallicity=0.014)
        """

        binary_list = []
        for idx in self._select(**params):
            row = {col: _to_python_scalar(value) for col, value in self.frame.loc[idx].items()}
            wave_arr, flux_arr = self._spectrum(idx)
            binary_list.append(Binary(**row, wavelength=wave_arr, flux=flux_arr))

        return binary_list


    def get_spectrum(self, age, metallicity, interpolate=False, plot_results=False, **params):

        """
        Wavelength and flux arrays of the stellar spectrum with the requested age (log yr) and metallicity (Z).

        With interpolate=False the values must exist in the frame. With interpolate=True, the spectrum is linearly
        interpolated from the four neighbouring grid nodes, in log(age) and log10(Z). A requested value on a grid
        node returns that node's spectrum. Values outside the grid coverage raise an error (no extrapolation).
        The extra keyword arguments (library, alpha, imf...) select the rows first, as in get_binaries().

        With plot_results=True, an interpolated spectrum is plotted alongside the grid node spectra used to
        compute it (no plot if the request falls on a grid node).
        """

        # Exact match
        if not interpolate:
            idcs = self._select(age=age, metallicity=metallicity, **params)
            return self._single_spectrum(idcs, age, metallicity)

        # Grid nodes available for the requested library, alpha, imf...
        subset = self.frame.loc[self._select(**params)]
        grid_ages = np.unique(subset['age'])
        grid_mets = np.unique(subset['metallicity'])

        # Security checks
        for name, value, grid in (('age', age, grid_ages), ('metallicity', metallicity, grid_mets)):
            if not (grid[0] <= value <= grid[-1]):
                raise SpecSyError(f'The requested {name}={value} is outside the grid coverage ({grid[0]} to '
                                  f'{grid[-1]}). The spectra are not extrapolated.')

        # Neighbouring nodes and weights: age is already log(yr), metallicity in log10(Z)
        age_i0, age_i1, age_w = self._bracket(grid_ages, age)
        met_i0, met_i1, met_w = self._bracket(grid_mets, metallicity, log=True)

        weights = {}
        for age_idx, a_w in ((age_i0, 1.0 - age_w), (age_i1, age_w)):
            for met_idx, z_w in ((met_i0, 1.0 - met_w), (met_i1, met_w)):
                if a_w * z_w > 0:
                    node = (grid_ages[age_idx], grid_mets[met_idx])
                    weights[node] = weights.get(node, 0.0) + a_w * z_w

        # Weighted sum of the node spectra
        wave_arr, flux_arr, nodes = None, 0.0, []
        for (node_age, node_met), weight in weights.items():

            node_mask = ((subset['age'] == node_age) & (subset['metallicity'] == node_met)).to_numpy()
            if not node_mask.any():
                raise SpecSyError(f'The grid node age={node_age}, metallicity={node_met} needed to interpolate '
                                  f'age={age}, metallicity={metallicity} is missing.')

            node_wave, node_flux = self._single_spectrum(subset.index[node_mask], node_age, node_met)
            if wave_arr is None:
                wave_arr = node_wave
            elif (node_wave is not wave_arr) and not np.array_equal(node_wave, wave_arr):
                raise SpecSyError(f'The spectra around age={age}, metallicity={metallicity} do not share the same '
                                  f'wavelength array. Please resample them with to_stellar_binaries(...).')

            flux_arr = flux_arr + weight * node_flux
            nodes.append((node_age, node_met, weight, node_flux))

        # Interpolated spectrum against the grid nodes used (only if an interpolation took place)
        if plot_results:

            fig, ax = plt.subplots(figsize=(10, 6))
            for i, (node_age, node_met, weight, node_flux) in enumerate(nodes):
                ax.step(wave_arr, node_flux, where='mid', linewidth=0.8, alpha=0.6, color=f'C{i}',
                        label=f'log(age)={node_age:.2f}, Z={node_met:.5f} (w={weight:.3f})', linestyle='--')

            ax.step(wave_arr, flux_arr, where='mid', linewidth=1.2, color='black',
                    label=f'Interpolated: log(age)={age:.3f}, Z={metallicity:.5f}')

            ax.set(xlabel='Wavelength (Å)', ylabel='Flux',
                   title=f'Stellar spectrum interpolation — source: {self.source}')
            ax.legend(fontsize='small')
            plt.tight_layout()
            plt.show()

        return wave_arr, flux_arr


    def to_stellar_binaries(self, libraries=None, wmin=None, wmax=None, met_limits=None, age_limits=None,
                            disp_intvl=None, pixel_width=None, pixel_number=None, constant_pixel_width=True,
                            idcs_ssp=None):

        """
        New StellarBinaries with the selected spectra, cropped to wmin-wmax and resampled if one of `disp_intvl`,
        `pixel_width` or `pixel_number` is provided (common output dispersion for all the spectra). The rows are
        selected with `idcs_ssp` (frame indices or boolean mask), which takes precedence over the `libraries`,
        `met_limits` and `age_limits` arguments. The new frame index is reset.
        """

        selection, wave_list, flux_list = self._process_selection(libraries, wmin, wmax, met_limits, age_limits,
                                                                  disp_intvl, pixel_width, pixel_number,
                                                                  constant_pixel_width, idcs_ssp)

        return self.from_frame(selection, wave_list, flux_list)


    def to_fits(self, fname, libraries=None, wmin=None, wmax=None, ext_map=None, met_limits=None, age_limits=None,
                disp_intvl=None, pixel_width=None, pixel_number=None, constant_pixel_width=True, idcs_ssp=None):

        """
        Export the selected spectra to a FITS file with one table extension per metallicity (one row per wavelength
        and one column per age). If several libraries share a metallicity, each one gets its own extension named
        '{metallicity key}_{library}'. The selection, crop and resampling arguments behave as in to_stellar_binaries.
        """

        # Check output folder exists
        output_path = Path(fname)
        if not output_path.parent.exists():
            raise SpecSyError(f'Output folder not found: {output_path.parent}')

        selection, wave_list, flux_list = self._process_selection(libraries, wmin, wmax, met_limits, age_limits,
                                                                  disp_intvl, pixel_width, pixel_number,
                                                                  constant_pixel_width, idcs_ssp)

        self._write_fits(output_path, selection, wave_list, flux_list, ext_map)

        return


    def ionizing_photons(self, edge='HI', flux_to_cgs=L_SUN_CGS, libraries=None, met_limits=None, age_limits=None, idcs_ssp=None):

        """
        log10 Q (photons s^-1) blueward of `edge` ('HI', 'HeI', 'HeII' or a wavelength in Å) for the selected
        spectra, as a Series with the frame index. The selection arguments behave as in to_stellar_binaries.
        Use the stellar spectra on their full wavelength range (before any crop).
        """

        edge_wave = Q_LIMITS_LAMBDA[edge] if isinstance(edge, str) else float(edge)
        idcs = self._selection_index(idcs_ssp, libraries, met_limits, age_limits)

        n_missing = int((~self._loaded_mask()[idcs.to_numpy()]).sum())
        if n_missing > 0:
            raise SpecSyError(f'{n_missing} of the {len(idcs)} selected spectra are not loaded. {_NOT_LOADED_MSG}')

        # The nebular files do not carry the bare stellar ionizing output
        neb_libs = [lib for lib in self.frame.loc[idcs, 'library'].unique() if str(lib).endswith('_neb')]
        if neb_libs:
            _logger.warning(f'The selection includes the nebular libraries {neb_libs}. Their flux below '
                            f'{edge_wave:.1f} Å may not be the stellar ionizing output.')

        log_q = [log_ionizing_photons(*self._spectrum(idx), edge=edge_wave, flux_to_cgs=flux_to_cgs) for idx in idcs]

        return pd.Series(log_q, index=idcs, name=f'logQ_{edge}')


    def to_ionization_table(self, fname, edge='HI', flux_to_cgs=L_SUN_CGS, ext_map=None, log_q_floor=0.0,
                            libraries=None, met_limits=None, age_limits=None, idcs_ssp=None):

        """
        SESAMME ionization table: one row per metallicity (named as the to_fits extensions) and one column per
        age (named as the to_fits columns) with log10 Q. The rows of different libraries (e.g. VMS below 2.5 Myr
        and BPASS above) are merged per metallicity. Returns the metallicity x age grid as a DataFrame.
        """

        log_q = self.ionizing_photons(edge, flux_to_cgs, libraries, met_limits, age_limits, idcs_ssp)
        selection = self.frame.loc[log_q.index, ['metallicity', 'age', 'library']].assign(log_q=log_q)

        # One spectrum per grid node
        dup_mask = selection.duplicated(['metallicity', 'age'], keep=False)
        if dup_mask.any():
            dup_libs = sorted(selection.loc[dup_mask, 'library'].unique(), key=str)
            raise SpecSyError(f'{int(dup_mask.sum())} spectra share metallicity and age across the libraries '
                              f'{dup_libs}. Please select one library per age.')

        # Old populations without ionizing flux
        n_floor = int((selection['log_q'] < log_q_floor).sum())
        if n_floor > 0:
            _logger.warning(f'{n_floor} of the {len(selection)} spectra have log Q < {log_q_floor}. '
                            f'They were set to the log_q_floor={log_q_floor}.')
            selection['log_q'] = selection['log_q'].clip(lower=log_q_floor)

        grid = selection.pivot(index='metallicity', columns='age', values='log_q').sort_index().sort_index(axis=1)
        n_empty = int(grid.isna().to_numpy().sum())
        if n_empty > 0:
            raise SpecSyError(f'{n_empty} metallicity-age nodes have no spectrum. SESAMME needs the same ages at '
                              f'every metallicity: please adjust the selection.')

        ext_map = {met: _metallicity_to_key(met) for met in grid.index} if ext_map is None else ext_map
        table = Table([[ext_map[met] for met in grid.index]] + [grid[age].to_numpy() for age in grid.columns],
                      names=['Z'] + [str(float(age)) for age in grid.columns])
        table.write(fname, format='ascii', overwrite=True)

        return grid


    def save_frame(self, fname, **kwargs):

        """Save the properties frame with lime.save_frame (the keyword arguments are passed to it)."""

        save_frame(fname, self.frame, **kwargs)

        return


    def plot_age_metallicity(self, libraries=None, fig_cfg=None, fname=None):

        # Group the selected rows by library
        selection = self.frame.loc[self._library_mask(libraries, 'plot')]
        library_groups = selection.groupby('library', sort=False, dropna=False)
        use_colors = library_groups.ngroups > 1

        with rc_context(fig_cfg):

            fig, ax = plt.subplots(figsize=(8, 6))
            for i, (lib, group) in enumerate(library_groups):
                ax.scatter(group['age'], group['metallicity'], alpha=0.6, s=15,
                           color=f'C{i}' if use_colors else 'C0',
                           label=lib if use_colors else None)

            sources = sorted(self.frame['source'].unique(), key=str)
            ax.set_title(f"Age-Metallicity coverage — source: {', '.join(map(str, sources))}")
            ax.set_xlabel("log(Age/yr)")
            ax.set_ylabel("Metallicity (Z)")
            # ax.set_yscale('log')

            if use_colors:
                ax.legend(title='Library', markerscale=2)

            plt.tight_layout()
            if fname is None:
                plt.show()
            else:
                plt.savefig(fname)

        return


    def plot_ionizing_hardness(self, ages_myr, libraries=None, energies=None, ion_lines=MIRI_ION_LINES, z_sun=None,
                               ylim=None, lib_labels=None, fig_cfg=None, fname=None):

        # Loaded spectra of the requested libraries and ages
        selection = self.frame.loc[self._library_mask(libraries, 'plot') & self._loaded_mask()]
        grid_ages = np.unique(selection['age'])

        # Nearest grid age to each request
        ages = np.log10(np.where(np.atleast_1d(ages_myr) == 0, VMS_ZERO_AGE_MYR, np.atleast_1d(ages_myr)) * 1e6)
        nearest = grid_ages[np.abs(grid_ages[:, None] - ages[None, :]).argmin(axis=0)]
        missing = ages[np.abs(nearest - ages) > 0.02]
        if missing.size > 0:
            raise SpecSyError(f'No grid age within 0.02 dex of log(age)={missing.tolist()}. Available ages (Myr): '
                              f'{[float(f"{10 ** a / 1e6:.3g}") for a in grid_ages]}')

        # Unpack the values
        ages = np.unique(nearest)
        mets = np.unique(selection['metallicity'])
        libs = sorted(selection['library'].unique(), key=str)

        # Check if the library has the nebular continuum
        neb_libs = [lib for lib in selection['library'].unique() if str(lib).endswith('_neb')]
        if neb_libs:
            _logger.warning(f'The selection includes the nebular libraries {neb_libs}. Their flux below '
                            f'{Q_LIMITS_LAMBDA["HI"]:.1f} Å may not be the stellar ionizing output.')

        # Compute the ionization steps
        energies = (13.6, 60.01) if energies is None else energies
        energies = np.arange(energies[0], energies[1], 0.2)
        edges = HC_EV_ANGSTROM / energies

        # Figure format variables
        n_cols = min(len(ages), 2)
        n_rows = int(np.ceil(len(ages) / n_cols))

        colors = {met: f'C{i}' for i, met in enumerate(mets)}
        styles = {lib: ('-', '--', ':', '-.')[i % 4] for i, lib in enumerate(libs)}
        lib_labels = {} if lib_labels is None else lib_labels
        z_label = (lambda met: rf'$Z = {met / z_sun:.2g}\,Z_\odot$') if z_sun is not None else (lambda met: f'Z = {met}')

        # Generate the figure
        with rc_context(fig_cfg):

            fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, squeeze=False)

            # Check if the spectrum can produce ionization photons
            without_Q_counter = 0
            for ax, age in zip(axes.ravel(), ages):
                group = selection.loc[np.isclose(selection['age'], age)]
                for idx, row in group.iterrows():

                    wave_arr, flux_arr = self._spectrum(idx)
                    log_q0 = log_ionizing_photons(wave_arr, flux_arr, units_wave=self.units_wave, units_flux=self.units_flux)
                    if not np.isfinite(log_q0):
                        without_Q_counter += 1
                        continue

                    # Energies without photons are NaN, so the curve stops there
                    log_q = np.array([log_ionizing_photons(wave_arr, flux_arr, self.units_wave, self.units_flux,
                                                           lambda_limit=edge) for edge in edges])
                    ratio = np.where(np.isfinite(log_q), 10 ** (log_q - log_q0), np.nan)
                    ax.plot(energies, ratio, color=colors[row['metallicity']], linestyle=styles[row['library']])

                # Ionization energies of the diagnostic lines
                for label, energy in (ion_lines or {}).items():
                    if energies[0] <= energy <= energies[-1]:
                        ax.axvline(energy, color='black', linewidth=0.7, linestyle=':')
                        ax.text(energy, 0.03, label, rotation=90, ha='right', va='bottom', fontsize=8,
                                transform=ax.get_xaxis_transform())

                ax.set_yscale('log')
                ax.set_xlim(energies[0], energies[-1])
                if ylim is not None:
                    ax.set_ylim(*ylim)
                ax.set_title(f'{10 ** age / 1e6:.3g} Myr')

            if without_Q_counter > 0:
                _logger.warning(f'{without_Q_counter} spectra have no flux below {Q_LIMITS_LAMBDA["HI"]:.1f} Å and were not plotted.')

            # Unused panels
            for ax in axes.ravel()[len(ages):]: ax.set_visible(False)

            # Label the axes
            for ax in axes[-1, :]: ax.set_xlabel('Energy, E [eV]')
            for ax in axes[:, 0]: ax.set_ylabel(r'$Q(>E)\,/\,Q(>13.6\,\mathrm{eV})$')

            # Colors for the metallicity and line styles for the library (if several)
            handles = [Line2D([], [], color=colors[met], label=z_label(met)) for met in mets]
            if len(libs) > 1:
                handles += [Line2D([], [], color='0.3', linestyle=styles[lib], label=lib_labels.get(lib, lib))
                            for lib in libs]

            # Legend in first plot
            axes[0, 0].legend(handles=handles, loc='upper right', fontsize='small')

            plt.tight_layout()
            if fname is None:
                plt.show()
            else:
                plt.savefig(fname, bbox_inches='tight')
                plt.close(fig)

        return

    def plot_spectra(self, fname=None, metallicity=None, age=None, metallicity_range=None, age_range=None,
                     libraries=None, log_scale=False, in_fig=_NO_FIG, fig_cfg=None, ax_cfg=None, maximize=False):

        # Display check for input figures
        display_check = True if in_fig is _NO_FIG else False

        # Select loaded rows matching the filters (exact values and inclusive ranges)
        mask = self._library_mask(libraries, 'plot') & self._loaded_mask()
        if metallicity is not None:
            mask &= self.frame['metallicity'].isin(np.atleast_1d(metallicity)).to_numpy()
        if metallicity_range is not None:
            mask &= self.frame['metallicity'].between(*metallicity_range).to_numpy()
        if age is not None:
            mask &= self.frame['age'].isin(np.atleast_1d(age)).to_numpy()
        if age_range is not None:
            mask &= self.frame['age'].between(*age_range).to_numpy()

        if not mask.any():
            raise SpecSyError('No loaded binaries match the requested filters.')

        selection = self.frame.loc[mask].sort_values(['metallicity', 'age'])

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
            unique_mets = np.unique(selection['metallicity'])
            met_color = {m: f'C{i}' for i, m in enumerate(unique_mets)}

            legend_entries = {}
            for idx, row in selection.iterrows():
                wave_arr, flux_arr = self._spectrum(idx)
                label_key = f'Z={row["metallicity"]:.5f} | {row["library"]}'
                self.ax.step(wave_arr, flux_arr, where='mid', alpha=0.7,
                             linewidth=0.8, color=met_color[row['metallicity']],
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


    def _single_spectrum(self, idcs, age, metallicity):

        """Wavelength and flux of the single row in `idcs`, with the multiple-match and loaded checks."""

        if len(idcs) > 1:
            subset = self.frame.loc[idcs]
            varying = {col: subset[col].unique().tolist() for col in subset.columns
                       if (col != 'fpath') and (subset[col].nunique(dropna=False) > 1)}
            raise SpecSyError(f'Multiple stellar spectra ({len(idcs)}) match age={age}, metallicity={metallicity}. '
                              f'Please refine with: {varying if varying else "none (duplicated rows)"}.')

        wave_arr, flux_arr = self._spectrum(idcs[0])
        if flux_arr is None:
            raise SpecSyError(f'The stellar spectrum with age={age}, metallicity={metallicity} is not loaded. '
                              f'{_NOT_LOADED_MSG}')

        return wave_arr, flux_arr


    def _set_data(self, frame, wave_list, flux_list):

        # Attributes
        self.frame = frame.reset_index(drop=True)
        self.flux_series = None
        self.wave_series = None
        self.units_wave = au.Unit('AA')
        self.units_flux = au.Unit('LLAM')

        # Spectra series sharing the frame index (None if no spectrum is loaded)
        if (flux_list is not None) and (wave_list is not None):
            n_rows = len(self.frame)
            n_loaded = sum(flux is not None for flux in flux_list)
            if n_loaded > 0:
                if n_loaded < n_rows:
                    _logger.warning(f'{n_rows - n_loaded} of the {n_rows} stellar spectra are not loaded.')
                self.flux_series = pd.Series(flux_list, index=self.frame.index, dtype=object)
                self.wave_series = pd.Series(wave_list, index=self.frame.index, dtype=object)
        else:
            if not ((flux_list is None) and (wave_list is None)):
                msg = "flux_list has values but the wave_list does not" if flux_list is not None else "wave_list has values but the flux_list does not."
                raise SpecSyError(f'The input {msg}')

        # Unique metallicities and ages
        self.metallicities = np.unique(self.frame['metallicity'])
        self.ages = np.unique(self.frame['age'])

        # Source and library assignment
        unique_sources = sorted(self.frame['source'].unique(), key=str)
        unique_libraries = sorted(self.frame['library'].unique(), key=str)

        self.source = unique_sources[0] if len(unique_sources) == 1 else unique_sources
        self.library = unique_libraries[0] if len(unique_libraries) == 1 else unique_libraries

        if len(unique_sources) > 1:
            _logger.warning(f'StellarBinaries contains multiple sources: {unique_sources}')
        if len(unique_libraries) > 1:
            _logger.warning(f'StellarBinaries contains multiple libraries: {unique_libraries}')

        # Rows with the same value in every column
        n_duplicated = int(self.frame.duplicated().sum())
        if n_duplicated > 0:
            _logger.warning(f'{n_duplicated} of the {n_rows} rows in the StellarBinaries frame repeat all the values '
                            f'of another row. These spectra cannot be told apart in a selection.')

        # Check wavelength uniformity across loaded spectra
        self._check_dispersion_uniformity()

        return


    @staticmethod
    def _bracket(grid, value, log=False):

        """Indices of the grid nodes around `value` and the weight of the upper one (linear or log10 scale)."""

        if grid.size == 1:
            return 0, 0, 0.0

        x_grid, x = (np.log10(grid), np.log10(value)) if log else (grid, value)
        i = int(np.clip(np.searchsorted(x_grid, x, side='right') - 1, 0, grid.size - 2))

        return i, i + 1, float((x - x_grid[i]) / (x_grid[i + 1] - x_grid[i]))


    def _check_dispersion_uniformity(self):

        self.uniform_wmin = None
        self.uniform_wmax = None
        self.uniform_deltalambda = None
        self.uniform_dispersion = False

        if self.wave_series is None:
            return

        # Wavelength uniformity
        loaded = [wave for wave in self.wave_series if wave is not None]
        wmins = np.array([wave.min() for wave in loaded])
        wmaxs = np.array([wave.max() for wave in loaded])
        deltas = np.array([np.mean(np.diff(wave)) for wave in loaded])

        self.uniform_wmin = float(wmins[0]) if np.all(np.isclose(wmins, wmins[0])) else None
        self.uniform_wmax = float(wmaxs[0]) if np.all(np.isclose(wmaxs, wmaxs[0])) else None
        self.uniform_deltalambda = float(deltas[0]) if np.all(np.isclose(deltas, deltas[0])) else None

        self.uniform_dispersion = ((self.uniform_wmin is not None) and (self.uniform_wmax is not None)
                                   and (self.uniform_deltalambda is not None))

        non_uniform = [k for k, v in [('wmin', self.uniform_wmin), ('wmax', self.uniform_wmax),
                                      ('delta_lambda', self.uniform_deltalambda)] if v is None]
        if non_uniform:
            _logger.warning(f'Wavelength grid is not uniform across loaded binaries for: {non_uniform}')

        # A uniform dispersion is stored once, with index 0
        if self.uniform_dispersion:
            self.wave_series = pd.Series([loaded[0]], index=[0], dtype=object)

        return


    def _loaded_mask(self):
        """Boolean array with the frame rows whose spectrum is loaded."""
        if self.flux_series is None:
            return np.zeros(len(self.frame), dtype=bool)
        return np.array([flux is not None for flux in self.flux_series])


    def _spectrum(self, idx):

        """Wavelength and flux arrays for the frame row `idx` (None, None if not loaded)."""

        if (self.flux_series is None) or (self.flux_series.at[idx] is None):
            return None, None

        wave_arr = self.wave_series.at[0] if self.uniform_dispersion else self.wave_series.at[idx]

        return wave_arr, self.flux_series.at[idx]


    def _select(self, **params):

        """
        Frame index of the rows matching the requested column values. The values must exist in the column
        (exact match). A list selects several values and None skips the column.
        """

        unknown = set(params) - set(self.frame.columns)
        if unknown:
            raise SpecSyError(f'Unrecognised parameters {sorted(unknown)}. Available: {list(self.frame.columns)}')

        mask = np.ones(len(self.frame), dtype=bool)
        for col, value in params.items():
            if value is None:
                continue

            values = [value] if np.ndim(value) == 0 else list(value)
            column = self.frame[col]

            missing = [v for v in values if not (column == v).any()]
            if missing:
                raise SpecSyError(f'{col}={missing} not found. Available values: {column.unique().tolist()}')

            mask &= column.isin(values).to_numpy()

        if not mask.any():
            requested = {k: v for k, v in params.items() if v is not None}
            raise SpecSyError(f'No stellar spectra match the combination {requested}.')

        return self.frame.index[mask]


    def _library_mask(self, libraries, action):

        """Boolean array with the frame rows of the requested libraries (all if None)."""

        all_libraries = sorted(self.frame['library'].unique(), key=str)
        if libraries is None:
            return np.ones(len(self.frame), dtype=bool)

        requested = [libraries] if isinstance(libraries, str) else list(libraries)
        unrecognised = set(requested) - set(all_libraries)
        if unrecognised:
            _logger.warning(f'Libraries not found in StellarBinaries: {unrecognised}. '
                            f'Available libraries: {all_libraries}')

        requested = [lib for lib in requested if lib in all_libraries]
        if len(requested) == 0:
            raise SpecSyError(f'No valid libraries to {action}. Available libraries: {all_libraries}')

        return self.frame['library'].isin(requested).to_numpy()


    def _limits_mask(self, column, limits):

        """Boolean array with the frame rows where low < value < high (a None bound is unbounded)."""

        values = self.frame[column].to_numpy()
        mask = np.ones(values.size, dtype=bool)
        if limits is not None:
            low, high = limits
            if low is not None:
                mask &= values > low
            if high is not None:
                mask &= values < high

        return mask


    def _selection_index(self, idcs_ssp, libraries, met_limits, age_limits):

        """Frame index of the selected rows. The `idcs_ssp` argument takes precedence over the other ones."""

        # Direct indexing with frame indices or a boolean mask
        if idcs_ssp is not None:

            ignored = [name for name, value in (('libraries', libraries), ('met_limits', met_limits),
                                                ('age_limits', age_limits)) if value is not None]
            if ignored:
                _logger.warning(f'The "idcs_ssp" argument takes precedence: {ignored} are ignored.')

            idcs = np.atleast_1d(np.asarray(idcs_ssp))
            if idcs.dtype == bool:
                if idcs.size != len(self.frame):
                    raise SpecSyError(f'A boolean "idcs_ssp" must have one entry per frame row ({len(self.frame)}), '
                                      f'received {idcs.size}.')
                return self.frame.index[idcs]

            if (idcs.size > 0) and (idcs.dtype.kind not in 'iu'):
                raise SpecSyError(f'The "idcs_ssp" must be integer frame indices or a boolean mask, received dtype '
                                  f'{idcs.dtype}.')

            missing = idcs[~np.isin(idcs, self.frame.index)]
            if missing.size > 0:
                raise SpecSyError(f'The "idcs_ssp" indices {missing.tolist()} are not in the frame '
                                  f'(0 to {len(self.frame) - 1}).')

            return pd.Index(idcs)

        # Check the age and metallicity limits
        for name, limits in (('met_limits', met_limits), ('age_limits', age_limits)):
            if limits is None:
                continue
            if len(limits) != 2:
                raise SpecSyError(f'The "{name}" argument must be a (low, high) tuple, received: {limits}')
            low, high = limits
            if (low is not None) and (high is not None) and (low >= high):
                raise SpecSyError(f'The "{name}" lower value must be below the higher one, received: {limits}')

        # Select the rows by library and by the age/metallicity limits
        mask = (self._library_mask(libraries, 'select') & self._limits_mask('metallicity', met_limits)
                & self._limits_mask('age', age_limits))

        return self.frame.index[mask]


    def _process_selection(self, libraries, wmin, wmax, met_limits, age_limits, disp_intvl, pixel_width,
                           pixel_number, constant_pixel_width, idcs_ssp):

        """
        Selected frame rows with their wavelength and flux arrays, cropped to the wmin-wmax interval and resampled
        with LiMe's spectrum_resampling if requested. Common to to_fits and to_stellar_binaries.
        """

        # At most one resampling criterion
        rebin_args = {'disp_intvl': disp_intvl, 'pixel_width': pixel_width, 'pixel_number': pixel_number}
        active_rebin = [k for k, v in rebin_args.items() if v is not None]
        if len(active_rebin) > 1:
            raise SpecSyError(f'Rebinning arguments {active_rebin} are mutually exclusive. Please provide only one.')

        # Rows selection
        idcs = self._selection_index(idcs_ssp, libraries, met_limits, age_limits)
        if len(idcs) == 0:
            selection_str = ('idcs_ssp' if idcs_ssp is not None else
                             f'libraries={libraries}, met_limits={met_limits}, age_limits={age_limits}')
            raise SpecSyError(f'No spectra left after the selection ({selection_str}). Available metallicities: '
                              f'{list(self.metallicities)}, available ages: {list(self.ages)}.')

        # The spectra must be in memory
        n_missing = int((~self._loaded_mask()[idcs.to_numpy()]).sum())
        if n_missing > 0:
            raise SpecSyError(f'{n_missing} of the {len(idcs)} selected spectra are not loaded. {_NOT_LOADED_MSG}')

        # Crop to the wavelength interval
        wave_list, flux_list = [], []
        for idx in idcs:
            wave_arr, flux_arr = self._crop_spectrum(*self._spectrum(idx), wmin, wmax)
            if wave_arr.size < 2:
                row = self.frame.loc[idx]
                raise SpecSyError(f'Fewer than 2 pixels left for age={row["age"]}, metallicity={row["metallicity"]} '
                                  f'after the wmin={wmin}, wmax={wmax} crop.')
            wave_list.append(wave_arr)
            flux_list.append(flux_arr)

        # Resample on a common output dispersion, computed once from the first spectrum
        if len(active_rebin) == 1:

            target_grid = None if disp_intvl is None else np.asarray(disp_intvl, dtype=float)
            if target_grid is None:
                target_grid, _, _ = spectrum_resampling(None, pixel_width, pixel_number, constant_pixel_width,
                                                        wave_list[0], flux_list[0], None, None)

            empty_bins = 0
            for i, (wave_arr, flux_arr) in enumerate(zip(wave_list, flux_list)):
                new_wave, new_flux, _ = spectrum_resampling(target_grid, pixel_width, pixel_number,
                                                            constant_pixel_width, wave_arr, flux_arr, None, None)
                wave_list[i] = np.asarray(new_wave, dtype=float)
                flux_list[i] = np.asarray(new_flux, dtype=float)
                empty_bins += int(np.isnan(flux_list[i]).sum())

            if empty_bins > 0:
                _logger.warning(f'{empty_bins} output bins received no input pixels and were filled with NaN. '
                                f'The output dispersion is finer than the input one over part of the range.')

        return self.frame.loc[idcs], wave_list, flux_list