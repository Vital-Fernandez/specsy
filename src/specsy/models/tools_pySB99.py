import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import pandas as pd
from pathlib import Path
from astropy.table import Table
from specsy import cfg as sy_cfg


def find_ionizing_files(root_folder, mets=None, ion_file='ion_flux_HI.txt', timesteps_file='timesteps.txt'):

    """
    Locate the ionizing output and age step files in a compiled spectra root folder.

    Parameters
    ----------
    root_folder : str or Path
        Folder holding one subfolder per metallicity.
    mets : list, optional
        Metallicity labels declaring the output order. If None the subfolder names are
        sorted alphabetically. Otherwise every label must match a subfolder name.
    ion_file : str
        Name of the log(Q(HI)) file within each metallicity folder.
    timesteps_file : str
        Name of the age steps file within each metallicity folder.

    Returns
    -------
    file_list, mets, timesteps_path
    """

    root_folder = Path(root_folder)

    # Metallicity folders containing the ionizing output
    folder_dict = {folder.name: folder for folder in sorted(root_folder.iterdir())
                   if folder.is_dir() and (folder/ion_file).is_file()}

    if len(folder_dict) == 0:
        raise ValueError(f'No subfolder with a "{ion_file}" file at {root_folder}')

    # The requested order, or the alphabetical one
    list_folders = list(folder_dict)
    if mets is None:
        mets = [sy_cfg['stellar']['ssp']['pystarburst99']['z_keys'][folder.split('_')[0]]
                for folder in list_folders]
    else:
        missing = [metal for metal in mets if metal not in folder_dict]
        if len(missing) > 0:
            raise ValueError(f'Metallicities without a "{ion_file}" file: {missing}')
        list_folders = list(mets)
        mets = [sy_cfg['stellar']['ssp']['pystarburst99']['z_keys'][folder.split('_')[0]]
                for folder in list_folders]

    # The age steps must be common to every metallicity
    file_list = [root_folder/folder/ion_file for folder in list_folders]
    timesteps_list = [root_folder/folder/timesteps_file for folder in list_folders]

    missing = [path.parent.name for path in timesteps_list if not path.is_file()]
    if len(missing) > 0:
        raise ValueError(f'Metallicities without a "{timesteps_file}" file: {missing}')

    ages_ref = np.loadtxt(timesteps_list[0])
    for metal, path in zip(mets[1:], timesteps_list[1:]):
        ages_i = np.loadtxt(path)
        if (ages_i.shape != ages_ref.shape) or not np.allclose(ages_i, ages_ref):
            raise ValueError(f'The "{timesteps_file}" of {metal} does not match {mets[0]}')

    return file_list, mets, timesteps_list[0]


def ionizing_table(file_list, mets, timesteps_file, output_file=None):

    """
    Combine single-column log(Q(HI)) files into one table, with one metallicity per row
    and one column per age step. Rows are sorted by increasing metallicity.

    Parameters
    ----------
    file_list : list
        Paths to the ionizing output files, one per metallicity, each a single ascii
        column of log(Q(HI)) values ordered by age step.
    mets : list
        Numerical metallicity values, in the same order as ``file_list``.
    timesteps_file : str or Path
        Ascii file with the grid age steps in logarithmic scale, one per line.
    output_file : str or Path, optional
        Destination for the ascii table. The table is only written if this is provided.
    """

    if len(file_list) != len(mets):
        raise ValueError(f'Received {len(file_list)} files for {len(mets)} metallicities')

    # Age steps labelling the columns
    ages = np.loadtxt(timesteps_file)

    # One row per metallicity, one column per age step
    mets = np.asarray(mets, dtype=float)
    q_matrix = np.array([np.loadtxt(file) for file in file_list])

    if q_matrix.shape != (mets.size, ages.size):
        raise ValueError(f'Ionizing data has shape {q_matrix.shape}, expected {(mets.size, ages.size)}')

    # Sort the rows by increasing metallicity
    idcs_sort = np.argsort(mets)
    mets, q_matrix = mets[idcs_sort], q_matrix[idcs_sort, :]

    # Metallicity column followed by the age columns
    nion_table = Table()
    nion_table['Z'] = mets
    for i, age in enumerate(ages):
        nion_table[str(age)] = q_matrix[:, i]

    if output_file is not None:
        nion_table.write(output_file, format='ascii', overwrite=True)

    return nion_table

def calc_Nostars(IMF_masses, IMF_exponents, IMF_mass_limits, M_total):
    if len(IMF_exponents) == 1:
        N_IMF_intervals = 1
    elif len(IMF_exponents) > 1:
        N_IMF_intervals = 2

    if N_IMF_intervals > 1:
        N_IMF_intervals = len(IMF_exponents)
        A_ic = np.zeros_like(IMF_exponents)
        A_ic[0] = 1.

        for exponent_index in range(len(IMF_exponents)):
            exponent_index = exponent_index + 1
            if exponent_index == len(IMF_exponents):
                break
            A_i = A_ic[exponent_index-1] * (IMF_mass_limits[exponent_index]**(IMF_exponents[exponent_index] - IMF_exponents[exponent_index-1]))
            A_ic[exponent_index] = A_i
        k_ic = []
        for exponent_index in range(len(IMF_exponents)):
            if IMF_exponents[exponent_index] == 2.0:
                k_i = np.log(IMF_mass_limits[exponent_index+1]) - np.log(IMF_mass_limits[exponent_index])
            else:
                k_i = (IMF_mass_limits[exponent_index+1]**(2.0 - IMF_exponents[exponent_index]) - IMF_mass_limits[exponent_index]**(2.0 - IMF_exponents[exponent_index])) / (2 - IMF_exponents[exponent_index])
            k_ic.append(k_i)

        Ak_ic = A_ic * k_ic
        Ak = sum(Ak_ic)

    else:
        A_ic = 1
        k_ic = (IMF_mass_limits[1]**(2.0 - IMF_exponents[0]) - IMF_mass_limits[0]**(2.0 - IMF_exponents[0])) / (2 - IMF_exponents[0])
        Ak = A_ic * k_ic

    S = M_total / Ak
    A = A_ic * S

    xmhigh = np.full_like(IMF_masses, 0)
    xmlow = np.full_like(IMF_masses, 0)

    for mass_index in range(len(IMF_masses)):
        if IMF_masses[mass_index] == min(IMF_masses):
            xmhigh[mass_index] = 0.5 * (IMF_masses[mass_index] + IMF_masses[mass_index+1])
            xmlow[mass_index] = IMF_masses[mass_index]
        elif IMF_masses[mass_index] == max(IMF_masses):
            xmhigh[mass_index] = IMF_masses[mass_index]
            xmlow[mass_index] = 0.5 * (IMF_masses[mass_index-1] + IMF_masses[mass_index])
        else:
            xmhigh[mass_index] = 0.5 * (IMF_masses[mass_index] + IMF_masses[mass_index+1])
            xmlow[mass_index] = 0.5 * (IMF_masses[mass_index] + IMF_masses[mass_index-1])

    dens = []
    if N_IMF_intervals == 2:
        for mass_index in range(len(IMF_masses)):
            dens_i = A[1] * (xmhigh[mass_index]**(1-IMF_exponents[1]) - xmlow[mass_index]**(1-IMF_exponents[1])) / (1 - IMF_exponents[1])
            dens.append(dens_i)
    else:
        for mass_index in range(len(IMF_masses)):
            dens_i = A * (xmhigh[mass_index]**(1-IMF_exponents[0]) - xmlow[mass_index]**(1-IMF_exponents[0])) / (1 - IMF_exponents[0])
            dens.append(dens_i)

    N_stars = np.array(dens)
    total_mass = N_stars * IMF_masses
    return N_stars, IMF_masses, dens, xmhigh, xmlow


def read_spectra_grid(spec_file):
    df = pd.read_csv(spec_file)
    ind_spec = [i for i in range(len(df)) if df['start'][i] == ' CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC']
    spec_params = []
    spectra = []
    for i in range(len(ind_spec)-1):
        spec_param = df['start'][ind_spec[i]+1]
        spec_params.append(spec_param)
    for i in range(len(ind_spec)-1):
        spec_param = df['start'][ind_spec[i]+2: ind_spec[i+1]]
        spectra.append(spec_param)
    return spec_params, spectra


def reformat_spec(spectrum):
    wave_reform = []
    flux_reform = []
    for i in spectrum:
        spec_wavep, spec_fluxp = i.split()
        wave_reform.append(float(spec_wavep))
        flux_reform.append(float(spec_fluxp) * 4.0 * 3.14142)
    flux_reform = np.array(flux_reform)
    return np.column_stack((wave_reform, flux_reform))


def reformat_spec_WR(spectrum):
    wave_reform = []
    flux_reform = []
    for i in spectrum:
        spec_wavep, spec_fluxp = i.split()
        wave_reform.append(float(spec_wavep))
        flux_reform.append(float(spec_fluxp))
    flux_reform = np.array(flux_reform)
    return np.column_stack((wave_reform, flux_reform))


def reform_spec_grid(spectra, spec_params):
    reformed_spec_grid = [reformat_spec(spectra[i]) for i in range(len(spectra))]
    spec_params_id, spec_params_teff, spec_params_logl, spec_params_logg = [], [], [], []
    for i in spec_params:
        spec_id, spec_teff, spec_logl, spec_logg = i.split()
        spec_params_id.append(float(spec_id))
        spec_params_teff.append(float(spec_teff))
        spec_params_logl.append(float(spec_logl))
        spec_params_logg.append(float(spec_logg))
    spec_params_teff = 10**np.array(spec_params_teff)
    spec_params_reform = np.column_stack((spec_params_id, spec_params_teff, spec_params_logl, spec_params_logg))
    return reformed_spec_grid, spec_params_reform, spec_params_teff, spec_params_logl, spec_params_logg


def reform_spec_grid_WR(spectra, spec_params):
    reformed_spec_grid = [reformat_spec_WR(spectra[i]) for i in range(len(spectra))]
    spec_params_id, spec_params_teff = [], []
    for i in spec_params:
        spec_id, spec_teff = i.split()
        spec_params_id.append(float(spec_id))
        spec_params_teff.append(float(spec_teff))
    spec_params_teff = np.array(spec_params_teff)
    spec_params_reform = np.column_stack((spec_params_id, spec_params_teff))
    return reformed_spec_grid, spec_params_reform, spec_params_teff


def reform_spec_grid_powr(spectra, spec_params):
    reformed_spec_grid = [reformat_spec_WR(spectra[i]) for i in range(len(spectra))]
    spec_params_id, spec_params_teff, spec_params_radius, spec_params_length = [], [], [], []
    for i in spec_params:
        spec_id, spec_teff, spec_radius, spec_length = i.split()
        spec_params_id.append(float(spec_id))
        spec_params_teff.append(float(spec_teff))
        spec_params_radius.append(float(spec_radius))
        spec_params_length.append(float(spec_length))
    spec_params_teff = np.array(spec_params_teff)
    spec_params_reform = np.column_stack((spec_params_id, spec_params_teff, spec_params_radius, spec_params_length))
    return reformed_spec_grid, spec_params_reform, spec_params_teff, spec_params_radius, spec_params_length


def integrate_spec_grid(reformed_spec_grid):
    return [np.trapezoid(reformed_spec_grid[i][:, 1], reformed_spec_grid[i][:, 0]) for i in range(len(reformed_spec_grid))]


def integrate_spec_grid_powr(wave_grid, reformed_spec_grid):
    return [np.trapz(wave_grid, reformed_spec_grid[i, :]) for i in range(len(reformed_spec_grid))]


def interpolate_param(tracks_parameter, track_masses, run_speed_mode):
    grid_masses_adjinterp_total = []
    grid_params_adjinterp_total = []
    for i in range(len(track_masses)-1):
        track_mass_upper = np.flip(track_masses[i])
        track_param_upper = np.flip(tracks_parameter[i])
        track_mass_lower = np.flip(track_masses[i+1])
        track_param_lower = np.flip(tracks_parameter[i+1])

        initial_mass_upper = track_mass_upper[-1]
        initial_mass_lower = track_mass_lower[-1]

        if run_speed_mode == 'DEFAULT':
            if round(initial_mass_upper) > 9 and round(initial_mass_upper) < 15:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 15
            elif round(initial_mass_upper) > 15 and round(initial_mass_upper) < 35:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 5
            elif round(initial_mass_upper) < 9:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 10
            elif round(initial_mass_upper) > 35 and round(initial_mass_upper) < 60:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 5
            elif round(initial_mass_upper) > 60:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1)
        if run_speed_mode == 'FAST':
            inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 1
        if run_speed_mode == 'HIGH_RES':
            if round(initial_mass_upper) > 7 and round(initial_mass_upper) < 35:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 100
            else:
                inter_track_sampling = (initial_mass_upper - initial_mass_lower + 1) * 50

        track_masses_adjinterp = np.column_stack((track_mass_lower, track_mass_upper))
        track_param_adjinterp = np.column_stack((track_param_lower, track_param_upper))

        grid_masses_adjinterp = [np.linspace(track_mass_lower[j], track_mass_upper[j], round(inter_track_sampling))
                                  for j in range(len(track_mass_upper))]

        grid_params_adjinterp = []
        for track_param, track_mass, grid_mass in zip(track_param_adjinterp, track_masses_adjinterp, grid_masses_adjinterp):
            grid_params_adjinterp_fn = interp1d(track_mass, track_param, kind='linear', fill_value='extrapolate')
            grid_params_adjinterp.append(grid_params_adjinterp_fn(grid_mass))

        grid_masses_adjinterp_total.append(grid_masses_adjinterp)
        grid_params_adjinterp_total.append(grid_params_adjinterp)

    return grid_masses_adjinterp_total, grid_params_adjinterp_total


def rearrange_grid_array(grid_array):
    lengths = [len(grid_array[i][0]) for i in range(len(grid_array))]
    total_grid_size = sum(lengths)
    rearranged = []
    for k in range(len(grid_array)):
        for j in range(len(grid_array[k][0])-1, -1, -1):
            for i in range(len(grid_array[k])):
                rearranged.append(grid_array[k][i][j])
    return np.array_split(rearranged, total_grid_size)


def get_timestep_0_ind(timestep, grid_ages, grid_masses, times_steps, IMF_mass_limits):
    timestep_masses = []
    for i in range(len(grid_ages)):
        ind_nearest_age = np.argmin((grid_ages[i] - times_steps[0])**2)
        timestep_masses.append(grid_masses[i][ind_nearest_age])
    timestep_masses_unique = np.unique(timestep_masses, return_index=True)
    timestep_masses_IMF_ind = [i for i in range(len(timestep_masses_unique[0])) if timestep_masses_unique[0][i] < IMF_mass_limits[-1]]
    return timestep_masses_unique[1][timestep_masses_IMF_ind]


def get_timestep_params(timestep, timestep_mass_ind, grid_ages, grid_masses, grid_temps, grid_lums,
                        grid_H_abundances, grid_He_abundances, grid_12C_abundances, grid_14N_abundances,
                        grid_16O_abundances, grid_core_temps, grid_mass_loss_rates):
    timestep_ages, timestep_masses, timestep_temps, timestep_lums = [], [], [], []
    timestep_H_abundances, timestep_He_abundances = [], []
    timestep_12C_abundances, timestep_14N_abundances, timestep_16O_abundances = [], [], []
    timestep_core_temps, timestep_mass_loss_rates = [], []

    for i in range(len(grid_ages)):
        ind_nearest_age = np.argmin(abs(grid_ages[i] - timestep))
        timestep_ages.append(grid_ages[i][ind_nearest_age])
        timestep_masses.append(grid_masses[i][ind_nearest_age])
        timestep_temps.append(grid_temps[i][ind_nearest_age])
        timestep_lums.append(grid_lums[i][ind_nearest_age])
        timestep_H_abundances.append(grid_H_abundances[i][ind_nearest_age])
        timestep_He_abundances.append(grid_He_abundances[i][ind_nearest_age])
        timestep_12C_abundances.append(grid_12C_abundances[i][ind_nearest_age])
        timestep_14N_abundances.append(grid_14N_abundances[i][ind_nearest_age])
        timestep_16O_abundances.append(grid_16O_abundances[i][ind_nearest_age])
        timestep_core_temps.append(grid_core_temps[i][ind_nearest_age])
        timestep_mass_loss_rates.append(grid_mass_loss_rates[i][ind_nearest_age])

    timestep_ages_final, timestep_masses_final, timestep_temps_final, timestep_lums_final = [], [], [], []
    timestep_H_abundances_final, timestep_He_abundances_final = [], []
    timestep_12C_abundances_final, timestep_14N_abundances_final, timestep_16O_abundances_final = [], [], []
    timestep_core_temps_final, timestep_mass_loss_rates_final = [], []

    for i in timestep_mass_ind:
        timestep_ages_final.append(timestep_ages[i])
        timestep_temps_final.append(timestep_temps[i])
        timestep_lums_final.append(timestep_lums[i])
        timestep_masses_final.append(timestep_masses[i])
        timestep_H_abundances_final.append(timestep_H_abundances[i])
        timestep_He_abundances_final.append(timestep_He_abundances[i])
        timestep_12C_abundances_final.append(timestep_12C_abundances[i])
        timestep_14N_abundances_final.append(timestep_14N_abundances[i])
        timestep_16O_abundances_final.append(timestep_16O_abundances[i])
        timestep_core_temps_final.append(timestep_core_temps[i])
        timestep_mass_loss_rates_final.append(timestep_mass_loss_rates[i])

    timestep_teffs_calc  = 10**np.array(timestep_temps_final)
    timestep_lum_calc    = 10**np.array(timestep_lums_final) * 3.839e33
    radii_denom          = 12.566 * 5.670e-5 * timestep_teffs_calc**4
    timestep_radii_calc  = ((timestep_lum_calc / radii_denom)**0.5) * 100
    logg_num             = 6.67e-12 * np.array(timestep_masses_final) * 1.9891e33 * 1e8
    timestep_loggs_final = np.log10(logg_num / timestep_radii_calc**2)

    timestep_H_abundances_final  = np.array(timestep_H_abundances_final)
    timestep_He_abundances_final = np.array(timestep_He_abundances_final)
    timestep_12C_abundances_final= np.array(timestep_12C_abundances_final)
    timestep_14N_abundances_final= np.array(timestep_14N_abundances_final)
    timestep_16O_abundances_final= np.array(timestep_16O_abundances_final)

    timestep_cnr  = timestep_12C_abundances_final / timestep_14N_abundances_final
    timestep_coher= ((timestep_12C_abundances_final/12.) + (timestep_16O_abundances_final/16.)) / (timestep_He_abundances_final/4.)

    WR_correction_factor = 0.6
    for i in range(len(timestep_masses_final)):
        if timestep_H_abundances_final[i] < 0.1:
            timestep_core_teff_final = 10**timestep_core_temps_final[i]
            corrected_teff = timestep_core_teff_final + (WR_correction_factor - 1.0) * (timestep_core_teff_final - timestep_teffs_calc[i])
            timestep_temps_final[i] = np.log10(corrected_teff)

    timestep_temps_final = np.array(timestep_temps_final)
    return (timestep_ages_final, timestep_temps_final, timestep_lums_final, timestep_masses_final,
            timestep_H_abundances_final, timestep_loggs_final, timestep_mass_loss_rates_final,
            timestep_12C_abundances_final, timestep_14N_abundances_final, timestep_16O_abundances_final,
            timestep_masses, timestep_cnr, timestep_coher)


def get_specsyn_params(timestep_temps_final, timestep_masses_final, timestep_lums_final):
    specsyn_teffs, specsyn_loggs, specsyn_radii, specsyn_cotests, specsyn_bbfluxes = [], [], [], [], []
    for i in range(len(timestep_temps_final)):
        teff_specsyn   = 10**(timestep_temps_final[i])
        logg_specsyn   = np.log10(timestep_masses_final[i]) + 4.*timestep_temps_final[i] - timestep_lums_final[i] - 10.6
        radius_specsyn = 10.**(10.8426 + 0.5*timestep_lums_final[i] - 2.*timestep_temps_final[i] + 7.52)
        cotest         = 5.71*np.log10(teff_specsyn) - 21.95
        specsyn_bbflux = 0.0 if teff_specsyn == 0.0 else 5.6696196e-05 * teff_specsyn**4.
        specsyn_teffs.append(teff_specsyn)
        specsyn_loggs.append(logg_specsyn)
        specsyn_radii.append(radius_specsyn)
        specsyn_cotests.append(cotest)
        specsyn_bbfluxes.append(specsyn_bbflux)
    specsyn_teffs = np.array(specsyn_teffs)
    ind_SN_stars  = [i for i in range(len(specsyn_teffs)) if specsyn_teffs[i] == 1.0]
    specsyn_teffs[ind_SN_stars] = 0.0
    return specsyn_teffs, specsyn_loggs, specsyn_radii, specsyn_cotests, specsyn_bbfluxes


def assign_spectra_to_grid_WR(wave_grid, empty_flux, minimum_wr_mass, WC_spec_params_reform, WN_spec_params_reform,
                              WN_spec_params_teff, WN_spec_params_logl, WN_reformed_spec_grid, WN_integrated_spectra,
                              WC_spec_params_teff, WC_spec_params_logl, WC_reformed_spec_grid, WC_integrated_spectra,
                              No_stars, reformed_spec_grid, integrated_spectra, lowmass_params, timestep_loggs_final,
                              lowmass_teffs, lowmass_loggs, lm_spec_params_logl, lowmass_spec, lowmass_int_spec,

                              timestep_temps_final,spec_params_reform, spec_params_teff, spec_params_logl, timestep_lums_final,
                              timestep_H_abundances_final, spec_params_logg, timestep_masses_final, initial_masses,
                              specsyn_cotests, specsyn_loggs, timestep_cnr, timestep_coher):

    assigned_integrated_spectra=[]
    assigned_spec_teff=[]
    assigned_spec_logl=[]
    assigned_spec_logg=[]
    assigned_spectra=[]
    timestep_teffs_final = 10**timestep_temps_final
    population_choice = []

    for j in range(len(timestep_temps_final)):
        distance_to_spec = []

        if timestep_lums_final[j] == -20.0:
            population_choice.append('nope')
            assigned_spec_teff.append(0.0)
            assigned_spectra.append(np.column_stack((wave_grid, empty_flux)))
            assigned_integrated_spectra.append(1e-30)

        elif timestep_temps_final[j] > 4.4 and timestep_H_abundances_final[j] < 0.4 and initial_masses[j] > minimum_wr_mass:# e.g. 20 for solar but can change
        #elif timestep_temps_final[j] > 4.4 and timestep_H_abundances_final[j] > 0.9 and initial_masses[j] > minimum_wr_mass:# e.g. 20 for solar but can change
            population_choice.append('WR')
            if timestep_H_abundances_final[j] > 0.1:
                for i in range(len(WN_spec_params_reform)):
                    distance_to_spec.append((WN_spec_params_teff[i] - 10**timestep_temps_final[j])**2)
                nearest_spec_ind = np.argmin(distance_to_spec)
                near_spec_temp = WN_spec_params_teff[nearest_spec_ind]
                assigned_spec_teff.append(near_spec_temp)
                near_spec_logl = WN_spec_params_logl[nearest_spec_ind]
                assigned_spec_logl.append(near_spec_logl)
                assigned_spec_logg.append(-1.)
                assigned_spectra.append(WN_reformed_spec_grid[nearest_spec_ind])
                assigned_integrated_spectra.append(WN_integrated_spectra[nearest_spec_ind])
            elif timestep_cnr[j] < 10.:
                for i in range(len(WN_spec_params_reform)):
                    distance_to_spec.append((WN_spec_params_teff[i] - 10**timestep_temps_final[j])**2)
                nearest_spec_ind = np.argmin(distance_to_spec)
                near_spec_temp = WN_spec_params_teff[nearest_spec_ind]
                assigned_spec_teff.append(near_spec_temp)
                near_spec_logl = WN_spec_params_logl[nearest_spec_ind]
                assigned_spec_logl.append(near_spec_logl)
                assigned_spec_logg.append(-1.)
                assigned_spectra.append(WN_reformed_spec_grid[nearest_spec_ind])
                assigned_integrated_spectra.append(WN_integrated_spectra[nearest_spec_ind])
            elif timestep_coher[j] < 0.5:
                for i in range(len(WC_spec_params_reform)):
                    distance_to_spec.append((WC_spec_params_teff[i] - 10**timestep_temps_final[j])**2)
                nearest_spec_ind = np.argmin(distance_to_spec)
                near_spec_temp = WC_spec_params_teff[nearest_spec_ind]
                assigned_spec_teff.append(near_spec_temp)
                near_spec_logl = WC_spec_params_logl[nearest_spec_ind]
                assigned_spec_logl.append(near_spec_logl)
                assigned_spec_logg.append(-1.)
                assigned_spectra.append(WC_reformed_spec_grid[nearest_spec_ind])
                assigned_integrated_spectra.append(WC_integrated_spectra[nearest_spec_ind])
#                vinfs.append(vinf)
            elif timestep_coher[j] < 1.:
                for i in range(len(WC_spec_params_reform)):
                    distance_to_spec.append((WC_spec_params_teff[i] - 10**timestep_temps_final[j])**2)
                nearest_spec_ind = np.argmin(distance_to_spec)
                near_spec_temp = WC_spec_params_teff[nearest_spec_ind]
                assigned_spec_teff.append(near_spec_temp)
                near_spec_logl = WC_spec_params_logl[nearest_spec_ind]
                assigned_spec_logl.append(near_spec_logl)
                assigned_spec_logg.append(-1.)
                assigned_spectra.append(WC_reformed_spec_grid[nearest_spec_ind])
                assigned_integrated_spectra.append(WC_integrated_spectra[nearest_spec_ind])
            elif timestep_coher[j] >= 1.:
                for i in range(len(WC_spec_params_reform)):
                    distance_to_spec.append((WC_spec_params_teff[i] - 10**timestep_temps_final[j])**2)
                nearest_spec_ind = np.argmin(distance_to_spec)
                near_spec_temp = WC_spec_params_teff[nearest_spec_ind]
                assigned_spec_teff.append(near_spec_temp)
                near_spec_logl = WC_spec_params_logl[nearest_spec_ind]
                assigned_spec_logl.append(near_spec_logl)
                assigned_spec_logg.append(-1.)
                assigned_spectra.append(WC_reformed_spec_grid[nearest_spec_ind])
                assigned_integrated_spectra.append(WC_integrated_spectra[nearest_spec_ind])

        elif timestep_temps_final[j] < 3.65 and initial_masses[j] > minimum_wr_mass:
            population_choice.append(No_stars[j])

        elif timestep_teffs_final[j] > 20000. and specsyn_loggs[j] > 2.2 and specsyn_loggs[j] < specsyn_cotests[j]:

            population_choice.append('OB')
            for i in range(len(spec_params_reform)):
                #distance_to_spec.append((spec_params_teff[i] - 10**timestep_temps_final[j])**2 + 1e5*(spec_params_logl[i] - timestep_lums_final[j])**2)
                #distance_to_spec.append((spec_params_teff[i] - 10**timestep_temps_final[j])**2)# + 1e5*(spec_params_logl[i] - timestep_lums_final[j])**2)
                #distance_to_spec.append((spec_params_teff[i] - 10**timestep_temps_final[j])**2 + (spec_params_logl[i] - timestep_lums_final[j])**2)
                #distance_to_spec.append((10**timestep_temps_final[j] - spec_params_teff[i])**2 + 5*(timestep_loggs_final[j] - spec_params_logg[i])**2)
                #distance_to_spec.append((spec_params_teff[i] - 10**timestep_temps_final[j])**2 + (spec_params_logg[i] - timestep_loggs_final[j])**2)
                distance_to_spec.append(abs(spec_params_teff[i] - 10**timestep_temps_final[j]) + abs(spec_params_logg[i] - timestep_loggs_final[j]))
            nearest_spec_ind = np.argmin(distance_to_spec)
            near_spec_temp = spec_params_teff[nearest_spec_ind]
            assigned_spec_teff.append(near_spec_temp)
            near_spec_logl = spec_params_logl[nearest_spec_ind]
            assigned_spec_logl.append(near_spec_logl)
            near_spec_logg = spec_params_logg[nearest_spec_ind]
            assigned_spec_logg.append(near_spec_logg)
            assigned_spectra.append(reformed_spec_grid[nearest_spec_ind])
            assigned_integrated_spectra.append(integrated_spectra[nearest_spec_ind])
        else:
            population_choice.append('lowmass')
            for i in range(len(lowmass_params)):
                distance_to_spec.append((lowmass_teffs[i] - 10**timestep_temps_final[j])**2 + (lowmass_loggs[i] - timestep_loggs_final[j])**2)
            nearest_spec_ind = np.argmin(distance_to_spec)
            near_spec_temp = lowmass_teffs[nearest_spec_ind]
            assigned_spec_teff.append(near_spec_temp)
            near_spec_logl = lm_spec_params_logl[nearest_spec_ind]
            assigned_spec_logl.append(timestep_lums_final[j])
            near_spec_logg = lowmass_loggs[nearest_spec_ind]
            assigned_spec_logg.append(near_spec_logg)
            #assigned_spec_logl.append(near_spec_logl)
            assigned_spectra.append(lowmass_spec[nearest_spec_ind])
            assigned_integrated_spectra.append(lowmass_int_spec[nearest_spec_ind])

    return(assigned_integrated_spectra,assigned_spec_teff,assigned_spectra, population_choice, assigned_spec_logl, assigned_spec_logg)


def specsyn(assigned_integrated_spectra, specsyn_bbfluxes, assigned_spectra, radii, No_stars):
    template_wave = assigned_spectra[0][:, 0]
    assigned_flux_renormed = []
    for i in range(len(assigned_integrated_spectra)):
        xinte = assigned_integrated_spectra[i] / specsyn_bbfluxes[i]
        assigned_flux_renormed.append(assigned_spectra[i][:, 1] / xinte)
    assigned_flux_scaled = np.array([12.566 * radii[i]**2 / 1e20 * assigned_flux_renormed[i] * No_stars[i]
                                     for i in range(len(assigned_flux_renormed))])
    population_flux = assigned_flux_scaled.sum(axis=0)
    return population_flux, []


def specsyn_hires(assigned_integrated_spectra, specsyn_bbfluxes, assigned_spectra, assigned_cont,
                  radii, No_stars, population_nebular, hires_wave_grid):
    population_nebular_resampled = np.interp(hires_wave_grid, population_nebular[:, 0], population_nebular[:, 1])
    assigned_flux        = [assigned_spectra[i][:, 1] for i in range(len(assigned_spectra))]
    assigned_cont_flux   = [assigned_cont[i][:, 1]    for i in range(len(assigned_cont))]
    assigned_flux_scaled = np.array([12.566 * radii[i]**2 / 1e20 * assigned_flux[i]      * No_stars[i] for i in range(len(assigned_flux))])
    assigned_cont_scaled = np.array([12.566 * radii[i]**2 / 1e20 * assigned_cont_flux[i] * No_stars[i] for i in range(len(assigned_cont_flux))])
    population_flux      = assigned_flux_scaled.sum(axis=0)
    population_cont      = assigned_cont_scaled.sum(axis=0)
    population_flux_norm = population_flux / (population_cont + 1.0e-35)
    population_flux      = np.add(population_flux, population_nebular_resampled)
    return population_flux, population_flux_norm, population_nebular_resampled


def compute_ion_flux(wave, freq, flux_wave, flux_freq, limit, code):
    ind_ION      = [i for i in range(len(wave)) if wave[i] <= limit]
    wave_ION     = wave[ind_ION]
    freq_ION     = freq[ind_ION]
    flux_wave_ION= flux_wave[ind_ION]
    flux_freq_ION= flux_freq[ind_ION]
    integral_photons_ION = np.trapezoid(flux_freq_ION, freq_ION)
    integral_flux_ION    = np.trapezoid(flux_wave_ION, wave_ION)
    int_ION = [flux_wave_ION[i] / (6.626e-27 * wave_ION[i]) for i in range(len(wave_ION))]
    Qi      = np.sum(int_ION) * 4 * np.pi * 18**2
    Qilog   = np.log10(Qi)
    if code == 1:
        No_photons_ION = np.log10(integral_photons_ION + 1.0e-30) + 20.
        ion_flux_ION   = np.log10(-1.0*integral_flux_ION + 1.0e-30) + 20.
    if code == 2:
        No_photons_ION = np.log10(-1.0*integral_photons_ION + 1.0e-30) + 20.
        ion_flux_ION   = np.log10(integral_flux_ION + 1.0e-30) + 20.
    return No_photons_ION, ion_flux_ION, Qilog


def ionise(wave, pop_flux, code):
    spectrum_flux      = np.array(pop_flux)
    spectrum_freq      = 2.997925e18 / wave
    spectrum_flux_conv = 3.33564e-19 * spectrum_flux * wave**2 / spectrum_freq / 6.6262e-27
    full_lum           = np.trapezoid(np.flip(spectrum_flux), np.flip(wave))
    bolo_lum = (np.log10(full_lum + 1.0e-30) + 20.) if code == 1 else (np.log10(-1.0*full_lum + 1.0e-30) + 20.)
    ionising_flux_HI   = compute_ion_flux(wave, spectrum_freq, spectrum_flux, spectrum_flux_conv, 912.,  code)
    ionising_flux_HeI  = compute_ion_flux(wave, spectrum_freq, spectrum_flux, spectrum_flux_conv, 504.,  code)
    ionising_flux_HeII = compute_ion_flux(wave, spectrum_freq, spectrum_flux, spectrum_flux_conv, 228.,  code)
    return bolo_lum, ionising_flux_HI, ionising_flux_HeI, ionising_flux_HeII


def continuum(pop_ionising_flux):
    xrange = np.array([10.,912.,913.,1300.,1500.,1800.,2200.,2855.,3331.,3421.,3422.,3642.,3648.,
                       5700.,7000.,8207.,8209.,14583.,14585.,22787.,22789.,32813.,32815.,44680.,44682.,2000000.])
    gamma  = np.array([0.,0.,2.11e-4,5.647,9.35,9.847,10.582,16.101,24.681,26.736,24.883,29.979,6.519,
                       8.773,11.545,13.585,6.333,10.444,7.023,9.361,7.59,9.35,8.32,9.53,8.87,0.])
    cont   = 2.998e18 / xrange / xrange / 2.60e-13 * 1.0e-30 * gamma * 10**(pop_ionising_flux - 30.)
    return np.column_stack((xrange, cont))


def calc_wind(timestep_temps, timestep_lums, timestep_masses, timestep_mdot, timestep_H, timestep_12C,
              timestep_14N, timestep_16O, initial_masses, No_stars, timestep_radii, timestep_cnr,
              minimum_wr_mass, Z):
    timestep_mdot   = np.array(timestep_mdot)
    timestep_masses = np.array(timestep_masses)
    timestep_radii  = np.array(timestep_radii)
    timestep_lums   = np.array(timestep_lums)
    timestep_H      = np.array(timestep_H)
    timestep_12C    = np.array(timestep_12C)
    timestep_14N    = np.array(timestep_14N)
    timestep_16O    = np.array(timestep_16O)
    timestep_coher  = ((timestep_12C/12.) + (timestep_16O/16.)) / (timestep_H/4.)
    timestep_teffs  = 10**timestep_temps
    timestep_vesc   = ((2*6.67e-11*(timestep_masses*2e30)) / (timestep_radii*7e8))**0.5 / 1e3
    vinfs, vinf_vink21, mdot_vink00, mdot_leuven, vinf_xshootu = [], [], [], [], []

    for i in range(len(timestep_temps)):
        if timestep_lums[i] == -20.0:
            vinf = vinf_vink21i = vinf_xshootui = 0.
            mdot_vink00i = mdot_leuveni = timestep_mdot[i]
        elif timestep_temps[i] < 4.4 and timestep_temps[i] > 3.75 and timestep_mdot[i] > -3.5:
            vinf = vinf_vink21i = vinf_xshootui = 200.
            mdot_vink00i = mdot_leuveni = timestep_mdot[i]
        elif timestep_temps[i] < 3.9:
            vinf = vinf_vink21i = vinf_xshootui = 30.
            mdot_vink00i = mdot_leuveni = timestep_mdot[i]
        elif timestep_temps[i] > 4.4 and timestep_H[i] < 0.4 and initial_masses[i] > minimum_wr_mass:
            vinf = 1650. if timestep_H[i] > 0.1 else (1900. if timestep_cnr[i] < 10. else (1800. if timestep_coher[i] < 0.5 else (2800. if timestep_coher[i] < 1. else 3500.)))
            vinf_vink21i = vinf_xshootui = vinf
            mdot_vink00i = mdot_leuveni = timestep_mdot[i]
        elif timestep_temps[i] > 3.9:
            gamma0 = max(1. - (2.7e-5 * 10.**(timestep_lums[i]) / timestep_masses[i]), 1e-10)
            vinf   = 618. * np.sqrt(timestep_masses[i] / (10.**(0.5*timestep_lums[i] - 2.*timestep_temps[i] + 7.52)) * gamma0) * (0.58 + 2.04*(0.5*timestep_lums[i] - 2.*timestep_temps[i] + 7.52))
            timestep_mdot[i] = np.log10(10**timestep_mdot[i])
            mdot_leuveni   = -5.52 + 2.39*np.log10(10**timestep_lums[i]/1e6) - 1.48*np.log10(timestep_masses[i]/45.) + 2.12*np.log10(timestep_teffs[i]/45000.) + (0.75 - 1.87*np.log10(timestep_teffs[i]/40000.))*np.log10(1./1.)
            vinf_xshootui  = 0.092*timestep_teffs[i] - 1040.0
            if timestep_teffs[i] < 25883.:
                vinf_vink21i = 10**(-7.79 - 0.07*timestep_lums[i] + 2.57*timestep_temps[i])
                mdot_vink00i = -6.688 + 2.210*np.log10(10**timestep_lums[i]/1e5) - 1.339*np.log10(timestep_masses[i]/30.) - 1.601*np.log10((vinf_vink21i/timestep_vesc[i])/2.) + 1.07*np.log10(timestep_teffs[i]/20000.)
            else:
                vinf_vink21i = 10**(0.39 - 0.04*timestep_lums[i] + 0.74*timestep_temps[i])
                mdot_vink00i = -6.697 + 2.194*np.log10(10**timestep_lums[i]/1e5) - 1.313*np.log10(timestep_masses[i]/30.) - 1.226*np.log10((vinf_vink21i/timestep_vesc[i])/2.) + 0.933*np.log10(timestep_teffs[i]/40000.) - 10.92*(np.log10(timestep_teffs[i]/40000.))**2
        else:
            vinf = vinf_vink21i = vinf_xshootui = 0.
            mdot_vink00i = mdot_leuveni = timestep_mdot[i]

        vinfs.append(vinf)
        vinf_vink21.append(vinf_vink21i)
        mdot_vink00.append(mdot_vink00i)
        mdot_leuven.append(mdot_leuveni)
        vinf_xshootu.append(vinf_xshootui)

    Z_values = {'MWC': 1.4285, 'MW': 1.0, 'LMC': 0.4285, 'SMC': 0.1428, 'IZw18': 0.0285, 'XMP': 0.00071428, 'Z0': 1e-5}
    Z_value  = Z_values[Z]
    vinfs_Z  = np.array(vinfs) * (Z_value**0.13)
    vinf_xshootu = np.array(vinf_xshootu) * (Z_value**0.22)
    vinf_vink21  = np.array(vinf_vink21)
    mdot_vink00  = np.array(mdot_vink00)
    mdot_leuven  = np.array(mdot_leuven)

    wind_mom_sum         = np.sum(2 * 10.**timestep_mdot * vinfs_Z * 3.155e-5 * No_stars)
    wind_power_sum       = np.sum(10.**timestep_mdot * vinfs_Z**2 * 3.155 * No_stars)
    wind_mom_vink_sum    = np.sum(2 * 10.**mdot_vink00 * vinf_vink21 * 3.155e-5 * No_stars)
    wind_mom_leuven_sum  = np.sum(2 * 10.**mdot_leuven * vinfs_Z * 3.155e-5 * No_stars)
    wind_mom_xshootu_sum = np.sum(2 * 10.**timestep_mdot * vinf_xshootu * 3.155e-5 * No_stars)
    wind_power_xshootu_sum = np.sum(10.**timestep_mdot * vinf_xshootu**2 * 3.155 * No_stars)

    wpower          = np.log10(wind_power_sum + 1e-30) + 35.
    wmom            = np.log10(wind_mom_sum + 1e-30) + 35.
    wmom_vink       = np.log10(wind_mom_vink_sum + 1e-30) + 35.
    wmom_leuven     = np.log10(wind_mom_leuven_sum + 1e-30) + 35.
    wmom_xshootu    = np.log10(wind_mom_xshootu_sum + 1e-30) + 35.
    wpower_xshootu  = np.log10(wind_power_xshootu_sum + 1e-30) + 35.

    return (vinfs_Z, 10.**timestep_mdot * vinfs_Z**2 * 3.155 * No_stars, wpower,
            2 * 10.**timestep_mdot * vinfs_Z * 3.155e-5 * No_stars, wmom, wmom_vink,
            vinf_vink21, mdot_vink00, timestep_vesc, wmom_leuven, wmom_xshootu, wpower_xshootu)


def compute_radii(temps, lums):
    sigma    = 5.670e-8
    lums_sol = 10**np.array(lums) * 3.828e26
    teffs    = 10**temps
    radii    = []
    for i in range(len(teffs)):
        R = 0.1 if temps[i] == 0.0 else (lums_sol[i] / (4 * np.pi * sigma * teffs[i]**4))**0.5 / (7e8)
        radii.append(R)
    return radii


def planck(wave, teff):
    h, c, boltz = 6.626196e-27, 2.997925e10, 1.380622e-16
    cgswave = 1.0e-8 * wave
    c1      = h * c / boltz / teff / cgswave
    c1      = np.where(c1 > 80, 1.0e-30, c1)
    flux    = 2 * np.pi * h * c**2 / cgswave**5 / (np.exp(c1) - 1) * 1e-8
    return np.column_stack((wave, flux))


def f(x, A, B):
    return A*x + B