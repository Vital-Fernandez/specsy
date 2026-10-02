from logging import getLogger

_logger = getLogger('SpecSy')

import os
import numpy as np
import pandas as pd
from scipy import interpolate

from specsy.io import SpecSyError

from astropy import constants
from lime.tools import au


# Physical constants
h_cgs = constants.h.cgs.value  # erg s
c_cgs = constants.c.cgs.value  # cm / s
c_angs = constants.c.to(au.AA / au.s).value  # ang / s
eV2_erg = au.eV.to(au.erg)
masseCGS = constants.m_e.cgs.value  # g
e_proton = constants.e.esu.value  # statCoulomb = 1 erg^1/2 cm^1/2
k_cgs = constants.k_B.cgs.value  # erg / K
_rydberg = constants.h * constants.c * constants.Ryd  # Rydberg energy
H0_ion_Energy = _rydberg.to(au.eV).value  # eV
Ryd2erg = _rydberg.to(au.erg).value  # Rydberg to erg

# Browns and Seaton FF methodology
pi = np.pi
nu_0 = H0_ion_Energy * eV2_erg / h_cgs  # Hz
H_Ryd_Energy = h_cgs * nu_0

# Coefficients for calculating A_2q The total radiative probability 2s -> 1s (s^-1)
alpha_A = 0.88
beta_A = 1.53
gamma_A = 0.8
lambda_2q = 1215.7  # Angstroms
C_A = 202.0  # (s^-1)
A2q = 8.2249  # (s^-1) Transition probability at lambda = 1215.7

# ----------------------------------------------------------------------------------------------------------------------
# Starburst99 nebular continuum constants. Deliberately NOT read from astropy: they are the values used by
# sesamme.models.nebular_continuum, kept identical for parity (e.g. LSUN_ERG = 3.83e33 vs astropy's 3.828e33)
# ----------------------------------------------------------------------------------------------------------------------
C_ANG = 2.998e18        # speed of light [A s-1]
ALPHA_B = 2.6e-13       # Case B H recombination coefficient at 1e4 K [cm3 s-1]
LSUN_ERG = 3.83e33      # same constant used for the SSP erg/s/A -> Lsun/A conversion

# HI (free-free, bound-free, 2-photon) + HeI (He/H = 0.1) emission coefficients, Aller (1984) and Ferland (1980),
# Case B, T = 1e4 K, f_esc = 0 [erg cm3 s-1 Hz-1] at NEB_GAMMA_WAVE [A]
NEB_GAMMA_WAVE = np.array([912., 913., 1300., 1500., 1800., 2200., 2855., 3331., 3421., 3422., 3642., 3648., 5700.,
                           7000., 8207., 8209., 14583., 14585., 22787., 22789., 32813., 32815., 44680., 44682.])
NEB_GAMMA_COEF = np.array([0., 2.11e-4, 5.647, 9.35, 9.847, 10.582, 16.101, 24.681, 26.736, 24.883, 29.979, 6.519,
                           8.773, 11.545, 13.585, 6.333, 10.444, 7.023, 9.361, 7.59, 9.35, 8.32, 9.53, 8.87]) * 1e-40


def importErcolanoTables(file_address):

    """
    This function imports the atomic data from Ercolano et al.(2006) to compute the Free-Bound nebular continuum
    :param file_address:
    :return: dictionary

    """

    dict_ion = {}

    # Reading the text files
    with open(file_address, 'r') as f:

        a = f.readlines()

        dict_ion['nTe'] = int(str.split(a[0])[0])  # number of Te columns
        dict_ion['nEner'] = int(str.split(a[0])[1])  # number of energy points rows
        dict_ion['skip'] = int(1 + np.ceil(dict_ion['nTe'] / 8.))  # 8 es el numero de valores de Te por fila.
        dict_ion['temps'] = np.zeros(dict_ion['nTe'])
        dict_ion['lines'] = a

        # Storing temperature range
        for i in range(1, dict_ion['skip']):
            tt = str.split(a[i])
            for j in range(0, len(tt)):
                dict_ion['temps'][8 * (i - 1) + j] = tt[j]

        # Storing gamma_cross grids
        dict_ion['matrix'] = np.loadtxt(file_address, skiprows=dict_ion['skip'])

    # Get wavelengths corresponding to table threshold and remove zero entries
    wave_thres = dict_ion['matrix'][:, 0] * dict_ion['matrix'][:, 1]
    idx_zero = (wave_thres == 0)
    dict_ion['wave_thres'] = wave_thres[~idx_zero]

    return dict_ion


class NebularContinua:

    def __init__(self, biblio_folder=None):

        # Use default folder if data if no folder is declared
        if biblio_folder is None:
            _dir_path = os.path.dirname(os.path.realpath(__file__))
            biblio_folder = os.path.abspath(os.path.join(_dir_path, os.path.join(os.pardir, 'resources')))

        # Load files
        self.HI_fb_dict = importErcolanoTables(os.path.join(biblio_folder, 'HI_t3_elec.ascii'))
        self.HeI_fb_dict = importErcolanoTables(os.path.join(biblio_folder, 'HeI_t5_elec.ascii'))
        self.HeII_fb_dict = importErcolanoTables(os.path.join(biblio_folder, 'HeII_t4_elec.ascii'))

        return

    def flux_spectrum_backup(self, wave_rest, Te, Halpha_Flux, He1_abund, He2_abund, cHbeta=None, flambda=None):

        neb_gamma = self.gamma_spectrum(wave_rest, Te, He1_abund, He2_abund)

        neb_flux = self.zanstra_calibration(wave_rest, Te, Halpha_Flux, neb_gamma)

        # Apply nebular  corruction if available
        if (cHbeta is not None) and (flambda is not None):
            return neb_flux * np.power(10, -1 * np.flambda_neb * cHbeta)

        else:
            return neb_flux

    def flux_spectrum(self, wave, Te, Halpha_Flux, He1_abund=0.1, He2_abund=0.001):

        neb_gamma = self.gamma_spectrum(wave, Te, He1_abund, He2_abund)

        return self.zanstra_calibration(wave, Te, Halpha_Flux, neb_gamma)

    def components(self, wave, Te, He1_abund=0.1, He2_abund=0.001):

        H_He_frac = 1 + He1_abund * 4 + He2_abund * 4

        # Bound bound continuum
        gamma_2q = self.boundbound_gamma(wave, Te)

        # Free-Free continuum
        gamma_ff = self.freefree_gamma(wave, Te, Z_ion=1.0)

        # Free-Bound continuum
        gamma_fb_HI = self.freebound_gamma(wave, Te, self.HI_fb_dict)

        return gamma_2q, gamma_ff, gamma_fb_HI

    def gamma_spectrum(self, wave, Te, HeII_HII, HeIII_HII):

        H_He_frac = 1 + HeII_HII * 4 + HeIII_HII * 4

        # Bound bound continuum
        gamma_2q = self.boundbound_gamma(wave, Te)

        # Free-Free continuum
        gamma_ff = H_He_frac * self.freefree_gamma(wave, Te, Z_ion=1.0)

        # Free-Bound continuum
        gamma_fb_HI = self.freebound_gamma(wave, Te, self.HI_fb_dict)
        gamma_fb_HeI = self.freebound_gamma(wave, Te, self.HeI_fb_dict)
        gamma_fb_HeII = self.freebound_gamma(wave, Te, self.HeII_fb_dict)
        gamma_fb = gamma_fb_HI + HeII_HII * gamma_fb_HeI + HeIII_HII * gamma_fb_HeII

        return gamma_2q + gamma_ff + gamma_fb

    def boundbound_gamma(self, wave, Te):

        # Prepare arrays
        idx_limit = (wave > lambda_2q)
        gamma_array = np.zeros(wave.size)

        # Params
        q2 = 5.92e-4 - 6.1e-9 * Te  # (cm^3 s^-1) Collisional transition rate coefficient for protons and electrons
        alpha_eff_2q = 6.5346e-11 * np.power(Te, -0.72315) # (cm^3 s^-1) Effective Recombination coefficient

        nu_array = c_angs / wave[idx_limit]
        nu_limit = c_angs / lambda_2q

        y = nu_array / nu_limit

        A_y = C_A * (y * (1 - y) * (1 - (4 * y * (1 - y)) ** gamma_A) + alpha_A *
                    (y * (1 - y)) ** beta_A * (4 * y * (1 - y)) ** gamma_A)

        g_nu1 = h_cgs * nu_array / nu_limit / A2q * A_y
        g_nu2 = alpha_eff_2q * g_nu1 / (1 + q2 / A2q)

        gamma_array[idx_limit] = g_nu2[:]

        return gamma_array

    def freefree_gamma(self, wave, Te, Z_ion=1.0):

        cte_A = (32 * (Z_ion ** 2) * (e_proton ** 4) * h_cgs) / (
                    3 * (masseCGS ** 2) * (c_cgs ** 3))

        cte_B = ((np.pi * H_Ryd_Energy / (3 * k_cgs * Te)) ** 0.5)

        cte_Total = cte_A * cte_B

        nu_array = c_angs / wave

        gamma_Comp1 = np.exp(((-1 * h_cgs * nu_array) / (k_cgs * Te)))
        gamma_Comp2 = (h_cgs * nu_array) / ((Z_ion ** 2) * e_proton * 13.6057)
        gamma_Comp3 = k_cgs * Te / (h_cgs * nu_array)

        gff = 1 + 0.1728 * np.power(gamma_Comp2, 0.33333) * (1 + 2 * gamma_Comp3) - 0.0496 * np.power(gamma_Comp2, 0.66667) * (
                          1 + 0.66667 * gamma_Comp3 + 1.33333 * np.power(gamma_Comp3, 2))

        gamma_array = cte_Total * gamma_Comp1 * gff

        return gamma_array

    def freebound_gamma(self, wave, Te, data_dict):

        wave_ryd = (h_cgs * c_angs) / (Ryd2erg * wave)

        # Temperature entry
        t4 = Te / 10000.0
        logTe = np.log10(Te)

        # Interpolating the grid for the right wavelength
        thres_idxbin = np.digitize(wave_ryd, data_dict['wave_thres'])
        interpol_grid = interpolate.interp2d(data_dict['temps'], data_dict['matrix'][:, 1], data_dict['matrix'][:, 2:], kind='linear')

        ener_low = data_dict['wave_thres'][thres_idxbin - 1] # WARNING: This one could be an issue for wavelength range table limits

        # Interpolate table for the right temperature
        gamma_inter_Te = interpol_grid(logTe, data_dict['matrix'][:, 1])[:, 0]
        gamma_inter_Te_Ryd = np.interp(wave_ryd, data_dict['matrix'][:, 1], gamma_inter_Te)

        Gamma_fb_f = gamma_inter_Te_Ryd * 1e-40 * np.power(t4, -1.5) * np.exp(-15.7887 * (wave_ryd - ener_low) / t4)

        return Gamma_fb_f

    def zanstra_calibration(self, wave, Te, flux_Emline, gNeb_cont_nu, lambda_EmLine=6562.819):

        # Zanstra like calibration for the continuum
        t4 = Te / 10000.0

        # Pequignot et al. 1991
        # alfa_eff_alpha = 2.708e-13 * t4 ** -0.648 / (1 + 1.315 * t4 ** 0.523)
        alfa_eff_beta = 0.668e-13 * t4**-0.507 / (1 + 1.221*t4**0.653)

        fNeb_cont_lambda = gNeb_cont_nu * lambda_EmLine * flux_Emline / (alfa_eff_beta * h_cgs * wave * wave)

        return fNeb_cont_lambda

    def zanstra_calibration_tt(self, wave, Te, flux_Emline, gNeb_cont_nu, lambda_EmLine=6562.819):

        # Zanstra like calibration for the continuum
        t4 = Te / 10000.0

        # Pequignot et al. 1991
        alfa_eff_alpha = 2.708e-13 * np.power(t4, -0.648) / (1 + 1.315 * np.power(t4, 0.523))
        fNeb_cont_lambda = gNeb_cont_nu * lambda_EmLine * flux_Emline / (alfa_eff_alpha * h_cgs * wave * wave)

        return fNeb_cont_lambda


# ----------------------------------------------------------------------------------------------------------------------
# SSP sampling constants (tabulated Aller 1984 / Ferland 1980 coefficients, scaled by Q)
# ----------------------------------------------------------------------------------------------------------------------

def load_ionization_frame(file_name):

    """
    Loads in the table of ionizing photon outputs per SSP associated with a model cube.

    Parameters
    ----------
    file_name : str
        File path and name

    Returns
    -------
    ion_table : pandas.DataFrame
        DataFrame containing ionizing fluxes per SSP
    """

    return pd.read_csv(file_name, sep=r'\s+', header=0, comment='#')


def _met_key_value(key):
    """Linear Z from an ionization-table metallicity key: 'Z004' -> 0.004, 'Zem5' -> 1e-5 (numbers pass through)."""
    try:
        return float(key)
    except ValueError:
        digits = str(key).strip().lstrip('Zz')
        return 10.0 ** -int(digits[2:]) if digits.startswith('em') else int(digits) / 1000.0


def ion_table_grid(ion_table, z_col='Z'):
    """
    Axes and values of a SESAMME ionization table (one row per metallicity key in `z_col`, one column per log(age),
    e.g. '6.0'): log10 Z (ascending), log(age) (ascending) and the log10 Q matrix with shape (n_Z, n_age).
    """
    age_cols = {}
    for col in ion_table.columns:
        if col != z_col:
            try:
                age_cols[float(col)] = col
            except ValueError:
                pass

    tab_age = np.array(sorted(age_cols))
    tab_z = np.log10(ion_table[z_col].map(_met_key_value).to_numpy(dtype=float))
    q_mat = ion_table[[age_cols[a] for a in tab_age]].to_numpy(dtype=float)

    order = np.argsort(tab_z)
    return tab_z[order], tab_age, q_mat[order]


def ion_table_log_q(ion_table, age_node, z_node, z_col='Z', atol=1e-3):
    """
    log10 Q of each (log age, log10 Z) node from a SESAMME ionization table. Raises if a node has no entry.
    """
    tab_z, tab_age, q_mat = ion_table_grid(ion_table, z_col)

    age_node, z_node = np.asarray(age_node, dtype=float), np.asarray(z_node, dtype=float)
    match_z = np.isclose(z_node[:, None], tab_z[None, :], atol=atol, rtol=0)
    match_age = np.isclose(age_node[:, None], tab_age[None, :], atol=atol, rtol=0)

    missing = ~(match_z.any(axis=1) & match_age.any(axis=1))
    if missing.any():
        pairs = [f'({t:.2f}, {10 ** z:.5f})' for t, z in zip(age_node[missing][:5], z_node[missing][:5])]
        raise SpecSyError(f'{int(missing.sum())} of the {age_node.size} (log(age), Z) nodes have no entry in the '
                          f'ionization table (atol={atol}), e.g. {", ".join(pairs)}. Table covers log(age) '
                          f'{tab_age.min():.2f}-{tab_age.max():.2f} and Z keys {ion_table[z_col].tolist()}.')

    return q_mat[match_z.argmax(axis=1), match_age.argmax(axis=1)]


def nebular_continuum(wave, log_q, gamma_wave=NEB_GAMMA_WAVE, gamma_coef=NEB_GAMMA_COEF, alpha_b=ALPHA_B,
                      lsun=LSUN_ERG):
    """
    Case B nebular continuum (H I + He I free-free, bound-free and 2-photon; He/H = 0.1) for each SSP node,
    following the Starburst99 CONTINUUM recipe: L_nu = gamma_nu * Q / alpha_B, L_lambda = L_nu * c / lambda^2.

    gamma_coef: emission coefficients in erg cm3 s-1 Hz-1 at gamma_wave [A]. Default: the SESAMME/SB99 table.
    log_q: log10 H-ionizing photon rate of each node, with the same mass normalization as the SSP fluxes.
    Returns an (n_node, n_pixel) array in Lsun/A.
    """
    wave = np.asarray(wave, dtype=float)
    log_q = np.atleast_1d(np.asarray(log_q, dtype=float))
    order = np.argsort(gamma_wave)
    g_wave, g_coef = np.asarray(gamma_wave, float)[order], np.asarray(gamma_coef, float)[order]

    out = (wave < g_wave[0]) | (wave > g_wave[-1])
    if out.any():
        _logger.warning(f'{int(out.sum())} of the {wave.size} pixels fall outside the gamma table '
                        f'({g_wave[0]:.1f}-{g_wave[-1]:.1f} A). Their nebular continuum is set to zero.')

    bad_q = ~np.isfinite(log_q)
    if bad_q.any():
        _logger.warning(f'{int(bad_q.sum())} of the {log_q.size} nodes have a non-finite log(Q). '
                        f'Their nebular continuum is set to zero.')

    # L_lambda per ionizing photon at the table nodes, then interpolated (the SESAMME order of operations)
    node_lum = g_coef / alpha_b * C_ANG / g_wave ** 2 / lsun             # Lsun/A per (photon s-1)
    lum_per_photon = np.interp(wave, g_wave, node_lum, left=0.0, right=0.0)
    q = np.where(bad_q, 0.0, np.power(10.0, np.where(bad_q, 0.0, log_q)))

    return q[:, None] * lum_per_photon[None, :]


def nebular_grid_sesamme(wave, age_node, z_node, ion_table):
    """Same grid computed with sesamme.models.nebular_continuum (E(B-V) = 0, log(A) = 0), for cross-checking."""
    from astropy.table import Table
    from sesamme import models as ses_models
    if isinstance(ion_table, pd.DataFrame):
        ion_table = Table.from_pandas(ion_table)
    return np.stack([np.asarray(ses_models.nebular_continuum(wave, [t, z, 0.0, 0.0], ion_table), dtype=float)
                     for t, z in zip(age_node, z_node)])



# def importErcolanoTables(file_address):
#
#     """
#     This function imports the atomic data from Ercolano et al.(2006) to compute the Free-Bound nebular continuum
#     :param file_address:
#     :return: dictionary
#
#     """
#
#     dict_ion = {}
#
#     # Reading the text files
#     with open(file_address, 'r') as f:
#
#         a = f.readlines()
#
#         dict_ion['nTe'] = int(str.split(a[0])[0])  # number of Te columns
#         dict_ion['nEner'] = int(str.split(a[0])[1])  # number of energy points rows
#         dict_ion['skip'] = int(1 + np.ceil(dict_ion['nTe'] / 8.))  # 8 es el numero de valores de Te por fila.
#         dict_ion['temps'] = np.zeros(dict_ion['nTe'])
#         dict_ion['lines'] = a
#
#         # Storing temperature range
#         for i in range(1, dict_ion['skip']):
#             tt = str.split(a[i])
#             for j in range(0, len(tt)):
#                 dict_ion['temps'][8 * (i - 1) + j] = tt[j]
#
#         # Storing gamma_cross grids
#         dict_ion['matrix'] = np.loadtxt(file_address, skiprows=dict_ion['skip'])
#
#     # Get wavelengths corresponding to table threshold and remove zero entries
#     wave_thres = dict_ion['matrix'][:, 0] * dict_ion['matrix'][:, 1]
#     idx_zero = (wave_thres == 0)
#     dict_ion['wave_thres'] = wave_thres[~idx_zero]
#
#     return dict_ion
#
#
# class NebularContinua:
#
#     def __init__(self, biblio_folder=None):
#
#         # Use default folder if data if no folder is declared
#         if biblio_folder is None:
#             _dir_path = os.path.dirname(os.path.realpath(__file__))
#             biblio_folder = os.path.abspath(os.path.join(_dir_path, os.path.join(os.pardir, 'resources')))
#
#         # Load files
#         self.HI_fb_dict = importErcolanoTables(os.path.join(biblio_folder, 'HI_t3_elec.ascii'))
#         self.HeI_fb_dict = importErcolanoTables(os.path.join(biblio_folder, 'HeI_t5_elec.ascii'))
#         self.HeII_fb_dict = importErcolanoTables(os.path.join(biblio_folder, 'HeII_t4_elec.ascii'))
#
#         return
#
#     def flux_spectrum_backup(self, wave_rest, Te, Halpha_Flux, He1_abund, He2_abund, cHbeta=None, flambda=None):
#
#         neb_gamma = self.gamma_spectrum(wave_rest, Te, He1_abund, He2_abund)
#
#         neb_flux = self.zanstra_calibration(wave_rest, Te, Halpha_Flux, neb_gamma)
#
#         # Apply nebular  corruction if available
#         if (cHbeta is not None) and (flambda is not None):
#             return neb_flux * np.power(10, -1 * np.flambda_neb * cHbeta)
#
#         else:
#             return neb_flux
#
#     def flux_spectrum(self, wave, Te, Halpha_Flux, He1_abund=0.1, He2_abund=0.001):
#
#         neb_gamma = self.gamma_spectrum(wave, Te, He1_abund, He2_abund)
#
#         return self.zanstra_calibration(wave, Te, Halpha_Flux, neb_gamma)
#
#     def components(self, wave, Te, He1_abund=0.1, He2_abund=0.001):
#
#         H_He_frac = 1 + He1_abund * 4 + He2_abund * 4
#
#         # Bound bound continuum
#         gamma_2q = self.boundbound_gamma(wave, Te)
#
#         # Free-Free continuum
#         gamma_ff = self.freefree_gamma(wave, Te, Z_ion=1.0)
#
#         # Free-Bound continuum
#         gamma_fb_HI = self.freebound_gamma(wave, Te, self.HI_fb_dict)
#
#         return gamma_2q, gamma_ff, gamma_fb_HI
#
#     def gamma_spectrum(self, wave, Te, HeII_HII, HeIII_HII):
#
#         H_He_frac = 1 + HeII_HII * 4 + HeIII_HII * 4
#
#         # Bound bound continuum
#         gamma_2q = self.boundbound_gamma(wave, Te)
#
#         # Free-Free continuum
#         gamma_ff = H_He_frac * self.freefree_gamma(wave, Te, Z_ion=1.0)
#
#         # Free-Bound continuum
#         gamma_fb_HI = self.freebound_gamma(wave, Te, self.HI_fb_dict)
#         gamma_fb_HeI = self.freebound_gamma(wave, Te, self.HeI_fb_dict)
#         gamma_fb_HeII = self.freebound_gamma(wave, Te, self.HeII_fb_dict)
#         gamma_fb = gamma_fb_HI + HeII_HII * gamma_fb_HeI + HeIII_HII * gamma_fb_HeII
#
#         return gamma_2q + gamma_ff + gamma_fb
#
#     def boundbound_gamma(self, wave, Te):
#
#         # Prepare arrays
#         idx_limit = (wave > lambda_2q)
#         gamma_array = np.zeros(wave.size)
#
#         # Params
#         q2 = 5.92e-4 - 6.1e-9 * Te  # (cm^3 s^-1) Collisional transition rate coefficient for protons and electrons
#         alpha_eff_2q = 6.5346e-11 * np.power(Te, -0.72315) # (cm^3 s^-1) Effective Recombination coefficient
#
#         nu_array = c_angs / wave[idx_limit]
#         nu_limit = c_angs / lambda_2q
#
#         y = nu_array / nu_limit
#
#         A_y = C_A * (y * (1 - y) * (1 - (4 * y * (1 - y)) ** gamma_A) + alpha_A *
#                     (y * (1 - y)) ** beta_A * (4 * y * (1 - y)) ** gamma_A)
#
#         g_nu1 = h_cgs * nu_array / nu_limit / A2q * A_y
#         g_nu2 = alpha_eff_2q * g_nu1 / (1 + q2 / A2q)
#
#         gamma_array[idx_limit] = g_nu2[:]
#
#         return gamma_array
#
#     def freefree_gamma(self, wave, Te, Z_ion=1.0):
#
#         cte_A = (32 * (Z_ion ** 2) * (e_proton ** 4) * h_cgs) / (
#                     3 * (masseCGS ** 2) * (c_cgs ** 3))
#
#         cte_B = ((np.pi * H_Ryd_Energy / (3 * k_cgs * Te)) ** 0.5)
#
#         cte_Total = cte_A * cte_B
#
#         nu_array = c_angs / wave
#
#         gamma_Comp1 = np.exp(((-1 * h_cgs * nu_array) / (k_cgs * Te)))
#         gamma_Comp2 = (h_cgs * nu_array) / ((Z_ion ** 2) * e_proton * 13.6057)
#         gamma_Comp3 = k_cgs * Te / (h_cgs * nu_array)
#
#         gff = 1 + 0.1728 * np.power(gamma_Comp2, 0.33333) * (1 + 2 * gamma_Comp3) - 0.0496 * np.power(gamma_Comp2, 0.66667) * (
#                           1 + 0.66667 * gamma_Comp3 + 1.33333 * np.power(gamma_Comp3, 2))
#
#         gamma_array = cte_Total * gamma_Comp1 * gff
#
#         return gamma_array
#
#     def freebound_gamma(self, wave, Te, data_dict):
#
#         wave_ryd = (h_cgs * c_angs) / (Ryd2erg * wave)
#
#         # Temperature entry
#         t4 = Te / 10000.0
#         logTe = np.log10(Te)
#
#         # Interpolating the grid for the right wavelength
#         thres_idxbin = np.digitize(wave_ryd, data_dict['wave_thres'])
#         interpol_grid = interpolate.interp2d(data_dict['temps'], data_dict['matrix'][:, 1], data_dict['matrix'][:, 2:], kind='linear')
#
#         ener_low = data_dict['wave_thres'][thres_idxbin - 1] # WARNING: This one could be an issue for wavelength range table limits
#
#         # Interpolate table for the right temperature
#         gamma_inter_Te = interpol_grid(logTe, data_dict['matrix'][:, 1])[:, 0]
#         gamma_inter_Te_Ryd = np.interp(wave_ryd, data_dict['matrix'][:, 1], gamma_inter_Te)
#
#         Gamma_fb_f = gamma_inter_Te_Ryd * 1e-40 * np.power(t4, -1.5) * np.exp(-15.7887 * (wave_ryd - ener_low) / t4)
#
#         return Gamma_fb_f
#
#     def zanstra_calibration(self, wave, Te, flux_Emline, gNeb_cont_nu, lambda_EmLine=6562.819):
#
#         # Zanstra like calibration for the continuum
#         t4 = Te / 10000.0
#
#         # Pequignot et al. 1991
#         # alfa_eff_alpha = 2.708e-13 * t4 ** -0.648 / (1 + 1.315 * t4 ** 0.523)
#         alfa_eff_beta = 0.668e-13 * t4**-0.507 / (1 + 1.221*t4**0.653)
#
#         fNeb_cont_lambda = gNeb_cont_nu * lambda_EmLine * flux_Emline / (alfa_eff_beta * h_cgs * wave * wave)
#
#         return fNeb_cont_lambda
#
#     def zanstra_calibration_tt(self, wave, Te, flux_Emline, gNeb_cont_nu, lambda_EmLine=6562.819):
#
#         # Zanstra like calibration for the continuum
#         t4 = Te / 10000.0
#
#         # Pequignot et al. 1991
#         alfa_eff_alpha = 2.708e-13 * np.power(t4, -0.648) / (1 + 1.315 * np.power(t4, 0.523))
#         fNeb_cont_lambda = gNeb_cont_nu * lambda_EmLine * flux_Emline / (alfa_eff_alpha * h_cgs * wave * wave)
#
#         return fNeb_cont_lambda

