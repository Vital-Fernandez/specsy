import os
import numpy as np
from scipy import interpolate

# Browns and Seaton FF methodology
h_cgs = 6.626068e-27  # erg s,  cm2 g / s  ==  erg s
c_cgs = 2.99792458e10  # cm / s
c_angs = 2.99792458e18  # ang / s
eV2_erg = 1.602177e-12
pi = 3.141592
masseCGS = 9.1096e-28
e_proton = 4.80320425e-10  # statCoulomb = 1 erg^1/2 cm^1/2 # Eperez definition electronCGS = 1.60217646e-19 * 3.e9 # Coulomb  # 1eV = 1.60217646e-19 Jul
k_cgs = 1.3806503e-16  # erg / K
H0_ion_Energy = 13.6057  # eV
nu_0 = H0_ion_Energy * eV2_erg / h_cgs  # Hz
H_Ryd_Energy = h_cgs * nu_0

# Coefficients for calculating A_2q The total radiative probability 2s -> 1s (s^-1)
alpha_A = 0.88
beta_A = 1.53
gamma_A = 0.8
lambda_2q = 1215.7  # Angstroms
C_A = 202.0  # (s^-1)
A2q = 8.2249  # (s^-1) Transition probability at lambda = 1215.7

# Free Bound constants
Ryd2erg = 2.1798723e-11  # Rydberg to erg   # (s^-1) Transition probability at lambda = 1215.7


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

