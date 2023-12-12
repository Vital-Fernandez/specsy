import logging
import numpy as np

import lime
from ..tools import truncated_gaussian, flux_distribution
from lime import label_decomposition

try:
    import pyneb as pn
    pyneb_check = True
except ImportError:
    pyneb_check = False

_logger = logging.getLogger('SpecSy')


def TOIII_from_TSIII_relation(T_low):
    # From Hagele et al 2006
    return (0.8403 * T_low / 10000.0 + 0.2689) * 10000.0

def TOII_from_TOIII_relation(T_high, n_e):
    # From Epm and Cotini 2009
    return ((1.2 + 0.002*n_e + 4.2/n_e) / (10000.0/T_high + 0.08 + 0.003*n_e + 2.5/n_e)) * 10000.0

def truncated_SII_density_dist(log=None, SII_lines=('S2_6716A', 'S2_6731A'), temp=10000, S2_pyneb=None, flux_dict=None,
                               n_steps=1000):

    '''

    This function computes the electron density from the [SII]6716,6731A doublet. The line label must adhere to the
    LiMe format.

    The user can input a pandas dataframe lines log. This log should adhere to LiMe formatting. Alternatively, the
    user can provide a dictionary with the lines flux distributions. The keys in this dictionary should be the same as
    in the "SII_lines" argument.

    The emissivity calculation is done using PyNeb. The user can provide its own "S2" Atom object. Otherwise, one is created
    with the default PyNeb atomic data.

    The output density distribution is truncated to avoid values outside the physical emissivity ratios.

    :param log: Lines log with the input fluxes. The pandas dataframe must adhere to LiMe formatting
    :type log: pd.DataFrame

    :param SII_lines: Tupple with the label for the [SII] lines. The default values are ('S2_6716A','S2_6731A')
    :type SII_lines: tuple, optional

    :param temp: Temperature for the emissivity calculation in degrees Kelvin. The default value is 10000 K.
    :type temp: float, optional

    :param S2_pyneb: Pyneb Atom, Atom for the S^+ ion.
    :type S2_pyneb: pyneb.Atom, optional

    :param flux_dict: Dictionary with the flux distribution for the [SII] lines.
    :type flux_dict: dict, optional

    :param n_steps: Number of steps in the Monte-Carlo sampling (only if flux_dict is not provided). The default value is 1000.
    :type n_steps: float, optional

    :return: [SII] electron density distribution.
    :rtype: np.array

    '''

    if flux_dict is None:
        flux_dict = flux_distribution(log, 'auto')

    # Compute the densities
    if (SII_lines[0] in flux_dict) and (SII_lines[1] in flux_dict):

        S2_ratio = flux_dict[SII_lines[0]]/flux_dict[SII_lines[1]]

        RSII, RSII_err = np.mean(S2_ratio), np.std(S2_ratio)
        RSII_dist = truncated_gaussian(RSII, RSII_err, n_steps, low_limit=0.28, up_limit=1.42)

        S2 = S2_pyneb if S2_pyneb is not None else pn.Atom('S', 2)
        neSII_dist = S2.getTemDen(RSII_dist, tem=temp, to_eval='L(6716)/L(6731)')

        if np.any(np.isnan(neSII_dist)):
            _logger.warning(f'ne_[SII] distribution contains nan entries')

    else:
        _logger.info('Both [SII] doublet not found in log, the density was not be calculated')
        neSII_dist = None

    return neSII_dist


def ratio_S23(flux_dict, S2_lines=('S2_6716A', 'S2_6731A'), S3_lines=('S3_9068A', 'S3_9530A'),
              norm_lines=('H1_6563A', 'H1_9546A'), H1_pyneb=None, temp=10000, den=100):

    S_23 = None
    if (S2_lines[0] in flux_dict) and (S2_lines[1] in flux_dict):
        if (S3_lines[0] in flux_dict) and (S3_lines[1] in flux_dict):
            if (norm_lines[0] in flux_dict) and (norm_lines[1] in flux_dict):

                H1_pyneb = H1_pyneb if H1_pyneb is not None else pn.RecAtom('H', 1)

                ion_norm1, norm_lines1, latex_norm1 = label_decomposition(norm_lines[0])
                ion_norm2, norm_lines2, latex_norm2 = label_decomposition(norm_lines[1])

                Hbeta_emis = H1_pyneb.getEmissivity(temp, den, wave=4861)
                S2_norm = H1_pyneb.getEmissivity(temp, den, wave=norm_lines1[0])/Hbeta_emis
                S3_norm = H1_pyneb.getEmissivity(temp, den, wave=norm_lines2[0])/Hbeta_emis

                S_2 = (flux_dict[S2_lines[0]]+flux_dict[S2_lines[1]])/flux_dict[norm_lines[0]] * S2_norm
                S_3 = (flux_dict[S3_lines[0]]+flux_dict[S3_lines[1]])/flux_dict[norm_lines[0]] * S3_norm

                S_23 = S_2 + S_3

            else:
                warn_text = f'{norm_lines[0]} ' if norm_lines[0] not in flux_dict else ""
                warn_text += f'{norm_lines[1]}' if norm_lines[1] not in flux_dict else ""
                _logger.info(f'The normalization lines {warn_text} are missing. Please provide others for S23 calculation')
        else:
            _logger.info('[SIII] lines missing. S_23 could not be calculated')
    else:
        _logger.info('[SII] lines missing for S_23 could not be calculated')

    return S_23


def sulfur_diaz_2020(S_23):

    n_steps = S_23.size
    a_dist = np.random.normal(6.636, 0.010, size=n_steps)
    b_dist = np.random.normal(2.202, 0.050, size=n_steps)
    c_dist = np.random.normal(1.060, 0.098, size=n_steps)

    SH = a_dist + b_dist * np.log10(S_23) + c_dist * np.square(np.log10(S_23))

    return SH

def sufur_diaz_2022(lines_log, S2_lines=('S2_6717A', 'S2_6731A'), S3_lines=('S3_9069A', 'S3_9532A'),
                    S2_norm="H1_6563A", S3_norm="H1_6563A", flux_column=f'line_int', n_steps=5000, temp=10000, den=100):


    # Line distribution container
    dist_dict = {}

    # Check the lines are observed and create Monte-Carlo arrays
    for line in list(S2_lines) + list(S3_lines) + [S2_norm] + [S3_norm]:
        int_flux, err_flux = lines_log.loc[line, flux_column], lines_log.loc[line, f'{flux_column}_err']
        dist_dict[line] = np.random.normal(int_flux, err_flux, size=n_steps)

    # Compute the normalization emissivity
    H1 = pn.RecAtom('H', 1)
    line_H_S2 = lime.Line(S2_norm)
    line_H_S3 = lime.Line(S3_norm)

    # Theoretical emissivity
    Hbeta_emis = H1.getEmissivity(temp, den, wave=4861)
    S2_norm_emis = H1.getEmissivity(temp, den, wave=line_H_S2.wavelength[0])
    S3_norm_emis = H1.getEmissivity(temp, den, wave=line_H_S3.wavelength[0])

    # Calculation S_23
    S_2 = (dist_dict[S2_lines[0]] + dist_dict[S2_lines[1]]) / dist_dict[S2_norm]
    S_3 = (dist_dict[S3_lines[0]] + dist_dict[S3_lines[1]]) / dist_dict[S3_norm]
    S23 = S_2 * (S2_norm_emis/Hbeta_emis) + S_3 * (S3_norm_emis/Hbeta_emis)

    # Calculation abundance coefficients
    k1 = np.random.normal(6.636, 0.011, n_steps)
    k2 = np.random.normal(2.202, 0.050, n_steps)
    k3 = np.random.normal(1.060, 0.098, n_steps)

    # Abundance calculation
    S_H = k1 + k2 * np.log10(S23) + k3 * np.power(np.log10(S23), 2)

    return np.mean(S_H), np.std(S_H)