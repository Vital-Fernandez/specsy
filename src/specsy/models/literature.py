import numpy as np


# From Hagele et al 2006
def TSIII_from_TOIII_relation(T_high):
    return (1.19 * T_high / 10000.0 - 0.32) * 10000.0


# From Hagele et al 2006
def TOIII_from_TSIII_relation(T_low):
    return (0.8403 * T_low / 10000.0 + 0.2689) * 10000.0


# From Epm and Cotini 2009
def TOII_from_TOIII_relation(T_high, n_e):
    return ((1.2 + 0.002*n_e + 4.2/n_e) / (10000.0/T_high + 0.08 + 0.003*n_e + 2.5/n_e)) * 10000.0


def sulfur_diaz_2020(S_23):

    n_steps = S_23.size
    a_dist = np.random.normal(6.636, 0.010, size=n_steps)
    b_dist = np.random.normal(2.202, 0.050, size=n_steps)
    c_dist = np.random.normal(1.060, 0.098, size=n_steps)

    SH = a_dist + b_dist * np.log10(S_23) + c_dist * np.square(np.log10(S_23))

    return SH


_TEM_FUNC_DICT = {'TSIII_Hagele2006': TSIII_from_TOIII_relation,
                  'TOIII_Hagele2006': TOIII_from_TSIII_relation,
                  'TOII_Epm2009': TOII_from_TOIII_relation}

_DEN_FUNC_DICT = {}


"""
Pagel-Stasinska 1990
t_OII  = 2 / (t_OIII^-1 + 0.8)

Perez-montero_2009
tNII = 1.85 / (t_OIII^-1 + 0.72)

Campbell 1986
TNII = 0.7 *T_OIII + 3000

Cordoba_2020_P.5- (It has uncertainty)
tNII^-1 = 0.54*tNIII^-1 + 0.52
tOIII^-1 = 1.04*tNII^-1 + 0.05

Cordoba_2020_P.5+ (It has uncertainty)
tNII^-1 = 0.61*tNIII^-1 + 0.36
tOIII^-1 = 1.00*tNII^-1 + 0.03

Perez-Montero2003
t_SII = 0.71*t_OII + 0.12

Perez-Montero2005
t_SIII = 1.05*t_OIII + 0.08

Hagele_2006 (It has uncertainty)
t_SIII = 1.19*t_OIII - 0.32

Garnet_1992
tOII = 0.7 * t_OIII + 0.3

Izotov_1997
tOII = 0.243 + 1.031*t_OIII - 0.184*t_OIII^2

Pilyugin2006 (it has uncertainty)
tOII = 0.729*tOIII + 0.257
"""