# import numpy as np
#

# # From Hagele et al 2006
# def TSIII_from_TOIII_relation(T_high):
#     return (1.19 * T_high / 10000.0 - 0.32) * 10000.0
#
#
# # From Hagele et al 2006
# def TOIII_from_TSIII_relation(T_low):
#     return (0.8403 * T_low / 10000.0 + 0.2689) * 10000.0
#
#
# # From Epm and Cotini 2009
# def TOII_from_TOIII_relation(T_high, n_e):
#     return ((1.2 + 0.002*n_e + 4.2/n_e) / (10000.0/T_high + 0.08 + 0.003*n_e + 2.5/n_e)) * 10000.0
#
#
# def sulfur_diaz_2020(S_23):
#
#     n_steps = S_23.size
#     a_dist = np.random.normal(6.636, 0.010, size=n_steps)
#     b_dist = np.random.normal(2.202, 0.050, size=n_steps)
#     c_dist = np.random.normal(1.060, 0.098, size=n_steps)
#
#     SH = a_dist + b_dist * np.log10(S_23) + c_dist * np.square(np.log10(S_23))
#
#     return SH

def TOII_TOIII_PagelStasinska1990(TOIII):
    return 20000 / (TOIII**-1 + 0.8 / 10000)

def TNII_TOIII_PerezMontero2009(TOIII):
    return 18500 / (TOIII**-1 + 0.72 / 10000)

def TNII_TOIII_Campbell1986(TOIII):
    return 0.7 * TOIII + 3000

def TNII_TOIII_Cordova2020_Pminus(TNIII):
    return 1 / (0.54 * TNIII**-1 + 0.52 / 10000)

def TOIII_TNII_Cordova2020_Pminus(TNII):
    return 1 / (1.04 * TNII**-1 + 0.05 / 10000)

def TNII_TOIII_Cordova2020_Pplus(TNIII):
    return 1 / (0.61 * TNIII**-1 + 0.36 / 10000)

def TOIII_TNII_Cordova2020_Pplus(TNII):
    return 1 / (1.00 * TNII**-1 + 0.03 / 10000)

def TSII_TOII_PerezMontero2003(TOII):
    return 0.71 * TOII + 1200

def TSIII_TOIII_PerezMontero2005(TOIII):
    return 1.05 * TOIII + 800

def TSIII_TOIII_Hagele2006(TOIII):
    return 1.19 * TOIII - 3200

def TOII_TOIII_Garnett1992(TOIII):
    return 0.7 * TOIII + 3000

def TOII_TOIII_Izotov1997(TOIII):
    return (0.243 + 1.031 * (TOIII / 10000) - 0.184 * (TOIII / 10000)**2) * 10000

def TOII_TOIII_Pilyugin2006(TOIII):
    return 0.729 * TOIII + 2570


TEM_FUNC_DICT = {'TOII_TOIII_PagelStasinska1990':    TOII_TOIII_PagelStasinska1990,
                 'TNII_TOIII_PerezMontero2009':      TNII_TOIII_PerezMontero2009,
                 'TNII_TOIII_Campbell1986':          TNII_TOIII_Campbell1986,
                 'TNII_TNIII_Cordova2020_Pminus':    TNII_TOIII_Cordova2020_Pminus,
                 'TOIII_TNII_Cordova2020_Pminus':    TOIII_TNII_Cordova2020_Pminus,
                 'TNII_TNIII_Cordova2020_Pplus':     TNII_TOIII_Cordova2020_Pplus,
                 'TOIII_TNII_Cordova2020_Pplus':     TOIII_TNII_Cordova2020_Pplus,
                 'TSII_TOII_PerezMontero2003':       TSII_TOII_PerezMontero2003,
                 'TSIII_TOIII_PerezMontero2005':     TSIII_TOIII_PerezMontero2005,
                 'TSIII_TOIII_Hagele2006':           TSIII_TOIII_Hagele2006,
                 'TOII_TOIII_Garnett1992':           TOII_TOIII_Garnett1992,
                 'TOII_TOIII_Izotov1997':            TOII_TOIII_Izotov1997,
                 'TOII_TOIII_Pilyugin2006':          TOII_TOIII_Pilyugin2006}

DEN_FUNC_DICT = {}


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