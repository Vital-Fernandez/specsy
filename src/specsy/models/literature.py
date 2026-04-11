# From Hagele et al 2006
def TSIII_from_TOIII_relation(T_high):
    return (1.19 * T_high / 10000.0 - 0.32) * 10000.0


# From Hagele et al 2006
def TOIII_from_TSIII_relation(T_low):
    return (0.8403 * T_low / 10000.0 + 0.2689) * 10000.0


# From Epm and Cotini 2009
def TOII_from_TOIII_relation(T_high, n_e):
    return ((1.2 + 0.002*n_e + 4.2/n_e) / (10000.0/T_high + 0.08 + 0.003*n_e + 2.5/n_e)) * 10000.0


_TEM_FUNC_DICT = {'TSIII_Hagele2006': TSIII_from_TOIII_relation,
                  'TOIII_Hagele2006': TOIII_from_TSIII_relation,
                  'TOII_Epm2009': TOII_from_TOIII_relation}

_DEN_FUNC_DICT = {}
