from astropy import constants as const
from lime.tools import au
from specsy import _setup_cfg
import re

# Physical constants
L_SUN_CGS = const.L_sun.cgs.value                                       # erg s^-1
HC_CGS = (const.h * const.c).cgs.value                                  # erg cm
HC_EV_ANGSTROM = (const.h * const.c).to(au.eV * au.AA).value            # eV Å (lambda = hc / E)
Q_LIMITS_LAMBDA = {'HI': 911.753, 'HeI': 504.259, 'HeII': 227.838}      # Å

Z_SUN = 0.0134 # Asplund et al. (2009); Starburst99-like grids use 0.020

# Energy (eV) to create the ion of each MIRI line, i.e. the ionization potential of the previous stage
MIRI_ION_LINES = {r'[Ne II] 12.81 $\mu$m': 21.56,
                  r'[S IV] 10.51 $\mu$m': 34.86,
                  r'[Ne III] 15.56 $\mu$m': 40.96,
                  r'[O IV] 25.89 $\mu$m': 54.9,
                  r'[Ar V] 7.90|13.10 $\mu$m': 59.8,
                  r'[Ne V] 14.32|24.32 $\mu$m': 97.1,
                  r'[Ne VI] 7.65 $\mu$m': 126.2}

# ── VMS_BPASS (Martins+25) ───────────────────────────────────────────────────

_VMS_BPASS_DEFAULTS = {'z_keys': {'zsmc':  0.002, 'z0p1':  0.001, 'z0p01': 0.0001,      # pop_vms_ files
                                  'z002':  0.002, 'z001':  0.001, 'z1em4': 0.0001}}     # bpass_sing_ files

VMS_BPASS_CFG = _setup_cfg['stellar']['ssp'].get('VMS_BPASS', _VMS_BPASS_DEFAULTS)

VMS_ZERO_AGE_MYR = 0.1      # The "0myr" models have no logarithm: they are floored at this age

_VMS_PATTERN = re.compile(r'^pop_vms_bpass_(z[\w]+?)_([\dp]+)myr_M(\d+)_(gr|grzsc)(\.neb)?$')

_SING_PATTERN = re.compile(r'^bpass_(sing|bin)_(imf\d+_\d+)\.(z[\w]+)\.([\dp]+)myr(\.neb)?$')

