import numpy as np
import specsy as sy

from pathlib import Path
from specsy.models.chemistry import TOIII_from_TSIII_relation
from lime.transitions import log_from_line_list

# Load fitting configuration file
fit_cfg_file = Path(f'data/fit_confg.cfg')
fit_cfg = sy.load_cfg(fit_cfg_file)

# Set the paramter values
fit_cfg['true_values'] = {'n_e': 150.0,
                          'T_low': 12500.0,
                          'T_high': TOIII_from_TSIII_relation(12500.0),
                          'tau': 0.60 + 0.15,
                          'cHbeta': 0.08 + 0.02,
                          'H1': 0.0,
                          'He1': np.log10(0.070 + 0.005),
                          'He2': np.log10(0.00088 + 0.0002),
                          'O2': 7.805 + 0.15,
                          'O3': 8.055 + 0.15,
                          'N2': 5.845 + 0.15,
                          'S2': 5.485 + 0.15,
                          'S3': 6.94 + 0.15,
                          'Ne3': 7.065 + 0.15,
                          'Fe3': 5.055 + 0.15,
                          'Ar3': 5.725 + 0.15,
                          'Ar4': 5.065 + 0.15}

# Declare lines to simulate
input_lines = ['H1_4340A', 'H1_4861A', 'H1_6563A',
               'He1_4026A', 'He1_4471A', 'He1_5876A', 'He1_6678A', 'He1_7065A', 'He2_4685A',
               'O2_3726A_m', 'O2_7319A_m', 'O3_4363A', 'O3_4959A', 'O3_5007A',
               'Ne3_3869A',
               'N2_6548A', 'N2_6583A',
               'S2_6716A', 'S2_6731A', 'S3_6312A', 'S3_9068A', 'S3_9530A',
               'Ar3_7136A', 'Ar4_4740A',
               'Fe3_4658A']

merged_lines = {'O2_3726A_m': 'O2_3726A-O2_3729A',
                'O2_7319A_m': 'O2_7319A-O2_7330A'}

# Generate table with the lines data
line_log = log_from_line_list(input_lines, comps_dict=merged_lines)
