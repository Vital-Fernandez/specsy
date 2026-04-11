import numpyro
numpyro.set_host_device_count(8)

import jax
# Force JAX to use CPU
jax.config.update("jax_platform_name", "cpu")

import pyneb as pn
import lime

import numpy as np
from examples.previous_version.tools import MultiRegionModel, Region, compile_emis_interp, generate_emissivity_file, emission_flux_model, direct_method_multi_region

# Compile line fluxes
true_params = { "H1":           0,
                "O2":           7.805,
                "O3":           8.055,
                "S2":           5.485,
                "S3":           6.942,
                "Ar3":          5.875,
                "Ar4":          5.215,
                'Ne3':          7.265,
                "N2":           7.655,
                "cHBeta":       0.3,
                "temp_high":    12000.0,
                "temp_low":     10000.0,
                "den_low":      150.0,
                "den_med":      350.0,}

# Declare lines to simulate
input_lines = ['H1_4340A', 'H1_6563A',
               'Ne3_3869A',
               'Ar3_7136A', 'Ar3_7751A',
               'Ar4_4740A',
               'O2_3726A', 'O2_3729A',
               'N2_6548A', 'N2_6583A',
               'O3_4363A', 'O3_4959A', 'O3_5007A',
               'S2_4069A', 'S2_6716A', 'S2_6731A',
               'S3_6312A', 'S3_9068A', 'S3_9530A']
ref_line = 'H1_4861A'

# Input lines frame
lines_df = lime.lines_frame(line_list=input_lines, rejected_lines=[ref_line])
lines_df[['f_lambda', 'line_flux', 'line_err', 'region']] = np.nan, np.nan, np.nan, 'none'
input_lines = lines_df.index.to_numpy()

# Generate the emissivity file
emis_file = '../emissivity_grids_pyneb.nc'
generate_emissivity_file(emis_file, lines_df)

# Extinction
rc = pn.RedCorr()
rc.law = 'CCM89'       # or 'F99', 'S79 H83 CCM89', etc.
rc.R_V = 3.1           # only needed for laws that depend on R_V

X_lambda = rc.X(lines_df.wavelength.to_numpy())
X_Hbeta = rc.X(4861.0)
f_lambda = X_lambda / X_Hbeta - 1.0

# Load the emissivity grids
emis_file = '../emissivity_grids_pyneb.nc'
interp_dict_pc = compile_emis_interp(emis_file)
interp_dict_np = compile_emis_interp(emis_file, array_mode=True)

# Region temperature and density configuration
r0 = Region(name='low', species=['H1', 'S2', 'S3', 'O2', 'N2'],
            temp_mode='free',
            den_mode='free')

r1 = Region(name='high', species=['O3', 'Ar3', 'Ar4', 'Ne3'],
            temp_mode='free',
            den_mode='tied', den_ref='low', den_eq=None)

mult_region = MultiRegionModel([r0, r1])
struc_df = r_str.map_line_structure(input_lines)

# r0 = Region(name='low', species=['H1', 'S2', 'S3'],
#             temp_mode='free',
#             den_mode='free')
#
# r1 = Region(name='high', species=['O2', 'O3'],
#             temp_mode='free',
#             den_mode='tied', den_ref='low', den_eq=None)

# r0 = Region(name='low', species=['H1', 'S2', 'S3'],
#             temp_mode='free',
#             den_mode='free')
#
# r1 = Region(name='med', species=['O2'],
#             temp_mode='tied', temp_ref='high',
#             den_mode='free')
#
# r2 = Region(name='high', species=['O3'],
#             temp_mode='free',
#             den_mode='tied', den_ref='med', den_eq=None)

r_str = MultiRegionModel([r0, r1])
struc_df = r_str.map_line_structure(input_lines)

# Generate the synthetic fluxes
flux_arr = emission_flux_model(input_lines, true_params, struc_df, f_lambda_ref='H1_4861A', emis_interp=interp_dict_np)

struc_df.insert(0, 'line_flux', flux_arr)
struc_df.insert(1, 'line_err', flux_arr * 0.05)
struc_df.insert(2, 'f_lambda', f_lambda)
print(struc_df.to_string())

direct_method_multi_region(struc_df, interp_dict_pc, true_params)


