import numpy as np
import lime
from examples.previous_version.tools import compile_emis_interp
from scipy.interpolate import RegularGridInterpolator
import pyneb as pn
import innate

# Declare lines to simulate
input_lines = ['H1_4340A', 'H1_6563A',
               'O2_3726A', 'O2_3729A', 'O2_7319A', 'O2_7330A',
               'O3_4363A', 'O3_4959A', 'O3_5007A',
               'S2_4069A', 'S2_6716A', 'S2_6731A',
               'S3_6312A', 'S3_9068A', 'S3_9530A']
ref_line = 'H1_4861A'

H1 = pn.RecAtom('H', 1)

# Input lines frame
lines_df = lime.lines_frame(line_list=input_lines, rejected_lines=[ref_line])
print(lines_df)

# Compute the ranges:
temp_range = np.linspace(9000, 20000, 251)
den_range = np.linspace(1, 600, 101)

# # Normalization for the emissivities
# H1 = pn.RecAtom('H', 1)
# norm_emiss = H1.getEmissivity(temp_range, den_range, wave=4861.0)
#
# # Container output grids
# atom_dict, emiss_dict = {}, {}
#
# # Loop through the lines and compute the emissivities
# wavelength, particle, t_type = lines_df[['wavelength', 'particle', 'trans']].to_numpy().T
# for i, line_name in enumerate(lines_df.index):
#
#     # Get transition atom pyneb object
#     print(f'-- {line_name}')
#     elem, ionization = particle[i][:-1], particle[i][-1]
#     atom = pn.RecAtom(elem, ionization) if t_type[i] == 'rec' else pn.Atom(elem, ionization)
#
#     # Compute and normalize the emissivities
#     grid_i = atom.getEmissivity(temp_range, den_range, wave=np.round(wavelength[i]))
#     emiss_dict[line_name] = grid_i/norm_emiss
#
# # Data attributes
# data_conf = {'parameter': 'emissivity', 'approximation': ('rgi', 'eqn'), 'axes': ('temp', 'den'),
#              'temp_range': (9000, 20000, 251), 'den_range': (1, 600, 101)}
# trans_conf = {'H1_6563A': {'eqn': 'y~2.3x + 4.3x**2 + 6.3x**3'}}
#
# # Save the data into a dictionary
# output_file = f'./test_emissivity_grids.nc'
# innate.save_dataset(output_file, emiss_dict, data_conf, trans_conf)

output_file = f'../test_emissivity_grids.nc'
interp_dict = compile_emis_interp(output_file, array_mode=True)

emis_set = innate.load_dataset(output_file)

grid_dict = emis_set[0]
grid_params = emis_set[1]

temp_range = emis_set[1]['H1_6563A']['temp_range']
temp_range = np.linspace(*temp_range)

den_range = emis_set[1]['H1_6563A']['den_range']
den_range = np.linspace(*den_range)

emis_matrix = grid_dict['H1_4340A']

sy_interp = RegularGridInterpolator((temp_range, den_range), emis_matrix, method="linear")

# # Converting pymc interpolator to numerical
# pt_interp = make_bilinear_interp(temp_range, den_range, emis_matrix)
# x_sym = tensor.dscalar("x")
# y_sym = tensor.dscalar("y")
# linearTensor_interp = pt_function([x_sym, y_sym], pt_interp(x_sym, y_sym))

tem, den = 12000, 300
print(f'Pyneb result: {H1.getEmissivity(tem, den, wave=4340)/H1.getEmissivity(tem, den, wave=4861):0.4f}')
print(f'Scipy result: {sy_interp([tem, den])[0]:0.4f}')
print(f'Pymc result: {np.power(10, interp_dict['H1_4340A'](tem, den)):0.4f}')
0.4707