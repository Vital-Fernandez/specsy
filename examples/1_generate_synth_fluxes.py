import numpy as np
import lime
from innate import load_dataset
from specsy.models.extinction import flambda_calc
from specsy.tools import Region, MultiRegionModel
from specsy.operations.interpolation import compile_bilinear_interp
from specsy.models.literature import _TEM_FUNC_DICT, _DEN_FUNC_DICT


# Synthetic region base parameters
cfg_fname = f'synthetic_region_v0.toml'
lines_fname = f'./synthetic_lines_region_v0.txt'

spec_cfg = lime.load_cfg(cfg_fname)
true_params = spec_cfg['synth_spectrum']['true_params']

# Get the line list
lines_df = lime.lines_frame(line_list=spec_cfg['synth_spectrum']['lines']['line_list'])
input_lines = lines_df.index.to_numpy()

# Normalization with respect to HBeta
ref_line = spec_cfg['synth_spectrum']['lines']['ref_line']
ref_line = lime.Line.from_transition(ref_line)

# Load the emissivity grids
emis_file = load_dataset('./emissivity_grids_pyneb_1.1.30.nc')
interp_dict_np = compile_bilinear_interp(emis_file, array_mode=True)

# Extinction model
f_lambda = flambda_calc(lines_df.wavelength.to_numpy(), norm_wave=ref_line.wavelength,
                        **spec_cfg['synth_spectrum']['extinction'])

# Temperature/density structure
r0 = Region(**spec_cfg['synth_spectrum']['regions']['r0'])
r1 = Region(**spec_cfg['synth_spectrum']['regions']['r1'])
mult_region = MultiRegionModel([r0, r1])
struc_df = mult_region.map_line_structure(input_lines)

# Temperature/density structure arrays
line_arr = struc_df.index.to_numpy()
Tem_label_arr = struc_df.temp.to_numpy()
den_label_arr = struc_df.den.to_numpy()
tem_eq_arr = struc_df.eq_temp.to_numpy()
den_eq_arr = struc_df.eq_den.to_numpy()
particle_arr = struc_df.particle.to_numpy()

# Convenience arrays for the workflow
range_arr = np.arange(len(input_lines))
temp_eq_check = struc_df.eq_temp.isnull().to_numpy()
den_eq_check = struc_df.eq_den.isnull().to_numpy()

# Loop through the transitions and compile the fluxes
flux_arr = np.zeros(len(input_lines))
for i in range_arr:

    # Compute the emissivity
    tem = true_params[Tem_label_arr[i]] if temp_eq_check[i] else _TEM_FUNC_DICT[tem_eq_arr[i]](true_params[Tem_label_arr[i]])
    den = true_params[den_label_arr[i]] if den_eq_check[i] else _DEN_FUNC_DICT[den_eq_arr[i]](true_params[den_label_arr[i]])
    emis = interp_dict_np[line_arr[i]](tem, den)

    # Compute the flux
    if particle_arr[i] == 'H1':
        flux = emis - f_lambda[i] * true_params['cHBeta']
    else:
        flux = true_params[particle_arr[i]] + emis - f_lambda[i] * true_params['cHBeta'] - 12

    flux_arr[i] = flux

# Convert to linear scale
flux_arr = np.power(10, flux_arr)

# Generate the dataframe with the region inputs and structure
struc_df.insert(0, 'line_flux', flux_arr)
struc_df.insert(1, 'line_err', flux_arr * 0.05)
struc_df.insert(2, 'f_lambda', f_lambda)
print(struc_df.to_string())

# Save the synthetic region
lime.save_frame(lines_fname, struc_df)