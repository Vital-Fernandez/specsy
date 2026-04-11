import numpy as np
import pyneb as pn
from lime import lines_frame
from innate import save_dataset


# Output file name
fname = f'./emissivity_grids_pyneb_{pn.__version__}.nc'

# Dataframe with the target lines
line_selection = ['Ne5_3345A', 'Ne5_3426A', 'H1_3704A', 'H1_3722A', 'O2_3726A',
                  'O2_3729A', 'H1_3734A', 'H1_3750A', 'H1_3771A', 'H1_3798A', 'H1_3835A',
                  'Ne3_3869A', 'He1_3889A', 'H1_3889A', 'Ne3_3968A', 'H1_3970A',
                  'He1_4026A', 'S2_4069A', 'S2_4076A', 'H1_4102A',
                  'H1_4340A', 'Fe2_4358A', 'O3_4363A', 'He1_4388A',
                  'He1_4471A', 'Fe3_4658A', 'He2_4686A', 'Ar4_4711A', 'He1_4713A',
                  'Ar4_4740A', 'H1_4861A', 'He1_4922A', 'Fe3_4925A',
                  'O3_4959A', 'O3_5007A', 'He1_5016A', 'Ar3_5192A',
                  'N1_5198A', 'N1_5200A', 'N2_5755A', 'He1_5876A', 'O1_6300A', 'S3_6312A',
                  'N2_6548A', 'H1_6563A', 'N2_6583A', 'He1_6677A', 'S2_6716A', 'S2_6731A',
                  'He1_7065A', 'Ar3_7136A', 'O2_7319A', 'O2_7320A', 'O2_7330A',
                  'O2_7331A', 'Ar3_7751A', 'H1_8392A', 'H1_8413A', 'H1_8438A', 'H1_8467A',
                  'H1_8502A', 'H1_8545A', 'H1_8598A', 'H1_8665A', 'H1_8750A', 'H1_8863A',
                  'H1_9015A', 'S3_9068A', 'H1_9229A', 'S3_9530A', 'H1_9546A']
lines_df = lines_frame(wave_intvl=(3000, 10000), line_list=line_selection)

# Compute the ranges:
temp_range = np.linspace(6000, 20000, 251)
den_range = np.linspace(1, 600, 101)

# Normalization for the emissivities
H1 = pn.RecAtom('H', 1)
norm_emiss = H1.getEmissivity(temp_range, den_range, wave=4861.0)

# Container output grids
atom_dict, emiss_dict = {}, {}

# Loop through the lines and compute the emissivities
wavelength, particle, t_type = lines_df[['wavelength', 'particle', 'trans']].to_numpy().T
for i, line_name in enumerate(lines_df.index):

    # Get transition atom pyneb object
    print(f'-- {line_name}')
    elem, ionization = particle[i][:-1], particle[i][-1]
    atom = pn.RecAtom(elem, ionization) if t_type[i] == 'rec' else pn.Atom(elem, ionization)

    # Compute and normalize the emissivities
    grid_i = atom.getEmissivity(temp_range, den_range, wave=np.round(wavelength[i]))
    if not isinstance(grid_i, dict):
        emiss_dict[line_name] = grid_i / norm_emiss
    else:
        print(f'Issue computing emissivity for: {line_name} {wavelength[i]} angstroms')

# Data attributes
data_conf = {'parameter': 'emissivity', 'approximation': ('rgi', 'eqn'), 'axes': ('temp', 'den'),
             'temp_range': (9000, 20000, 251), 'den_range': (1, 600, 101)}

# Save the data into a dictionary
save_dataset(fname, emiss_dict, data_conf, custom_cfg=None)


# line_selection = ['Ne5_3345A', 'Ne5_3426A', 'H1_3704A', 'H1_3722A', 'O2_3726A',
#                   'O2_3729A', 'H1_3734A', 'H1_3750A', 'H1_3771A', 'H1_3798A', 'H1_3835A',
#                   'Ne3_3869A', 'He1_3889A', 'H1_3889A', 'Ne3_3968A', 'H1_3970A',
#                   'Fe3_4008A', 'He1_4026A', 'S2_4069A', 'S2_4076A', 'H1_4102A',
#                   'H1_4340A', 'Fe2_4358A', 'O3_4363A', 'He1_4388A',
#                   'He1_4471A', 'Fe3_4658A', 'He2_4686A', 'Ar4_4711A', 'He1_4713A',
#                   'Ar4_4740A', 'H1_4861A', 'Fe3_4881A', 'He1_4922A', 'Fe3_4925A',
#                   'O3_4959A', 'Fe3_4986A', 'O3_5007A', 'He1_5016A', 'Ar3_5192A',
#                   'N1_5198A', 'N1_5200A', 'N2_5755A', 'He1_5876A', 'O1_6300A', 'S3_6312A',
#                   'N2_6548A', 'H1_6563A', 'N2_6583A', 'He1_6677A', 'S2_6716A', 'S2_6731A',
#                   'He1_7065A', 'Ar3_7136A', 'O2_7319A', 'O2_7320A', 'O2_7330A',
#                   'O2_7331A', 'Ar3_7751A', 'H1_8392A', 'H1_8413A', 'H1_8438A', 'H1_8467A',
#                   'H1_8502A', 'H1_8545A', 'H1_8598A', 'H1_8665A', 'H1_8750A', 'H1_8863A',
#                   'H1_9015A', 'S3_9068A', 'H1_9229A', 'S3_9530A', 'H1_9546A']