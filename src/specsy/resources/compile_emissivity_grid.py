import numpy as np
import pyneb as pn
from pathlib import Path
from lime import lines_frame
from innate import save_dataset, load_dataset
from specsy import cfg, __version__
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import RegularGridInterpolator

# Settings for the grid
specsy_version =__version__
emis_cfg = cfg['emissivity_grid']

# Get the lines
line_selection = emis_cfg['line_selection']
lines_df = lines_frame(line_list=line_selection)

# Compute the temperature and density ranges
temp_min, temp_max, temp_steps = emis_cfg['temp_min'], emis_cfg['temp_max'], emis_cfg['temp_steps']
den_min, den_max, den_steps = emis_cfg['den_min'], emis_cfg['den_max'], emis_cfg['den_steps']

temp_range = np.linspace(temp_min, temp_max, temp_steps)
den_range = np.linspace(den_min, den_max, den_steps)

# Normalization wavelength
norm_wave = emis_cfg['norm_wave']
H1 = pn.RecAtom('H', 1)
norm_emiss = H1.getEmissivity(temp_range, den_range, wave=norm_wave)

# Folder for the individual emissivity grids
grids_folder = Path(__file__).parent/'emissivity_grids'
grids_folder.mkdir(exist_ok=True)

# # Loop through the lines and compute the emissivities
# atom_dict, emiss_dict = {}, {}
# wavelength, particle, t_type = lines_df[['wavelength', 'particle', 'trans']].to_numpy().T
# for i, line_name in enumerate(lines_df.index):
#
#     print(f'-- {line_name}')
#     grid_path = grids_folder/f'{line_name}.npy'
#
#     # Load the grid if it already exists
#     if grid_path.is_file():
#         emiss_dict[line_name] = np.load(grid_path)
#         continue
#
#     # Get transition atom pyneb object
#     elem, ionization = particle[i][:-1], particle[i][-1]
#     atom = pn.RecAtom(elem, ionization) if t_type[i] == 'rec' else pn.Atom(elem, ionization)
#
#     # Compute and normalize the emissivities
#     grid_i = atom.getEmissivity(temp_range, den_range, wave=np.round(wavelength[i]))
#     if not isinstance(grid_i, dict):
#         emiss_dict[line_name] = grid_i / norm_emiss
#         np.save(grid_path, emiss_dict[line_name])
#     else:
#         print(f'Issue computing emissivity for: {line_name} {wavelength[i]} angstroms')
#
# # Data attributes
# data_conf = {'parameter': 'emissivity', 'approximation': ('rgi', 'eqn'), 'axes': ('temp', 'den'),
#              'temp_range': (temp_min, temp_max, temp_steps), 'den_range': (den_min, den_max, den_steps)}
#
# # Output file
# fname = Path(__file__).parent/f'data/emissivity_grids_pyneb_{pn.__version__}_v{specsy_version}.nc'
# save_dataset(fname, emiss_dict, data_conf, custom_cfg=None)


# Load the saved dataset and validate against fresh PyNeb computations
fname = Path(__file__).parent/f'data/emissivity_grids_pyneb_{pn.__version__}_v{specsy_version}.nc'
emiss_dict, data_conf = load_dataset(fname)[:2]

n_points = 500
rng = np.random.default_rng(42)
temp_rand = rng.uniform(temp_min, temp_max, n_points)
den_rand = rng.uniform(den_min, den_max, n_points)
points_rand = np.column_stack((temp_rand, den_rand))

# Direct PyNeb normalization emissivity at the random points (paired, not product grid)
norm_emiss_rand = H1.getEmissivity(temp_rand, den_rand, wave=norm_wave, product=False)

wavelength, particle, t_type = lines_df[['wavelength', 'particle', 'trans']].to_numpy().T

for i, line_name in enumerate(lines_df.index):

    if line_name not in emiss_dict:
        continue

    print(f'-- Validation plot: {line_name}')

    # Interpolated values from the SAVED grid
    grid_saved = np.asarray(emiss_dict[line_name])
    rgi = RegularGridInterpolator((temp_range, den_range), grid_saved)
    interp_vals = rgi(points_rand)

    # Direct PyNeb values at the same points (the reference)
    elem, ionization = particle[i][:-1], particle[i][-1]
    atom = pn.RecAtom(elem, ionization) if t_type[i] == 'rec' else pn.Atom(elem, ionization)
    direct_vals = atom.getEmissivity(temp_rand, den_rand, wave=np.round(wavelength[i]),
                                     product=False) / norm_emiss_rand

    # Relative difference in percent
    rel_diff_pct = 100 * np.abs(interp_vals - direct_vals) / np.abs(direct_vals)

    # Scatter plot with fixed 1% - 100% log color scale
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(temp_rand, den_rand, c=rel_diff_pct, cmap='viridis',
                    norm=LogNorm(vmin=1, vmax=100), s=25)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('Relative difference (%)')
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'Density (cm$^{-3}$)')
    ax.set_title(f'{line_name} (max diff: {rel_diff_pct.max():.3f}%)')

    fig.savefig(grids_folder/f'{line_name}_validation.png', dpi=150, bbox_inches='tight')
    plt.close(fig)