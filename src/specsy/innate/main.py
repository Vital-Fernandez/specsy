import numpy as np
from astropy.io import fits
from pathlib import Path
from warnings import catch_warnings, simplefilter
from arviz import to_netcdf, from_netcdf
from xarray import Dataset
from os import PathLike
from .interpol_pytensor import interpolation_selection
import h5netcdf


def load_inference_data(fname):

    # Load the data
    if isinstance(fname, (str, PathLike, bytes)):
        with catch_warnings():
            simplefilter("ignore", UserWarning)
            inference_data = from_netcdf(fname)
    else:
        inference_data = fname

    return inference_data


def save_inference_data(fname, inference_data, **kwargs):

    if kwargs is not None:

        with catch_warnings():
            simplefilter("ignore", UserWarning)

            for dataset_name, xarrary_dict in kwargs.items():
                inference_data.add_groups({dataset_name: Dataset(xarrary_dict)})

    # Save the adjusted version
    to_netcdf(inference_data, fname)

    return


def load_grids(fname):

    # Container
    grid_dict = None

    # Check there is an input grid file
    if fname is not None:

        # Locate the file
        if fname.is_file():

            grid_dict = {}

            ext = fname.suffix
            print(f'\n- Loading emissivity grid at {fname}')
            if ext == '.fits':
                with fits.open(fname) as hdu_list:
                    for i in range(1, len(hdu_list)):
                        grid_dict[hdu_list[i].name] = hdu_list[i].data
            elif ext == '.nc':
                with h5netcdf.File(fname, 'r') as f:
                    for var_name in f.variables:
                        grid_dict[var_name] = f.variables[var_name][...]

    return grid_dict


def save_grids(fname, grid_dict):

    fname = Path(fname)

    # Case of a fits file:
    ext = fname.suffix
    if fname.suffix == '.fits':

        # Create a primary HDU
        hdu_list = fits.HDUList([fits.PrimaryHDU()])

        # Generate the fits file
        for key, grid in grid_dict.items():
            hdu_i = fits.ImageHDU(grid, name=key)
            hdu_list.append(hdu_i)

        # Write the fits file
        hdu_list.writeto(fname, overwrite=True)

    elif ext == '.nc':

        # Use the first item for the dimensions
        grid_0 = grid_dict[list(grid_dict.keys())[0]]
        m, n = grid_0.shape

        with h5netcdf.File(fname, 'w') as f:

            # Unique dimensions for all the datase
            f.dimensions['m'], f.dimensions['n'] = m, n

            # Iterate over the dictionary and create a variable for each array
            for key, grid in grid_dict.items():
                var = f.create_variable(key, ('m', 'n'), data=grid)
                var.attrs['description'] = f'{key} emissivity'

    else:
        raise KeyError(f'The extension "{ext}" is nto recognized, please use ".nc" or ".fits"')

    print(f'- Emissivity grid saved at: {fname}')

    return


class Innate:

    def __init__(self, grid=None, x_space=None, y_space=None, interpolators='pytensor'):

        # Object attributes
        self.grid = None
        self.interpl = None
        self.lib_interpl = None
        self.x_range, self.y_range = None, None

        # Initiate data and interpolators # TODO Better check
        grid_path = Path(grid)
        self.grid = load_grids(grid_path) if grid_path.is_file() else grid

        if interpolators is not None:

            print(f'\n- Compiling {interpolators} interpolators')

            if interpolators == 'pytensor':
                self.lib_interpl = interpolators
                self.x_range = np.linspace(x_space[0], x_space[1], x_space[2])
                self.y_range = np.linspace(y_space[0], y_space[1], y_space[2])
                self.interpl = interpolation_selection(self.grid, self.x_range, self.y_range, z_range=None,
                                                       interp_type='point')
                print('-- done')

        return
