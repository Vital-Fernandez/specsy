import itertools
import numpy as np
import pytensor.tensor as tt
from ..io import load_HII_CHI_MISTRY_grid, label_decomposition
# from ..models.emissivity import IonEmissivity
#
#
# def emissivity_grid_calc(lines_array, comp_dict, temp_grid_points=(9000, 20000, 251), den_grid_points=(1, 600, 101)):
#
#     print(f'- Computing emissivity grids for {len(lines_array)} lines\n')
#
#     # Compute the atomic data grids
#     objIons = IonEmissivity(tempGrid=temp_grid_points, denGrid=den_grid_points)
#
#     ion_array, wave_array, latex_array = label_decomposition(lines_array, fit_conf=comp_dict)
#
#     # Define the dictionary with the pyneb ion objects
#     ionDict = objIons.get_ions_dict(ion_array)
#
#     # Compute the emissivity surfaces for the observed emission lines
#     objIons.computeEmissivityGrids(lines_array, ionDict, combined_dict=comp_dict)
#
#     # Compile exoplanet interpolator functions so they can be used wit numpy
#     emisGridInterpFun = gridInterpolatorFunction(objIons.emisGridDict, objIons.tempRange, objIons.denRange)
#
#     print(f'-- completed')
#
#     return emisGridInterpFun


def gridInterpolatorFunction(interpolatorDict, x_range, y_range, z_range=None, interp_type='point'):

    emisInterpGrid = {}

    if interp_type == 'point':
        for line, emisGrid_i in interpolatorDict.items():
            emisInterp_i = RegularGridInterpolator([x_range, y_range], emisGrid_i[:, :, None], nout=1)
            emisInterpGrid[line] = emisInterp_i.evaluate

    elif interp_type == 'axis':
        for line, emisGrid_i in interpolatorDict.items():
            emisGrid_i_reshape = emisGrid_i.reshape((x_range.size, y_range.size, -1))
            emisInterp_i = RegularGridInterpolator([x_range, y_range], emisGrid_i_reshape)
            emisInterpGrid[line] = emisInterp_i.evaluate

    elif interp_type == 'cube':
        for line, grid_ndarray in interpolatorDict.items():
            xo_interp = RegularGridInterpolator([x_range, y_range, z_range], grid_ndarray)
            emisInterpGrid[line] = xo_interp.evaluate

    return emisInterpGrid


def as_tensor_variable(x, dtype="float64", **kwargs):
    t = tt.as_tensor_variable(x, **kwargs)
    if dtype is None:
        return t
    return t.astype(dtype)


def regular_grid_interp(points, values, coords, *, fill_value=None):
    """Perform a linear interpolation in N-dimensions w a regular grid

    The data must be defined on a filled regular grid, but the spacing may be
    uneven in any of the dimensions.

    This implementation is based on the implementation in the
    ``scipy.interpolate.RegularGridInterpolator`` class which, in turn, is
    based on the implementation from Johannes Buchner's ``regulargrid``
    package https://github.com/JohannesBuchner/regulargrid.


    Args:
        points: A list of vectors with shapes ``(m1,), ... (mn,)``. These
            define the grid points in each dimension.
        values: A tensor defining the values at each point in the grid
            defined by ``points``. This must have the shape
            ``(m1, ... mn, ..., nout)``.
        coords: A matrix defining the coordinates where the interpolation
            should be evaluated. This must have the shape ``(ntest, ndim)``.
    """
    points = [as_tensor_variable(p) for p in points]
    ndim = len(points)
    values = as_tensor_variable(values)
    coords = as_tensor_variable(coords)

    # Find where the points should be inserted
    indices = []
    norm_distances = []
    out_of_bounds = tt.zeros(coords.shape[:-1], dtype=bool)
    for n, grid in enumerate(points):
        x = coords[..., n]
        i = tt.extra_ops.searchsorted(grid, x) - 1
        out_of_bounds |= (i < 0) | (i >= grid.shape[0] - 1)
        i = tt.clip(i, 0, grid.shape[0] - 2)
        indices.append(i)
        norm_distances.append((x - grid[i]) / (grid[i + 1] - grid[i]))

    result = tt.zeros(tuple(coords.shape[:-1]) + tuple(values.shape[ndim:]))
    for edge_indices in itertools.product(*((i, i + 1) for i in indices)):
        weight = tt.ones(coords.shape[:-1])
        for ei, i, yi in zip(edge_indices, indices, norm_distances):
            weight *= tt.where(tt.eq(ei, i), 1 - yi, yi)
        result += values[edge_indices] * weight

    if fill_value is not None:
        result = tt.switch(out_of_bounds, fill_value, result)

    return result


class GridWrapper:

    def __init__(self, grid_address=None):

        self.grid_LineLabels = None
        self.grid_emissionFluxes = None
        self.grid_emissionFluxErrs = None

        self.gridInterp = None
        self.idx_analysis_lines = None
        self.grid_array = None

        return

    def ndarray_from_DF(self, grid_DF, axes_columns=None, data_columns='all', sort_axes=True, dict_output=True,
                        empty_value=np.nan):

        if sort_axes:
            assert set(axes_columns).issubset(set(grid_DF.columns.values)), f'- Error: Mesh grid does not include all' \
                                                                            f' input columns {axes_columns}'
            grid_DF.sort_values(axes_columns, inplace=True)

        # Compute axes coordinates for reshaping
        axes_cords = {}
        reshape_array = np.zeros(len(axes_columns)).astype(int)
        for i, ax_name in enumerate(axes_columns):
            axes_cords[ax_name] = np.unique(grid_DF[ax_name].values)
            print(ax_name, axes_cords[ax_name], f'length ({len(axes_cords[ax_name])})')
            reshape_array[i] = axes_cords[ax_name].size

        # Declare grid data columns
        if data_columns == 'all':
            data_columns = grid_DF.columns[~grid_DF.columns.isin(axes_columns)].values
        axes_cords['data'] = data_columns

        # mesh_dict
        if dict_output:
            output_container = {}
            for i, dataColumn in enumerate(data_columns):
                data_array_flatten = grid_DF[dataColumn].values
                output_container[dataColumn] = data_array_flatten.reshape(reshape_array.astype(int))

        # mesh_array
        else:
            output_container = np.full(np.hstack((reshape_array, len(data_columns))), np.nan)
            for i, dataColumn in enumerate(data_columns):
                data_array_flatten = grid_DF[dataColumn].values
                output_container[..., i] = data_array_flatten.reshape(reshape_array.astype(int))

        return output_container, axes_cords

    def generate_xo_interpolators(self, grid_dict, axes_list, axes_coords, interp_type='point', empty_value=np.nan):

        # Establish interpolation axes: (x_range, y_range, z_range,...)
        ax_range_container = [None] * len(axes_list)
        for i, ax in enumerate(axes_list):
            ax_range_container[i] = axes_coords[ax]

        if interp_type == 'point':

            output_container = {}

            for grid_key, grid_ndarray in grid_dict.items():
                xo_interp = RegularGridInterpolator(ax_range_container, grid_ndarray)
                output_container[grid_key] = xo_interp.evaluate

            return output_container

        if interp_type == 'axis':

            # Generate empty grid from first data element
            grid_shape = list(grid_dict[axes_coords['data'][0]].shape) + [len(axes_coords['data'])]
            data_grid = np.full(grid_shape, empty_value)
            for i, label_dataGrid in enumerate(axes_coords['data']):
                data_grid[..., i] = grid_dict[label_dataGrid]

            # Add additional dimension with -1 for interpolation along axis
            reShapeDataGrid_shape = [len(item) for item in ax_range_container] + [-1]
            xo_interp = RegularGridInterpolator(ax_range_container, data_grid.reshape(reShapeDataGrid_shape))

            return xo_interp.evaluate


    def HII_Teff_models(self, obsLines, obsFluxes, obsErr):

        gridLineDict, gridAxDict = load_HII_CHI_MISTRY_grid(log_scale=True)
        self.gridInterp = gridInterpolatorFunction(gridLineDict,
                                                   gridAxDict['logU'],
                                                   gridAxDict['Teff'],
                                                   gridAxDict['OH'],
                                                   interp_type='cube')

        # Add merged lines
        if ('S2_6716A' in obsLines) and ('S2_6731A' in obsLines) and ('S2_6716A_m' not in obsLines):

            # Rename the grid label to match observable
            self.gridInterp['S2_6716A'] = self.gridInterp.pop('S2_6716A_m')

            lines_Grid = np.array(list(self.gridInterp.keys()))
            self.idx_analysis_lines = np.in1d(obsLines, lines_Grid)

            # Use different set of fluxes for direct method and grids
            self.grid_LineLabels = obsLines.copy()
            self.grid_emissionFluxes = obsFluxes.copy()
            self.grid_emissionFluxErrs = obsErr.copy()

            # Compute the merged line
            i_S2_6716A, i_S2_6731A = obsLines == 'S2_6716A', obsLines == 'S2_6731A'
            S2_6716A_m_flux = obsFluxes[i_S2_6716A][0] + obsFluxes[i_S2_6731A][0]
            S2_6716A_m_err = np.sqrt(obsErr[i_S2_6716A][0]**2 + obsErr[i_S2_6731A][0]**2)

            # Replace conflicting flux
            self.grid_emissionFluxes[i_S2_6716A] = S2_6716A_m_flux
            self.grid_emissionFluxErrs[i_S2_6716A] = S2_6716A_m_err

        else:
            lines_Grid = np.array(list(gridLineDict.keys()))
            self.idx_analysis_lines = np.in1d(obsLines, lines_Grid)
            self.grid_LineLabels = obsLines.copy()
            self.grid_emissionFluxes = obsFluxes.copy()
            self.grid_emissionFluxErrs = obsErr.copy()

        return


class RegularGridInterpolator:

    """Linear interpolation on a regular grid in arbitrary dimensions

    The data must be defined on a filled regular grid, but the spacing may be
    uneven in any of the dimensions.

    Args:
        points: A list of vectors with shapes ``(m1,), ... (mn,)``. These
            define the grid points in each dimension.
        values: A tensor defining the values at each point in the grid
            defined by ``points``. This must have the shape
            ``(m1, ... mn, ..., nout)``.
    """

    def __init__(self, points, values, fill_value=None, **kwargs):
        self.ndim = len(points)
        self.points = points
        self.values = values
        self.fill_value = fill_value

    def evaluate(self, t):
        """Interpolate the data

        Args:
            t: A matrix defining the coordinates where the interpolation
                should be evaluated. This must have the shape
                ``(ntest, ndim)``.
        """
        return regular_grid_interp(self.points, self.values, t, fill_value=self.fill_value)