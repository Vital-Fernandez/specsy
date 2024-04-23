import logging
import numpy as np
import pyneb as pn
from lime import label_decomposition, Line
from inspect import getfullargspec
from scipy.optimize import curve_fit

_logger = logging.getLogger('SpecSy')


def get_pyneb_element(line):

    if line.transition_comp[0] == 'col':
        atom = pn.Atom(line.particle[0].symbol, line.particle[0].ionization)
    elif line.transition_comp[0] == 'rec':
        atom = pn.RecAtom(line.particle[0].symbol, line.particle[0].ionization)
    else:
        _logger.info(f'No emissivity calculation for {line.transition_comp} transition: {line}')
        atom = None

    return atom


def emissivity_grid_calculation(emis_dict, log, temp_grid, den_grid):

    # Container create for observation declaration
    emis_dict = {} if emis_dict is None else emis_dict

    # Grid generation
    temp_range = np.linspace(temp_grid[0], temp_grid[1], temp_grid[2])
    den_range = np.linspace(den_grid[0], den_grid[1], den_grid[2])
    X, Y = np.meshgrid(temp_range, den_range)
    temp_flatten_range, den_flatten_range = X.flatten(), Y.flatten()

    # Loop through the lines and add missing lines
    for line_name in log.index:

        if line_name not in emis_dict:

            line = Line(line_name)
            norm_line = Line(log.loc[line_name, 'norm_line'])

            # Normalization emissivity
            if norm_line not in emis_dict:
                atom_pn = get_pyneb_element(norm_line)
                grid_i = atom_pn.getEmissivity(temp_flatten_range, den_flatten_range, wave=norm_line.wavelength[0], product=False)
                emis_dict[norm_line] = grid_i

            if line.label == 'N2_6584A':
                print('bicho')

            # Calculate and store emissivity ratio
            atom_pn = get_pyneb_element(line)
            if (atom_pn is not None) and (norm_line in emis_dict):
                grid_i = atom_pn.getEmissivity(temp_flatten_range, den_flatten_range, wave=line.wavelength[0], product=False)
                emis_dict[line] = grid_i/emis_dict[norm_line]

    return


class EmissGridGen:

    dict = None

    def __init__(self, log, norm_header='norm_line', temp_range=(9000, 20000, 251), den_range=(1, 600, 101),
                 update_grids=False):

        # Initialize dictionary if first time
        if EmissGridGen.dict is None:
            EmissGridGen.dict = {}

        # Compute the temperature and density range
        self.temp_range = np.linspace(temp_range[0], temp_range[1], temp_range[2])
        self.den_range = np.linspace(den_range[0], den_range[1], den_range[2])

        # Compute 2D arrays
        X, Y = np.meshgrid(temp_range, den_range)
        temp_flatten_range, den_flatten_range = X.flatten(), Y.flatten()

        # Loop through the lines and add missing lines
        for line_name in log.index:

            if (line_name not in EmissGridGen.dict) or update_grids:

                # Target line emissivity
                line = Line(line_name)

                # Calculate and store emissivity ratio
                atom_pn = get_pyneb_element(line)
                grid_i = atom_pn.getEmissivity(temp_flatten_range, den_flatten_range, wave=line.wavelength[0],
                                               product=False)

                # Establish normalization
                norm_line = None
                if norm_header is not None:
                    if norm_header in log.columns:
                        norm_line = Line(log.loc[line_name, norm_header])
                        if (norm_line not in EmissGridGen.dict) or update_grids:
                            atom_pn = get_pyneb_element(norm_line)
                            EmissGridGen.dict[norm_line] = atom_pn.getEmissivity(temp_flatten_range, den_flatten_range,
                                                                                 wave=norm_line.wavelength[0],
                                                                                 product=False)
                    else:
                        raise (f'The input lines log does not have the input normalization header "{norm_header}". '
                               f'The emissivities cannot be not be normalized')

                # Normalize and save
                norm_flux = EmissGridGen.dict[norm_line] if norm_line is not None else 1
                EmissGridGen.dict[line_name] = grid_i/norm_flux

        return


class EmissivitySurfaceFitter:

    def __init__(self):

        self.tempGridFlatten = None
        self.denGridFlatten = None
        self.emisGridDict = None
        self.emisCoeffs = {} # TODO put this one in the superior one

        # Class with the tensor operations of this class
        # TODO we should read this from the xlsx file
        self.ionEmisEq_fit = {'S2_6716A': self.emisEquation_TeDe,
                              'S2_6731A': self.emisEquation_TeDe,
                              'S3_6312A': self.emisEquation_Te,
                              'S3_9069A': self.emisEquation_Te,
                              'S3_9531A': self.emisEquation_Te,
                              'Ar4_4740A': self.emisEquation_Te,
                              'Ar3_7136A': self.emisEquation_Te,
                              'Ar3_7751A': self.emisEquation_Te,
                              'O3_4363A': self.emisEquation_Te,
                              'O3_4959A': self.emisEquation_Te,
                              'O3_5007A': self.emisEquation_Te,
                              'O2_7319A': self.emisEquation_TeDe,
                              'O2_7330A': self.emisEquation_TeDe,
                              'O2_7319A_b': self.emisEquation_TeDe,
                              'N2_6548A': self.emisEquation_Te,
                              'N2_6584A': self.emisEquation_Te,
                              'H1_4102A': self.emisEquation_HI,
                              'H1_4341A': self.emisEquation_HI,
                              'H1_6563A': self.emisEquation_HI,
                              'He1_3889A': self.emisEquation_HeI_fit,
                              'He1_4026A': self.emisEquation_HeI_fit,
                              'He1_4471A': self.emisEquation_HeI_fit,
                              'He1_5876A': self.emisEquation_HeI_fit,
                              'He1_6678A': self.emisEquation_HeI_fit,
                              'He1_7065A': self.emisEquation_HeI_fit,
                              'He1_10830A': self.emisEquation_HeI_fit,
                              'He2_4686A': self.emisEquation_HeII_fit}

        # Initial coeffient values to help with the fitting
        self.epm2017_emisCoeffs = {'He1_3889A': np.array([0.173, 0.00054, 0.904, 1e-5]),
                                   'He1_4026A': np.array([-0.09, 0.0000063, 4.297, 1e-5]),
                                   'He1_4471A': np.array([-0.1463, 0.0005, 2.0301, 1.5e-5]),
                                   'He1_5876A': np.array([-0.226, 0.0011, 0.745, -5.1e-5]),
                                   'He1_6678A': np.array([-0.2355, 0.0016, 2.612, 0.000146]),
                                   'He1_7065A': np.array([0.368, 0.0017, 4.329, 0.0024]),
                                   'He1_10830A': np.array([0.14, 0.00189, 0.337, -0.00027])}

        return

    def emisEquation_Te(self, xy_space, a, b, c):
        temp_range, den_range = xy_space
        return a + b / (temp_range / 10000.0) + c * np.log10(temp_range / 10000)

    def emisEquation_TeDe(self, xy_space, a, b, c, d, e):
        temp_range, den_range = xy_space
        return a + b / (temp_range / 10000.0) + c * np.log10(temp_range / 10000) + np.log10(1 + e * den_range)

    def emisEquation_HI(self, xy_space, a, b, c):
        temp_range, den_range = xy_space
        return a + b * np.log10(temp_range) + c * np.log10(temp_range) * np.log10(temp_range)

    def emisEquation_HeI_fit(self, xy_space, a, b, c, d):
        temp_range, den_range = xy_space
        return (a + b * den_range) * np.log10(temp_range / 10000.0) - np.log10(c + d * den_range)

    def emisEquation_HeII_fit(self, xy_space, a, b):
        temp_range, den_range = xy_space
        return a + b * np.log(temp_range / 10000)

    def fitEmis(self, func_emis, xy_space, line_emis, p0=None):
        p1, p1_cov = curve_fit(func_emis, xy_space, line_emis, p0)
        return p1, p1_cov

    def fitEmissivityPlane(self, linesDF):

        labels_list = linesDF.index.values

        #TODO integrate this function inside the compute emissivity method

        # Dictionary to store the emissivity surface coeffients
        for i in range(labels_list.size):
            lineLabel = labels_list[i]

            # Get equation type to fit the emissivity
            line_func = self.ionEmisEq_fit[lineLabel]
            n_args = len(getfullargspec(
                line_func).args) - 2  # TODO Not working in python 2.7 https://stackoverflow.com/questions/847936/how-can-i-find-the-number-of-arguments-of-a-python-function

            # Compute emissivity functions coefficients
            emis_grid_i = self.emisGridDict[lineLabel]
            p0 = self.epm2017_emisCoeffs[lineLabel] if lineLabel in self.epm2017_emisCoeffs else np.zeros(n_args)
            p1, cov1 = self.fitEmis(line_func, (self.tempGridFlatten, self.denGridFlatten), emis_grid_i, p0=p0)
            self.emisCoeffs[lineLabel] = p1

        return


# TODO undo class and move methods to data reading
class IonEmissivity(EmissivitySurfaceFitter):

    def __init__(self, atomic_references=None, tempGrid=(9000, 20000, 251), denGrid=(1, 600, 101)):

        self.emisGridDict = {}
        self.ftau_coeffs = None
        self.ftau_interp = {}
        self.tempRange = None
        self.denRange = None
        self.tempGridFlatten = None
        self.denGridFlatten = None
        self.normLine = 'H1_4861A'

        EmissivitySurfaceFitter.__init__(self)

        # Defining temperature and density grids
        # FIXME There is a weird error in linspace, check if int(denGrid[0]) can be removed
        if (tempGrid is not None) and (denGrid is not None):
            self.tempRange = np.linspace(int(tempGrid[0]), int(tempGrid[1]), int(tempGrid[2]))
            self.denRange = np.linspace(int(denGrid[0]), int(denGrid[1]), int(denGrid[2]))
            X, Y = np.meshgrid(self.tempRange, self.denRange)
            self.tempGridFlatten, self.denGridFlatten = X.flatten(), Y.flatten()

        # Load user atomic data references
        # TODO upgrade the right method? check for pyneb abundances
        if atomic_references is not None:
            pn.atomicData.defaultDict = atomic_references
            pn.atomicData.resetDataFileDict()

    def get_ions_dict(self, ions_list, recomb_atoms=('H1', 'He1', 'He2'), atomic_references=pn.atomicData.defaultDict):

        # Check if the atomic dataset is the default one
        if atomic_references == pn.atomicData.defaultDict:
            pn.atomicData.resetDataFileDict()
            pn.atomicData.removeFitsPath()
        else:
            pn.atomicData.includeFitsPath()
            pn.atomicData.setDataFileDict(atomic_references)

        # Generate the dictionary with pyneb ions
        ionDict = pn.getAtomDict(ions_list)

        # Remove r extension from recom 'H1', 'He1', 'He2' emissions
        for ion_r in recomb_atoms:
            recomb_key = f'{ion_r}r'
            if recomb_key in ionDict:
                ionDict[ion_r] = ionDict.pop(recomb_key)

        return ionDict

    def load_emis_coeffs(self, line_list, objParams, verbose=True):

        for line in line_list:
            if line in objParams:
                self.emisCoeffs[line] = objParams[line]
            else:
                if verbose:
                    print(f'-- Warning: No emissivity coefficients available for line {line}')
        return

    def computeEmissivityGrids(self, line_labels, ionDict, grids_folder=None, load_grids=False, normLine='H1_4861A', combined_dict={}):

        ion_array, wave_array, latex_array = label_decomposition(line_labels, fit_conf=combined_dict)

        # Generate a grid with the default reference line
        if normLine == 'H1_4861A':
            self.normLine = 'H1_4861A'
            wave_line = 4861.0
            Hbeta_emis_grid = ionDict['H1'].getEmissivity(self.tempRange, self.denRange, wave=wave_line)

        self.emisGridDict = {}
        for i, line_label in enumerate(line_labels):

            # Line emissivity references
            if (grids_folder is not None) and load_grids:
                emis_grid_i = np.load(grids_folder, line_label)

            # Otherwise generate it (and save it)
            else:
                # Single line:
                if ('_m' not in line_label) and ('_b' not in line_label):
                    emis_grid_i = ionDict[ion_array[i]].getEmissivity(self.tempRange, self.denRange, wave=np.round(wave_array[i]))

                # Blended line
                else:
                    emis_grid_i = np.zeros(Hbeta_emis_grid.shape)
                    for component in combined_dict[line_label].split('+'):
                        ion, wave, latex_label = label_decomposition(component, scalar_output=True)
                        emis_grid_i += ionDict[ion_array[i]].getEmissivity(self.tempRange, self.denRange, wave=np.round(wave))

                if grids_folder is not None:
                    np.save(grids_folder, emis_grid_i)

            # Save grid dict
            self.emisGridDict[line_label] = np.log10(emis_grid_i/Hbeta_emis_grid)

        return

    # def compute_ftau_grids(self, ftau_file_path):
    #
    #     """
    #     Correction grids for fluorescence excitation in HeI transitions
    #
    #     :math:`a^2`
    #
    #     :math:`\\alpha > \\beta`
    #
    #
    #     :param ftau_file_path:
    #     :return:
    #     """
    #
    #     # TODO add default path
    #     # Import Optical depth function coefficients
    #     self.ftau_coeffs = import_optical_depth_coeff_table(ftau_file_path)
    #
    #     den_matrix, temp_matrix = np.meshgrid(self.denRange, self.tempRange)
    #
    #     # ftau grid
    #     for lineLabel in self.ftau_coeffs:
    #         if self.ftau_coeffs[lineLabel].sum() != 0.0:
    #
    #             # TODO check equation format for log scale
    #             a, b, c, d = self.ftau_coeffs[lineLabel]
    #             ftau_grid = (a + (b + c * den_matrix + d * np.power(temp_matrix, 2)) * temp_matrix / 10000.0)
    #
    #             ftau_interpolator = RegularGridInterpolator([self.tempRange, self.denRange], ftau_grid[:, :, None], nout=1)
    #
    #             self.ftau_interp[lineLabel] = ftau_interpolator.evaluate
    #
    #     return

    def computeEmissivityEquations(self, linesDF, ionDict, grids_folder=None, load_grids=False, norm_Ion='H1r',
                              norm_pynebCode=4861, linesDb=None):

        labels_list = linesDF.index.values
        ions_list = linesDF.ion.values
        pynebCode_list = linesDF.pynebCode.values
        blended_list = linesDF.blended.values

        # Generate a grid with the default reference line
        Hbeta_emis_grid = ionDict[norm_Ion].getEmissivity(self.tempGridFlatten, self.denGridFlatten, wave=norm_pynebCode,
                                                          product=False)

        for i in range(len(labels_list)):

            # Line emissivity references
            line_label = labels_list[i]

            if (grids_folder is not None) and load_grids:
                emis_grid_i = np.load(grids_folder, line_label)

            # Otherwise generate it (and save it)
            else:

                # Check if it is a blended line:
                if ('_m' not in line_label) and ('_b' not in line_label):
                    # TODO I should change wave by label
                    emis_grid_i = ionDict[ions_list[i]].getEmissivity(self.tempGridFlatten, self.denGridFlatten,
                                                                           wave=float(pynebCode_list[i]), product=False)
                else:
                    emis_grid_i = np.zeros(self.tempGridFlatten.size)
                    for component in blended_list[i].split(','):
                        component_wave = float(linesDb.loc[component].pynebCode)
                        emis_grid_i += ionDict[ions_list[i]].getEmissivity(self.tempGridFlatten, self.denGridFlatten,
                                                                                wave=component_wave, product=False)
                if (grids_folder is not None):
                    np.save(grids_folder, emis_grid_i)

            # Save along the number of points
            self.emisGridDict[line_label] = np.log10(emis_grid_i / Hbeta_emis_grid)

        return
