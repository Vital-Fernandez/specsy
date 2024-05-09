import os
import numpy as np
import configparser
from pathlib import Path

import pandas as pd

from lime.transitions import label_decomposition
import lime
from lime.io import load_cfg, save_cfg, save_frame, check_file_dataframe, check_fit_conf
from collections.abc import Sequence
from astropy.io import fits


FITS_INPUTS_EXTENSION = {'lines_list': '20A', 'line_fluxes': 'E', 'line_err': 'E'}

FITS_OUTPUTS_EXTENSION = {'parameter_list': '20A',
                          'mean': 'E',
                          'std': 'E',
                          'median': 'E',
                          'p16th': 'E',
                          'p84th': 'E',
                          'true': 'E'}


class SpecSyError(Exception):
    """SpecSy exception function"""


# Load log
def load_frame(file_address, page: str ='LINELOG', sample_levels: list =['id', 'line'], flux_type=None, lines_list=None,
               norm_line=None):

    # Return
    if isinstance(file_address, pd.DataFrame):
        log = file_address.copy()
    else:
        log = lime.load_frame(file_address, page, sample_levels)

    extra_column = 'line_extract' if norm_line is not None else 'line_flux'

    # Create new column for the lines flux with the requested type (None for user to introduce "line_flux")
    if flux_type is not None:
        lime.tools.extract_fluxes(log, flux_type=flux_type, column_name=extra_column)

    # Check for requested lines and their normalization
    if norm_line is not None:
        lime.tools.normalize_fluxes(log, lines_list, norm_line, flux_column=extra_column, column_name='line_flux')

    return log


def load_HII_CHI_MISTRY_grid(log_scale=False, log_zero_value = -1000):

    # grid_file = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/C17_bb_Teff_30-90_pp.dat'

    # TODO make an option to create the lines and
    grid_file = '/home/vital/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/C17_bb_Teff_30-90_pp.dat'
    lineConversionDict = dict(O2_3726A_m='OII_3727',
                              O3_5007A='OIII_5007',
                              S2_6716A_m='SII_6717,31',
                              S3_9069A='SIII_9069',
                              He1_4471A='HeI_4471',
                              He1_5876A='HeI_5876',
                              He2_4686A='HeII_4686')

    # Load the data and get axes range
    grid_array = np.loadtxt(grid_file)

    grid_axes = dict(OH=np.unique(grid_array[:, 0]),
                     Teff=np.unique(grid_array[:, 1]),
                     logU=np.unique(grid_array[:, 2]))

    # Sort the array according to 'logU', 'Teff', 'OH'
    idcs_sorted_grid = np.lexsort((grid_array[:, 1], grid_array[:, 2], grid_array[:, 0]))
    sorted_grid = grid_array[idcs_sorted_grid]

    # Loop throught the emission line and abundances and restore the grid
    grid_dict = {}
    for i, item in enumerate(lineConversionDict.items()):
        lineLabel, epmLabel = item

        grid_dict[lineLabel] = np.zeros((grid_axes['logU'].size,
                                         grid_axes['Teff'].size,
                                         grid_axes['OH'].size))

        for j, abund in enumerate(grid_axes['OH']):
            idcsSubGrid = sorted_grid[:, 0] == abund
            lineGrid = sorted_grid[idcsSubGrid, i + 3]
            lineMatrix = lineGrid.reshape((grid_axes['logU'].size, grid_axes['Teff'].size))
            grid_dict[lineLabel][:, :, j] = lineMatrix[:, :]

    if log_scale:
        for lineLabel, lineGrid in grid_dict.items():
            grid_logScale = np.log10(lineGrid)

            # Replace -inf entries by -1000
            idcs_0 = grid_logScale == -np.inf
            if np.any(idcs_0):
                grid_logScale[idcs_0] = log_zero_value

            grid_dict[lineLabel] = grid_logScale

    return grid_dict, grid_axes


# Function to save a parameter dictionary into a cfg dictionary
def parseConfDict(output_file, param_dict, section_name, clear_section=False):
    # TODO add logic to erase section previous results
    # TODO add logic to create a new file from dictionary of dictionaries

    # Check if file exists
    if os.path.isfile(output_file):
        output_cfg = configparser.ConfigParser()
        output_cfg.optionxform = str
        output_cfg.read(output_file)
    else:
        # Create new configuration object
        output_cfg = configparser.ConfigParser()
        output_cfg.optionxform = str

    # Clear the section upon request
    if clear_section:
        if output_cfg.has_section(section_name):
            output_cfg.remove_section(section_name)

    # Add new section if it is not there
    if not output_cfg.has_section(section_name):
        output_cfg.add_section(section_name)

    # Map key values to the expected format and store them
    for item in param_dict:
        value_formatted = formatConfEntry(param_dict[item])
        output_cfg.set(section_name, item, value_formatted)

    # Save the text file data
    with open(output_file, 'w') as f:
        output_cfg.write(f)

    return


# Function to map variables to strings
def formatConfEntry(entry_value, float_format=None, nan_format='nan'):
    # TODO this one should be replaced by formatStringEntry
    # Check None entry
    if entry_value is not None:

        # Check string entry
        if isinstance(entry_value, str):
            formatted_value = entry_value

        else:

            # Case of an array
            scalarVariable = True
            if isinstance(entry_value, (Sequence, np.ndarray)):

                # Confirm is not a single value array
                if len(entry_value) == 1:
                    entry_value = entry_value[0]

                # Case of an array
                else:
                    scalarVariable = False
                    formatted_value = ','.join([str(item) for item in entry_value])

            if scalarVariable:

                # Case single float
                print(entry_value)
                if np.isnan(entry_value):
                    formatted_value = nan_format

                else:
                    formatted_value = str(entry_value)

    else:
        formatted_value = 'None'

    return formatted_value


# Function to save the PYMC3 simulation as a fits log
def fits_db(fits_address, model_db, ext_name='', header=None):

    line_labels = model_db['inputs']['lines_list']
    params_traces = model_db['outputs']

    sec_label = 'synthetic_fluxes' if ext_name == '' else f'{ext_name}_synthetic_fluxes'

    # ---------------------------------- Input data

    # Data
    list_columns = []
    for data_label, data_format in FITS_INPUTS_EXTENSION.items():
        data_array = model_db['inputs'][data_label]
        data_col = fits.Column(name=data_label, format=data_format, array=data_array)
        list_columns.append(data_col)

    # Header
    hdr_dict = {}
    for i_line, lineLabel in enumerate(line_labels):
        hdr_dict[f'hierarch {lineLabel}'] = model_db['inputs']['line_fluxes'][i_line]
        hdr_dict[f'hierarch {lineLabel}_err'] = model_db['inputs']['line_err'][i_line]

    # User values:
    for key, value in header.items():
        if key not in ['logP_values', 'r_hat']:
            hdr_dict[f'hierarch {key}'] = value

    # Inputs extension
    cols = fits.ColDefs(list_columns)
    sec_label = 'inputs' if ext_name == '' else f'{ext_name}_inputs'
    hdu_inputs = fits.BinTableHDU.from_columns(cols, name=sec_label, header=fits.Header(hdr_dict))

    # ---------------------------------- Output data
    params_list = model_db['inputs']['parameter_list']
    param_matrix = np.array([params_traces[param] for param in params_list])
    param_col = fits.Column(name='parameters_list', format=FITS_OUTPUTS_EXTENSION['parameter_list'], array=params_list)
    param_val = fits.Column(name='parameters_fit', format='E', array=param_matrix.mean(axis=1))
    param_err = fits.Column(name='parameters_err', format='E', array=param_matrix.std(axis=1))
    list_columns = [param_col, param_val, param_err]

    # Header
    hdr_dict = {}
    for i, param in enumerate(params_list):
        param_trace = params_traces[param]
        hdr_dict[f'hierarch {param}'] = np.mean(param_trace)
        hdr_dict[f'hierarch {param}_err'] = np.std(param_trace)

    for lineLabel in line_labels:
        param_trace = params_traces[lineLabel]
        hdr_dict[f'hierarch {lineLabel}'] = np.mean(param_trace)
        hdr_dict[f'hierarch {lineLabel}_err'] = np.std(param_trace)

    # # Data
    # param_array = np.array(list(params_traces.keys()))
    # paramMatrix = np.array([params_traces[param] for param in param_array])
    #
    # list_columns.append(fits.Column(name='parameter', format='20A', array=param_array))
    # list_columns.append(fits.Column(name='mean', format='E', array=np.mean(paramMatrix, axis=0)))
    # list_columns.append(fits.Column(name='std', format='E', array=np.std(paramMatrix, axis=0)))
    # list_columns.append(fits.Column(name='median', format='E', array=np.median(paramMatrix, axis=0)))
    # list_columns.append(fits.Column(name='p16th', format='E', array=np.percentile(paramMatrix, 16, axis=0)))
    # list_columns.append(fits.Column(name='p84th', format='E', array=np.percentile(paramMatrix, 84, axis=0)))

    cols = fits.ColDefs(list_columns)
    sec_label = 'outputs' if ext_name == '' else f'{ext_name}_outputs'
    hdu_outputs = fits.BinTableHDU.from_columns(cols, name=sec_label, header=fits.Header(hdr_dict))

    # ---------------------------------- traces data
    list_columns = []

    # Data
    for param, trace_array in params_traces.items():
        col_trace = fits.Column(name=param, format='E', array=params_traces[param])
        list_columns.append(col_trace)

    cols = fits.ColDefs(list_columns)

    # Header fitting properties
    hdr_dict = {}
    for stats_dict in ['logP_values', 'r_hat']:
        if stats_dict in header:
            for key, value in header[stats_dict].items():
                hdr_dict[f'hierarch {key}_{stats_dict}'] = value

    sec_label = 'traces' if ext_name == '' else f'{ext_name}_traces'
    hdu_traces = fits.BinTableHDU.from_columns(cols, name=sec_label, header=fits.Header(hdr_dict))

    # ---------------------------------- Save fits files
    hdu_list = [hdu_inputs, hdu_outputs, hdu_traces]

    if fits_address.is_file():
        for hdu in hdu_list:
            try:
                fits.update(fits_address, data=hdu.data, header=hdu.header, extname=hdu.name, verify=True)
            except KeyError:
                fits.append(fits_address, data=hdu.data, header=hdu.header, extname=hdu.name)
    else:
        hdul = fits.HDUList([fits.PrimaryHDU()] + hdu_list)
        hdul.writeto(fits_address, overwrite=True, output_verify='fix')

    return


def save_trace(trace, prior_dict, line_labels, input_fluxes, input_err, inference_model):

    #  ---------------------------- Treat traces and store outputs
    model_params = []
    output_dict = {}
    traces_ref = np.array(trace.varnames)
    for param in traces_ref:

        # Exclude pymc3 variables
        if ('_log__' not in param) and ('interval' not in param):

            trace_array = np.squeeze(trace[param])

            # Restore prior parametrisation
            if param in prior_dict:

                reparam0, reparam1 = prior_dict[param][3], prior_dict[param][4]
                if 'logParams_list' in prior_dict:
                    if param not in prior_dict['logParams_list']:
                        trace_array = trace_array * reparam0 + reparam1
                    else:
                        trace_array = np.power(10, trace_array * reparam0 + reparam1)
                else:
                    trace_array = trace_array * reparam0 + reparam1

                model_params.append(param)
                output_dict[param] = trace_array
                trace.add_values({param: trace_array}, overwrite=True)

            # Line traces
            elif param.endswith('_Op'):

                # Convert to natural scale
                trace_array = np.power(10, trace_array)

                if param == 'calcFluxes_Op':  # Flux matrix case
                    for i in range(trace_array.shape[1]):
                        output_dict[line_labels[i]] = trace_array[:, i]
                    trace.add_values({'calcFluxes_Op': trace_array}, overwrite=True)
                else:  # Individual line
                    ref_line = param[:-3]  # Not saving '_Op' extension
                    output_dict[ref_line] = trace_array
                    trace.add_values({param: trace_array}, overwrite=True)

            # None physical
            else:
                model_params.append(param)
                output_dict[param] = trace_array

    # ---------------------------- Save inputs
    inputs = {'lines_list': line_labels,
              'line_fluxes': input_fluxes,
              'line_err': input_err,
              'parameter_list': model_params}

    # ---------------------------- Store fit
    fit_results = {'models': inference_model, 'trace': trace, 'inputs': inputs, 'outputs': output_dict}

    return fit_results
