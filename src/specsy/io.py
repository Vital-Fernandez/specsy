import logging
import os
import numpy as np
from pandas import DataFrame
import configparser
from pathlib import Path
from collections.abc import Sequence
from astropy.io import fits
from innate import load_dataset
from lime import load_cfg
from arviz import from_netcdf


_logger = logging.getLogger('SpecSy')


FITS_INPUTS_EXTENSION = {'lines_list': '20A', 'line_fluxes': 'E', 'line_err': 'E'}

FITS_OUTPUTS_EXTENSION = {'parameter_list': '20A',
                          'mean': 'E',
                          'std': 'E',
                          'median': 'E',
                          'p16th': 'E',
                          'p84th': 'E',
                          'true': 'E'}


# Load lime configuration
_cfg_file_path = Path(__file__).parent/'specsy.toml'
specsy_cfg = load_cfg(_cfg_file_path)


class SpecSyError(Exception):
    """SpecSy exception function"""


_PYNEB_EMIS_GRID = None  # leading underscore: "module-private"


def load_emis_grid(fname=None):
    global _PYNEB_EMIS_GRID

    if _PYNEB_EMIS_GRID is None:
        if fname is None:
            fname = Path(__file__).parent / f"resources/data/{specsy_cfg['metadata']['emissivity_fname']}"
        _PYNEB_EMIS_GRID = load_dataset(fname)

    return _PYNEB_EMIS_GRID


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


def pack_results(trace, prior_dict, line_labels, input_fluxes, input_err, inference_model):

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


def load_emissivity_interp(fname, array_mode=False):

    emis_set = load_dataset(fname)

    interp_dict = {}
    for trans, emis_matrix in emis_set[0].items():

        temp_range = np.linspace(*emis_set[1][trans]['temp_range'])
        den_range = np.linspace(*emis_set[1][trans]['den_range'])
        log_emis = np.log10(emis_matrix)

        interp = make_bilinear_interp(temp_range, den_range, log_emis)

        if array_mode:
            x_sym = tensor.dscalar("x")
            y_sym = tensor.dscalar("y")
            interp_dict[trans] = pt_function([x_sym, y_sym], interp(x_sym, y_sym))

        else:
            interp_dict[trans] = interp

    return interp_dict


def load_ssp_trace_results(idata, param_names):

    """
    theta, best_model, predictive band (p16, p50, p84), node table and marginal-posterior stats (percentiles,
    and ESS/tau if present) stored in the trace's fit_results group.
    """

    fr = idata.outputs
    theta = [float(v) for v in fr['theta'].values]
    band = tuple(fr[f'band_p{q}'].values for q in (16, 50, 84))
    nodes = DataFrame({'log_age': fr['node_log_age'].values, 'Z': fr['node_Z'].values,
                          'ebv_p50': fr['node_ebv_p50'].values, 'log_A_p50': fr['node_log_A_p50'].values,
                          'prob': fr['node_prob'].values})

    cols = ['p16', 'p50', 'p84'] + (['ess', 'tau'] if 'marginal_ess' in fr else [])
    data = {c: fr[f'marginal_{c}'].values for c in cols}
    stats = DataFrame(data, index=list(param_names))

    return theta, fr['best_model'].values, band, nodes, stats



# def update_fit_results(idata, band_draws=1000, seed=None):
#     """Recomputes the fit_results group (best fit, node table, predictive band, marginal stats) from the posterior."""
#
#     d =  idata.constant_data
#     old = idata.fit_results
#
#     theta_best, best_model, tab = best_fit(idata, d=d)
#     p16, p50, p84 = predictive_band(idata, n_draws=band_draws, seed=seed, d=d)
#     stats = marginal_stats(idata)
#
#     idata['fit_results'] = xr.Dataset(
#         {'theta': ('param', np.asarray(theta_best, dtype=float)),
#          **{key: (old[key].dims, old[key].values) for key in ('p0', 'p0_arr', 'truth') if key in old},
#          'best_model': ('pixel', best_model),
#          'band_p16': ('pixel', p16), 'band_p50': ('pixel', p50), 'band_p84': ('pixel', p84),
#          'node_log_age': ('rank', tab['log_age'].to_numpy()), 'node_Z': ('rank', tab['Z'].to_numpy()),
#          'node_ebv_p50': ('rank', tab['ebv_p50'].to_numpy()), 'node_log_A_p50': ('rank', tab['log_A_p50'].to_numpy()),
#          'node_prob': ('rank', tab['prob'].to_numpy()),
#          'marginal_p16': ('param', stats['p16'].to_numpy()), 'marginal_p50': ('param', stats['p50'].to_numpy()),
#          'marginal_p84': ('param', stats['p84'].to_numpy()),
#          **({'marginal_ess': ('param', stats['ess'].to_numpy()),
#              'marginal_tau': ('param', stats['tau'].to_numpy())} if 'ess' in stats else {})},
#         coords={'param': list(PARAM_NAMES)})
#
#     return idata
#
#
# def chain_chi2(idata, param_names, d=None):
#     """chi2 of every posterior sample (nearest-node model), shape (chain, draw)."""
#
#     d = idata.constant_data if d is None else d
#     theta = np.stack([idata.posterior[name].values for name in param_names], axis=-1)    # (chain, draw, param)
#     flat = theta.reshape(-1, len(param_names))
#     rows = d.rows(flat[:, 0], flat[:, 1])
#
#     good = d.mask & np.isfinite(d.lum) & (d.lum_err > 0)
#     y, err, flux, red = d.lum[good], d.lum_err[good], d.grid_flux[:, good], d.red_corr[good]
#
#     chi2 = np.full(len(flat), np.inf)
#     for i, (row, (_, _, ebv, log_amp)) in enumerate(zip(rows, flat)):
#         if row >= 0:
#             chi2[i] = np.sum(((y - 10.0 ** log_amp * flux[row] * red ** ebv) / err) ** 2)
#
#     return chi2.reshape(theta.shape[:2])
#
#
# def stuck_chains(idata, k=1.5, dchi2=20.0):
#     """
#     Boolean mask of the chains to keep. A chain is dropped only if its median chi2 is both more than `dchi2` above
#     the best chain and above the outlier fence Q3 + k * IQR of the chain medians.
#     """
#
#     med = np.median(chain_chi2(idata), axis=1)
#     q1, q3 = np.percentile(med, (25, 75))
#     keep = (med - med.min() <= dchi2) | (med <= q3 + k * (q3 - q1))
#
#     if not keep.all():
#         _logger.warning(f'{int((~keep).sum())} of the {keep.size} chains dropped as stuck: median chi2 above the '
#                         f'best chain by {np.round(np.sort(med[~keep] - med.min()), 1).tolist()}.')
#
#     return keep

def load_trace(trace_pname, drop_stuck=False, k=1.5, dchi2=20.0, band_draws=1000, seed=None):
    """
    Loads a SSP_sampler trace and applies the post-sampling corrections.

    drop_stuck: remove the chains stuck far from the others (see stuck_chains for k and dchi2). The fit results
    are recomputed from the remaining chains, using `band_draws` and `seed` for the predictive band.
    """

    idata = from_netcdf(trace_pname)

    # if drop_stuck:
    #     keep = stuck_chains(idata, k, dchi2)
    #     if not keep.all():
    #         if keep.sum() < 2:
    #             raise SpecSyError(f'Only {int(keep.sum())} chain(s) left after removing the stuck chains from '
    #                               f'{trace_pname}. Try a larger dchi2 or k, or drop_stuck=False.')
    #         idata = idata.isel(chain=np.flatnonzero(keep))
    #         update_fit_results(idata, band_draws=band_draws, seed=seed)

    return idata
