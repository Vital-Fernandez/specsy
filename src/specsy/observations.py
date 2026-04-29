import logging
from dataclasses import dataclass, field
from typing import Callable, Optional, Union
from pandas import DataFrame
from lime import Line
from lime.io import check_fit_conf, check_file_dataframe
from specsy.io import SpecSyError
from specsy.inference import ModelManager

@dataclass
class RegionParam:
    """
    Defines how a single variable (temperature or density) is specified.

    - "free":  provide distribution and kwargs. The variable is sampled.
    - "tied":  provide ref (name or index) to link to another region.
               Optionally provide an equation key (string) to transform the
               parent value. If no equation is provided, the value is copied as-is.
    """

    label: str

    mode: str  # "free" | "tied"

    # Distribution for free variable and its parameters
    distr: Optional[Callable] = None
    kwargs: dict = field(default_factory=dict)

    # Parent region for the parameter and key for the empirical/theoretical relation if necessary
    ref: Optional[Union[str, int]] = None
    eq: Optional[str] = None  # key into equations dict


class Region:

    def __init__(self, name, species, temp_mode, temp_dist=None, temp_kwargs=None, temp_ref=None, temp_eq=None,
                 den_mode=None, den_dist=None, den_kwargs=None, den_ref=None, den_eq=None):

        self.name = name
        self.species = species
        self.temp = RegionParam('temp', temp_mode, distr=temp_dist, kwargs=temp_kwargs, ref=temp_ref, eq=temp_eq)
        self.den = RegionParam('den', den_mode, distr=den_dist, kwargs=den_kwargs, ref=den_ref, eq=den_eq)

        return


class IonizationStructure:

    def __init__(self, cfg):

        # Attributes
        self.n_regions = None
        self.region_list = None
        self.region_map = None

        # Loop through the input file regions and list them
        self.region_list = []
        for i, (region_label, region_params) in enumerate(cfg.items()):
            region = Region(**region_params)
            self.region_list.append(region)

        self.n_regions = len(self.region_list)
        self.region_map: dict = {r.name: r for r in self.region_list}

        return

    def __len__(self):
        return self.n_regions

    def __iter__(self):
        return iter(self.region_list)

    def __getitem__(self, key):
        return self.region_map[key]

    def map_line_structure(self, lines_frame, temp_label='temp', den_label='den', norm_line=None):

        # Add extra columns
        lines_frame['region'] = '-'
        lines_frame[temp_label] = '-'
        lines_frame[den_label] = '-'
        lines_frame['eq_temp'] = '-'
        lines_frame['eq_den'] = '-'

        # Map the regions
        for region in self.region_list:
            idcs = lines_frame.particle.isin(region.species)
            lines_frame.loc[idcs, 'region'] = region.name

        # Recover the temperature and density structure
        for region in self.region_list:
            idcs = lines_frame.region == region.name
            for param_label in [temp_label, den_label]:
                param = getattr(region, param_label)
                if param.mode == "free":
                    lines_frame.loc[idcs, param_label] = f'{param.label}_{region.name}'
                else:
                    if param.ref is None:
                        raise SpecSyError(f'In region "{region.name}" the parameter "{param_label}" is defined as "{param.mode}" to '
                                          f'another region but this region is not defined using "{param_label}_ref" on '
                                          f'the "{region.name}" configuration entry.')
                    lines_frame.loc[idcs, param_label] = f'{param.label}_{param.ref}'
                    if param.eq is not None:
                        lines_frame.loc[idcs, f'eq_temp'] = param.eq

        # Remove the normalization line
        if (norm_line is not None) and (norm_line in lines_frame.index):
            lines_frame = lines_frame.drop(norm_line)

        return lines_frame

class Nebula:

    def __init__(self, lines_frame: DataFrame, ionization_structure: IonizationStructure):

        # Attributes
        self.lines_frame = lines_frame
        self.ion_struct = ionization_structure
        self.sampled: dict = {}

        # Limits for the number of ionization regions
        assert 1 <= len(self.ion_struct) <= 4, "Between 1 and 4 regions allowed"

        # Functionality attributes
        self.infer = ModelManager(self.lines_frame, self.ion_struct)


        return

    @classmethod
    def from_ionization_structure(cls, cfg_params, nebula_name, region_name='regions'):

        # Open the configuration
        cfg_params = check_fit_conf(cfg_params, default_key=None, obj_key=None, fit_cfg_suffix=None)

        # Loop through the input file regions and list them
        region_list = []
        for i, (region_label, region_params) in enumerate(cfg_params[nebula_name][region_name].items()):
            region = Region(**region_params)
            region_list.append(region)

        return cls(region_list)

    @classmethod
    def from_lines_frame(cls, lines_frame, fit_cfg=None, flux_entry='profile', default_cfg_prefix='default',
                         obj_cfg_prefix=None, update_default=False):

        # Check if the lines log is a dataframe or a file address
        bands = check_file_dataframe(lines_frame, copy_input=True)
        if 'line_flux' not in bands.columns:
            flux_col = f'{flux_entry}_flux'
            err_col = f'{flux_entry}_flux_err'
            if flux_col in bands.columns and err_col in bands.columns:
                bands.insert(0, 'line_flux', lines_frame[flux_col])
                bands.insert(1, 'line_flux_err', lines_frame[err_col])
            else:
                raise SpecSyError(f'The input lines frame does not contain "line_flux"/"line_flux_err" or "{flux_col}"/"{err_col}" columns')

        # Import the ionization structure configuration keep track of the global / local configuration
        input_conf = check_fit_conf(fit_cfg, default_cfg_prefix, obj_cfg_prefix, update_default, fit_cfg_suffix='_ionization_structure')

        if 'region' not in input_conf['default_ionization_structure']:
            msg = (f'The input "fit_cfg" argument does not contain a "region" entry for the ionization entry:'
                   f'\n - The queried sections are: f"{default_cfg_prefix}_ionization_structure" '
                   f'{"" if obj_cfg_prefix is None else f"{default_cfg_prefix}_ionization_structure"}')
            SpecSyError(msg)

        # Recover the ionization structure
        ionization_structure = IonizationStructure(input_conf['default_ionization_structure']['region'])

        return cls(bands, ionization_structure)

    def build(self):

        # Sample the distribution value if free
        for region in self.ion_struct:
            self.sampled[region.name] = {}

            if region.temp.mode == "free":
                self.sampled[region.name]["temp"] = region.temp.distr(**region.temp.kwargs)

            if region.den.mode == "free":
                self.sampled[region.name]["den"] = region.den.distr(**region.den.kwargs)

        return self

    def get(self, region_name: str, var_type: str):

        # Recover the requested parameter for the physical region
        region = self.region_map[region_name]
        var_spec = getattr(region, var_type)
        equations = getattr(self, f"{var_type}_eq_dict")

        match var_spec.mode:

            # Free variable
            case "free":
                return self.sampled[region_name][var_type]

            # Tied to another example
            case 'tied':
                ref_name = self.ion_struct[var_spec.ref].name if isinstance(var_spec.ref, int) else var_spec.ref
                parent_value = self.get(ref_name, var_type)  # recursive call

                if var_spec.eq is not None:
                    return equations[var_spec.eq](parent_value)
                else:
                    return parent_value

            case _:
                raise ValueError(f"Unknown mode '{var_spec.mode}'. Use 'free' or 'tied'.")


