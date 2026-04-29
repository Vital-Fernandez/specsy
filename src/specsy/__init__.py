import logging
import tomllib

from pathlib import Path

# Creating the lime logger
_logger = logging.getLogger("SpecSy")
_logger.setLevel(logging.INFO)

# Outputting format
consoleHandle = logging.StreamHandler()
consoleHandle.setFormatter(logging.Formatter('%(name)s %(levelname)s: %(message)s'))
_logger.addHandler(consoleHandle)

# Read lime configuration .toml
_inst_dir = Path(__file__).parent
_conf_path = _inst_dir/'specsy.toml'
with open(_conf_path, mode="rb") as fp:
    _setup_cfg = tomllib.load(fp)

__version__ = _setup_cfg['metadata']['version']

from lime import load_cfg as load_cfg, load_frame as load_frame, Line as Line, lines_frame as lines_frame, save_frame as save_frame
from lime import Spectrum as Spectrum, Cube as Cube
from innate import load_dataset as load_dataset
from specsy.observations import Nebula
from specsy.models.extinction import extinction_coeff_calc