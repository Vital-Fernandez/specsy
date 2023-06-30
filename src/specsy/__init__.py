"""
Specsy - A python package for the analysis of astronomical spectra
"""

import logging

from .io import label_decomposition, load_log, save_log, load_cfg, save_cfg
from .tools import flux_distribution
from .astro.extinction import cHbeta_from_log
from .astro.chemistry import truncated_SII_density_dist, ratio_S23, sulfur_diaz_2020
from .astro.nebular_continuum import NebularContinua
from .treatement import SpectraSynthesizer

# Creating the lime logger
_logger = logging.getLogger("SpecSy")
_logger.setLevel(logging.INFO)

# Outputting format
consoleHandle = logging.StreamHandler()
consoleHandle.setFormatter(logging.Formatter('%(name)s %(levelname)s: %(message)s'))
_logger.addHandler(consoleHandle)


# Error function SpecSy



