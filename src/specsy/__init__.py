"""
Specsy - A python package for the analysis of astronomical spectra
"""

import logging
from pathlib import Path


# Creating the lime logger
_logger = logging.getLogger("SpecSy")
_logger.setLevel(logging.INFO)

# Outputting format
consoleHandle = logging.StreamHandler()
consoleHandle.setFormatter(logging.Formatter('%(name)s %(levelname)s: %(message)s'))
_logger.addHandler(consoleHandle)

from .tools import flux_distribution
from .models.extinction import cHbeta_from_log
from .models.chemistry import truncated_SII_density_dist, ratio_S23, sulfur_diaz_2020


