from logging import getLogger

_logger = getLogger('SpecSy')

from pathlib import Path
from copy import deepcopy
from lime.io import load_cfg
from lime.tools import au
from lime.plotting.format import Themer


# SpecSy figure labels and color formatter
theme_file = Path(__file__).resolve().parent/'specsy_theme.toml'
theme = Themer.from_toml(theme_file)

# Plot formating variables
plot_vars = load_cfg(theme_file, fit_cfg_suffix=None)['display_variables']