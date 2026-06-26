import numpy as np
from pathlib import Path
from specsy.models.stellar import load_pySB99_config, pystarburst99_file_manager
import time

# ── User settings ──────────────────────────────────────────────────────────────
pySB99_source_path = Path('/home/vital/Astrodata/pystarburst99')  # adjust as needed
out_folder = Path('/home/vital/Astrodata/pystarburst99/compiled_spectra/tutorial_model')

M_total         = 1.0e6                         # Normalization
IMF_exponents   = [1.3, 2.3]
IMF_mass_limits = (0.1, 0.5, 120.)
run_speed_mode  = 'DEFAULT'
SED_library     = 'FW'                         # WMBasic low-res
spectra_library = 'WM'
rot             = False

# Target metallicities and log-ages
# Z_list     = ['IZw18', 'SMC', 'LMC', 'MW', 'MWC']   # Z = 0.0004, 0.002, 0.006, 0.014, 0.020
Z_list       = ['SMC']   # Z = 0.0004, 0.002, 0.006, 0.014, 0.020
log_ages_out = [6.0, 6.602, 6.845, 7.]   # log(yr)
times_out_yr = 10**np.array(log_ages_out)             # linear years

# All plot flags off — we only want the spectra arrays
plot_ion_flux = plot_wind = plot_uv_slope = plot_ew = False
plot_colours  = plot_SED_with_time = plot_hires_spectra = plot_SN_rate = False
save_output   = True

# Loop through the metallicities
for i, metallicity in enumerate(Z_list):
    metal_cfg = load_pySB99_config(metallicity, pySB99_source_path, rot=False)
    wave, fluxes = pystarburst99_file_manager(cfg=metal_cfg, M_total=M_total, IMF_exponents=IMF_exponents,
                                              IMF_mass_limits=IMF_mass_limits, run_speed_mode=run_speed_mode,
                                              times_out_yr=times_out_yr)
    np.savetxt(out_folder / f'pySB99_specsy_{metallicity}_v0.txt', np.array(fluxes).T)

    print(fluxes)

