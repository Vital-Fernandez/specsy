import matplotlib.pyplot as plt
import numpy as np

import specsy
import specsy as sy
from matplotlib import pyplot, rc_context


specsy.theme.set_style(style='high_res', fig_cfg={'figure.figsize': (10, 10), 'figure.dpi': 400,
                                                  'xtick.color': 'none', 'ytick.color': 'none',  # Hides ticks
                                                    'xtick.labelsize': 0, 'ytick.labelsize': 0,  # Hides tick labels
                                                    'axes.labelsize': 0})
fig_conf = specsy.theme.fig_defaults()

nebCalc = sy.NebularContinua()

# Physical parameters
Te = 10000
Halpha_flux = 1e-14
HeII_HII = 0.01
HeIII_HeII = 0.001

# Nebular continuum for an input wavelength
wave_array = np.linspace(1220, 10000, num=10000)
neb_int = nebCalc.flux_spectrum(wave_array, Te, Halpha_flux, HeII_HII, HeIII_HeII)
gamma_2q, gamma_ff, gamma_fb_HI = nebCalc.components(wave_array, Te, HeII_HII, HeIII_HeII)

# Access the Viridis colormap
cmap = plt.get_cmap('inferno')
mag_color = [0.3, 0.5, 1]

with rc_context(fig_conf):
    fig, ax = plt.subplots()
    ax.plot(wave_array, gamma_2q, color=cmap(mag_color[0]))
    ax.plot(wave_array, gamma_ff, color=cmap(mag_color[1]))
    ax.plot(wave_array, gamma_fb_HI, color=cmap(mag_color[2]))
    ax.set_yscale('log')
    ax.axis('off')
    plt.savefig('/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/nebular_transparent.png', bbox_inches='tight',
                transparent=True)
    plt.savefig('/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/nebular_white.png', bbox_inches='tight',
                transparent=False)