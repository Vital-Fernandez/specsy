import matplotlib.pyplot as plt
import numpy as np
import specsy as sy
from matplotlib import pyplot

nebCalc = sy.NebularContinua()

# Physical parameters
Te = 10000
Halpha_flux = 1e-14
HeII_HII = 0.01
HeIII_HeII = 0.001

# Nebular continuum for an input wavelength
wave_array = np.linspace(1220, 10000, num=10000)
neb_int = nebCalc.flux_spectrum(wave_array, Te, Halpha_flux, HeII_HII, HeIII_HeII)

fig, ax = plt.subplots()
ax.plot(wave_array, neb_int)
ax.set_yscale('log')
plt.show()
