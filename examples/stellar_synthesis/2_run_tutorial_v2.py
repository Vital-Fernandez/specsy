import sesamme.models as models
import sesamme.mcmc as stats
import sesamme.vis as vis
import emcee
import numpy as np
import corner
from astropy.table import Table
from astropy import units as u
from sesamme.models import load_ionization_table, load_ssp_cube


models.set_ext_law('CCM')

stats.set_walker_size(128)
stats.set_chain_size(5000)
stats.set_initial_positions([7., -2.1, 0.3, -2.])

prior_lowbounds = [6.0, np.log10(0.008), 0.01, -20.]
prior_highbounds = [7.5, np.log10(0.03), 1.0, 1.0]

fname = '/home/vital/PycharmProjects/SESAMME/docs/example_data/M83-8_preprocess.txt'
specfile = Table.read(fname, format='ascii')
wl = specfile['WL']
flux = specfile['FLUX']
flux_err = specfile['ERROR']

# But rescaling to luminosity density units (L_Sun A-1) allows SESAMME to infer a stellar mass for the cluster
lum = flux * 4*np.pi * (4.8*u.Mpc.to(u.cm))**2 / 3.83e33
lum_err = flux_err * 4*np.pi * (4.8*u.Mpc.to(u.cm))**2 / 3.83e33

windowlist = np.array([[np.min(wl), 1133], [1172, 1176.5],  [1188, 1202], [1203.8, 1222],
                       [1257, 1262], [1299, 1305], [1331, 1336.2], [1276, 1286], [1454, 1456], [1465, 1469],
                       [1390.5, 1393.5], [1399.5, 1403.3],
                       [1523.5, 1529], [1543, 1550.], [1608.5, 1621], [1656, 1659], [1666, 1674], [1795, np.max(wl)] ])

mask = models.get_mask(windowlist, wl)

# Original fitting
filename = "M83_test_v3.h5"
runname = 'Default_v3'
modelcube = load_ssp_cube("/home/vital/Astrodata/BPASS_v2.3/Demo_Met_Table.fits")
iontable = load_ionization_table("/home/vital/Astrodata/BPASS_v2.3/Demo_Q_Table.txt")

stats.run_sesamme(filename, runname, wl, lum, lum_err, modelcube, iontable, mask, True)

print(f'Loading {filename} for run {runname}')
reader = emcee.backends.HDFBackend(filename, name=runname)
samples = reader.get_chain()
flat_samples = reader.get_chain(discard=40, thin=10, flat=True)