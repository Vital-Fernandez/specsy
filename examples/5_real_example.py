import numpyro

# Number of cores available
numpyro.set_host_device_count(8)

import jax

# Selecting CPU sampling with JAX
jax.config.update("jax_platform_name", "cpu")

from astropy.io import fits
from pathlib import Path
import lime
import specsy as sy

# lime.theme.set_style('dark')
#
# import aspect
# fname = Path('/home/vital/Downloads/espectros_regionesHII_UVES (2)/fc_1D_DIC2_437_NGC3576_5400s.fits')
# fname = Path('/home/vital/Downloads/espectros_regionesHII_UVES (2)/fc_1D_DIC2_437_NGC3576_5400s.fits')
# fname = Path('/home/vital/Downloads/espectros_regionesHII_UVES (2)/fc_1D_DIC2_860_NGC3576_5400s.fits')
# fname = Path('/home/vital/Downloads/espectros_regionesHII_UVES (2)/fc_1D_DIC2_860_S311_5400s.fits')
# fname = Path('/home/vital/Downloads/espectros_regionesHII_UVES (2)/fc_1D_DIC1_580_NGC3576_1800s.fits')
#
# folder = Path('/home/vital/Downloads/espectros_regionesHII_UVES (2)')
# fits_files = list(folder.glob("*.fits"))
#
# for i, fname in enumerate(fits_files):
#     if i == 2:
#         print(i, fname.name)
#         spec = lime.Spectrum.from_file(fname, instrument='ISIS', redshift=0)
#         bands = spec.retrieve.lines_frame(band_vsigma=30)
#         spec.plot.spectrum(bands=None, log_scale=False)

cfg_fname = './synthetic_spectrum_region_v0.toml'
cfg = lime.load_cfg(cfg_fname)

linesdf_fname = '/home/vital/PycharmProjects/lime/examples/doc_notebooks/0_resources/results/manga_lines_log_from_spectrum.txt'
lines_df = lime.load_frame(linesdf_fname)

# Emissivity file
emis_fname = f'./emissivity_grids_pyneb_1.1.30.nc'

obj = sy.Nebula.from_lines_frame(lines_df, cfg)
line_list = ['H1_4861A', 'H1_4340A', 'H1_6563A',
             'O2_3726A', 'O2_3729A',
             'O3_4363A', 'O3_4959A', 'O3_5007A',
             'N2_6548A', 'N2_6583A',
             'Ar3_7136A', 'Ar3_7751A', 'Ar4_4740A'
             'S2_6716A', 'S2_6731A',
             'S3_6312A', 'S3_9068A', 'S3_9530A']

obj.infer.direct_method.prepare_inputs(emissivity_source=emis_fname, line_list=line_list, prior_cfg=cfg['direct_method_priors'])

obj.infer.direct_method.run()

obj.infer.direct_method.plot_trace()





