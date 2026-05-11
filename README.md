# Specsy

<p align="center">
  <img src="https://github.com/Vital-Fernandez/specsy/blob/7e35568f6d154486f5603e94fe39dd08e5e54834/src/specsy/resources/images/Specsy_logo_transparent_dark.PNG" alt="Specsy Logo" width="300"/>
</p>

A Python library for the analysis of astronomical spectra. Specsy includes a Bayesian sampler for the
direct method parameter space, tools to fit photoionization model grids, and utilities for the analysis
of stellar and gas continua.

> **Note:** This package is currently in an alpha release. The preliminary documentation can be found at [ReadTheDocs](https://specsy.readthedocs.io/).

## Installation

Install directly from [PyPI](https://pypi.org/project/specsy/):

```bash
pip install specsy
```

For the recommended conda environment with PyMC sampler backends:

```bash
conda create -c conda-forge -n specsy python=3.13 nutpie pymc numba numpyro blackjax
conda activate specsy
pip install specsy
```

To upgrade to the latest version:

```bash
pip install --upgrade specsy
```

## Development

SpecSy is currently in an alpha release. Please check the [GitHub repository](https://github.com/Vital-Fernandez/specsy)
for the latest version or to report any issues.

**Author:** Vital Fernández — [vgf@stsci.edu](mailto:vgf@stsci.edu)