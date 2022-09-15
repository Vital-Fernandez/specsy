Specsy: The Spectra Synthesis library
=====================================


This library provides a set of tools for the analysis of astronomical spectra. This is the alpha-release.

Most of its functions assume an input with the `LiMe <https://lime-stable.readthedocs.io/en/latest/>`_ formatting.

These are the features currently available:

* Calculation of the logarithmic extinction calculation coefficient :math:`c(H\beta)`.
* Calculation of the electron density from the :math:`[SII]6716,6731\AA` doublet via a Monte-Carlo sampling
* Calculation of sulphur abundance using the calibration from `Díaz et al (2022) <https://doi.org/10.1093/mnras/stac387>`_

.. toctree::
   :maxdepth: 1
   :caption: Documentation:
   :name: doc-tree

   documentation/installation
   documentation/api
