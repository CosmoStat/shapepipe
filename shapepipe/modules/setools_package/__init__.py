"""SETOOLS PACKAGE.

This package contains the module(s) for ``setools``.

:Author: Axel Guinot

:Parent module: ``sextractor_runner``

:Input: ``SExtractor`` catalogue

:Output: ``SExtractor`` catalogue

Description
===========

This module defines samples via selection criteria, or masks. A variety of statistical
mathematical functions are implemented, that can be applied to the input catalogue.

The ``setools`` module can create random sub-sample output catalogues, for example
for training and validation purposes.

In addition, the module can produce user-defined plots, and compute summary
statistics of the selected population(s). The figures include scatter plots and histograms,
the summaries are number, mean, mode, extrema, and standard deviation.


Module-specific config file entries
===================================
SETOOLS_CONFIG_PATH : str
    path to module specific configuration file

"""

__all__ = ['setools']
