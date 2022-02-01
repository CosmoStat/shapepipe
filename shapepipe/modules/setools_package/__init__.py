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
    path to setools configuration file

Setools configuration file
==========================

This file can contain an arbitrary number of section. Each section is initiated by a line
``[TYPE:name]``. The following mask types ``TYPE`` are valid:
    - ``MASK``: define a selection or mask
    - ``RAND_SPLIT``: split an input into random subsamples
    - ``PLOT``: create a plot
    - ``STATS``: output summary statistics
``name`` is the file name, in which the result of the section is saved.
This is a FITS file for ``MASK``, .png file for ``PLOT``, and an ASCII file for ``STATS``.
For ``TYPE = MASK`` and ``RAND_SPLIT, ``name`` is also the reference, by which the defined
sample can be addressed within the configuration file.

"""

__all__ = ['setools']
