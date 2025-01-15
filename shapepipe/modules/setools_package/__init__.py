"""SETOOLS PACKAGE.

This package contains the module(s) for ``setools``.

:Author: Axel Guinot

:Parent module: ``sextractor_runner``

:Input: SExtractor catalogue

:Output: SExtractor catalogue

Description
===========

This module defines samples via selection criteria, or masks. A variety of
statistical mathematical functions are implemented, that can be applied to the
input catalogue.

The ``setools`` module can create random sub-sample output catalogues, for
example for training and validation purposes.

In addition, the module can produce user-defined plots, and compute summary
statistics of the selected population(s). The figures include scatter plots and
histograms, the summaries are number, mean, mode, extrema, and standard
deviation.


Module-specific config file entries
===================================

SETOOLS_CONFIG_PATH : str
    Path to setools configuration file

Setools configuration file
==========================

This file can contain an arbitrary number of sections. Each section is
initiated by a line ``[TYPE:name]``. The following mask types ``TYPE``
are valid.

- ``MASK``: define a selection or mask
- ``RAND_SPLIT``: split an input into random subsamples
- ``PLOT``: create a plot
- ``STATS``: output summary statistics

``name`` is the file name, in which the result of the section is saved.
This is a FITS file for ``MASK``, ``.png`` file for ``PLOT``, and an ASCII file
for ``STATS``. For ``TYPE = MASK`` and ``RAND_SPLIT``, ``name`` is also the
reference, by which the defined sample can be addressed within the
configuration file.

The following is an example of a star selection using two ``MASK`` sections.
First, a pre-selection is defined:

.. code-block:: ini

   [MASK:preselect]
   MAG_AUTO > 0
   MAG_AUTO < 21
   FWHM_IMAGE > 0.3 / 0.187
   FWHM_IMAGE < 1.5 / 0.187
   FLAGS == 0
   IMAFLAGS_ISO == 0
   NO_SAVE

This selects objects within a magnitude (``MAG_AUTO``) and size
(``FWHM_IMAGE``) ranges. The size limits of 0.3" and 1.5" are transformed from
arcseconds to pixels. Additional flag criteria for ``FLAGS`` AND
``IMAFLAGS_ISO`` are specified. The keyword ``NO_SAVE`` indicates that this
selection is not to be saved to disk.

Next, the final star selection is defined:

.. code-block:: ini

   [MASK:star_selection]
   MAG_AUTO > 18.
   MAG_AUTO < 22.
   FWHM_IMAGE <= mode(FWHM_IMAGE{preselect}) + 0.2
   FWHM_IMAGE >= mode(FWHM_IMAGE{preselect}) - 0.2
   FLAGS == 0
   IMAFLAGS_ISO == 0

The size range is now refined using the mode of preselected objects. The
preselection removes outliers before the mode computation.

An example of the definition of random subsamples is as follows:

.. code-block:: ini

   [RAND_SPLIT:star_split]
   RATIO = 20
   MASK = star_selection

This uses as input the ``star_selection`` sample (given by the ``MASK`` entry).
It creates two random subsamples: one with 20% and one with 80% (100 - RATIO)
percent of input objects. The output file base names are given by the mask name
and the ratios; in this case two files ``star_split_ratio_20-<num>.fits`` and
``star_split_ratio_80-<num>.fits`` are created, where ``<num>`` is the pipeline
number of the (single-exposure single-HDU) image that is processed.

Entries in the ``MASK`` and ``STATS`` sections, lines can contain basic
relational operations, and functions.

The operands can be

- input column names, or ``SExtractor`` keywords,
- or names defined by a previous mask.

The following are illustrations of two plots. First, another mask
is defined.

.. code-block:: ini

   [MASK:fwhm_mag_cut]
   FWHM_IMAGE > 0
   FWHM_IMAGE < 40
   MAG_AUTO < 35
   FLAGS == 0
   IMAFLAGS_ISO == 0
   NO_SAVE

   [PLOT:size_mag]
   TYPE = plot
   FORMAT = png
   X_1 = FWHM_IMAGE{fwhm_mag_cut}
   Y_1 = MAG_AUTO{fwhm_mag_cut}
   X_2 = FWHM_IMAGE{star_selection}
   Y_2 = MAG_AUTO{star_selection}
   MARKER_1 = +
   MARKER_2 = .
   MARKERSIZE_1 = 3
   MARKERSIZE_2 = 3
   LABEL_1 = All
   LABEL_2 = "Stars, mean FWHM: @mean(FWHM_IMAGE{star_selection})*0.187@arcsec"
   TITLE = "Stellar locus"
   XLABEL = "FWHM (pix)"
   YLABEL = Mag

   [PLOT:hist_mag_stars]
   TYPE = hist
   FORMAT = png
   Y = MAG_AUTO{star_selection}
   BIN = 20
   LABEL = "stars"
   XLABEL = "Magnitude"
   YLABEL = "Number"
   TITLE = "Magnitude of stars"

Possible values for ``TYPE`` are ``plot``, ``histogram``, ``scatter``.
More than one sample can be plotted by adding ``_1``, ``_2``, ... to all of the
above given keywords.

Keywords for all plot types:

Y : str
    Functions of samples to plot on y-axis
XLABEL, YLABEL : str, optional
    Labels for x- and y-axes, respectively
LABEL : str, optional
    Plot label
TITLE : str, optional
    Plot title
ALPHA : float
    Alpha (transparency) parameter
COLOR : str, optional
    Plot color; the default is the standard matplotlib color sequence
FORMAT : str, optional
    Output file format; the default value is ``png``

Keywords for standard and scatter plots (``TYPE`` in ``plot``, ``scatter``):

X : str
    Functions of samples to plot on x-axis
MARKER : str, optional
    Point marker symbol; the default value is ``+``

Additional keywords for standard line/points plot (``TYPE = plot``):

MARKERSIZE : float, optional
    Marker size; the default value is ``1``
LINE : str, optional
    Line type

Additional keywords for histograms (``TYPE = hist``):

BIN : int, optional
    Number of bins; the default value is ``50``
HTYPE : str, optional
    Histogram type; the default value is ``bar``
LOG : bool, optional
    Logarithmic (linear) bins if ``True`` (``False``); the defaul value is
    ``False``

Additional keywords for scatter plots (``TYPE = scatter``):

MARKER : str
    Point marker
SCATTER : str, optional
    Function of some input sample acting as third variable that is color-coded
    in the point colors.

"""

__all__ = ["setools"]
