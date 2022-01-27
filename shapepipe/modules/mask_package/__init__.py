"""MASK MODULE.

This package contains the module for ``mask``.

:Author: Axel Guinot

:Parent module: ``split_exp_runner`` or None

:Input: single-exposure single-CCD image, weight file, flag file (optional), star catalogue (optional)

:Output: single-exposure single-CCD flag files

Description
===========

This module creates masks for bright stars, diffraction spikes, Messier
objects, borders, and other artifacts. If a flag file is given as input,
for example from pre-processing, the mask that is created by this module
is joined with the mask from this external flag file. In this case
the config flag ``USE_EXT_FLAG`` needs to be set to ``True``. To distinguish
the newly created output flag file from the input ones, a prefix can added as
specificed by the config entry ``SUFFIX``.

Masked pixels of different mask types are indicated by integers, which
conveniently are powers of two such that they can be combined bit-wise.

To mask bright stars, this module either creates a star catalogue from the
online [guide star
catalogue](https://heasarc.gsfc.nasa.gov/W3Browse/all/gsc.html) database
relevant to the the footprint. Note that this requires online access, which in
some cases is not granted on compute nodes of a cluster. In this case, set the
config flag ``USE_EXT_STAR`` to False. Alternatively, a star catalogue can be
created before running this module via the script ``create_star_cat``. During
the processing of this module, this star catalogue is read from disk, with
``USE_SET_STAR = True``.

Module-specific config file entries
===================================

USE_EXT_FLAG : bool
    use external flag file to join with the mask created here;
    if ``True`` flag file needs to be given on input
USE_EXT_STAR : bool
    read external star catalogue instead of creating one during the
    call of this module;
    if ``True`` star catalogue file needs to be given on input
MASK_CONFIG_PATH : str
    path to mask config file
HDU : int, optional, default=0
    HDU of external flag FITS file
SUFFIX : str, optional, default=''
    prefix to be appended to output file name ``flag``;
    helps to distinguish the file patterns of newly created and external
    mask files

"""

__all__ = ['mask']
