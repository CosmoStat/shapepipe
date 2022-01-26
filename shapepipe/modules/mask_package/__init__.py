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
is joined with the mask from this external flag file.

Masked pixels of different mask types are indicated by integers, which
conveniently are powers of two such that they can be combined bit-wise.

Exposures, unlike tile images, come with external flag files on input. This is
specified by the key USE_EXT_FLAG. To distinguish the newly created output flag
file from the input ones, a suffix is added:

To mask bright stars, this module either creates a star catalogue from the online
[guide star catalogue](https://heasarc.gsfc.nasa.gov/W3Browse/all/gsc.html) data
base relevant to the the footprint. Note that this requires online access, which in
some cases is not granted on compute nodes of a cluster.
Alternatively, the star catalogue can be read from disk.

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
    prefix to be appended to output file name 'flag';
    helps to distinguish the file patterns of newly created and external
    mask files

"""

__all__ = ['mask']
