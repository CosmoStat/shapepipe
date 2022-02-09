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
specificed by the config entry ``PREFIX``.

Masked pixels of different mask types are indicated by integers, which
conveniently are powers of two such that they can be combined bit-wise.

To mask bright stars, this module either creates a star catalogue from the
online
`guide star catalogue <https://heasarc.gsfc.nasa.gov/W3Browse/all/gsc.html>`_
database relevant to the the footprint. This is done by calling a CDs
(Centre de DonneÃ©s astronomique de Strasbourg)
`client <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_
program. Note that this requires online access,
which in some cases is not granted on compute nodes of a cluster. In this case,
set the config flag ``USE_EXT_STAR`` to False. Alternatively, a star catalogue
can be created before running this module via the script ``create_star_cat``.
During the processing of this module, this star catalogue is read from disk,
with ``USE_SET_STAR = True``.

The masking is done with the software ``WeightWatcher``, which is installed by ``ShapePipe``
by default.

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
PREFIX : str, optional, default=''
    prefix to be appended to output file name ``flag``;
    helps to distinguish the file patterns of newly created and external
    mask files

Mask config file
================

An additional configuration file is used by the mask module, its path is ``MASK_CONFIG_PATH``
in the module config section, see above. The following describes the config file sections
and their entries.

[PROGRAM_PATH]
--------------

WW_PATH : str
    path to the ``WeightWatcher`` executable
WW_CONFIG_FILE : str
    path to the ``WeightWatcher`` configuration file
CDSCLIENT_PATH : str, optional
    path to cds client executable; required if ``USE_EXT_STAR`` is ``False``

[BORDER_PARAMETERS]
-------------------

BORDER_MAKE : bool
    create mask around borders if True
BORDER_WIDTH : int
    width of border mask in pixels
BORDER_FLAG_VALUE : int (power of 2)
    border mask pixel value

[HALO_PARAMETERS]

HALO_MAKE : bool
    create mask for halos of bright stars if True
HALO_MASKMODEL_PATH : str
    .path to halo mask geometry (``.reg`` file)
HALO_MAG_LIM : float
    faint stellar magnitude limit for halo mask
HALO_SCALE_FACTOR : float
    factor to scale between magnitude (relative to pivot) and halo mask sise
HALO_MAG_PIVOT : float
    pivot stellar magnitude
HALO_FLAG_VALUE : int (power of 2)
    halo mask pixel value
HALO_REG_FILE : str
    output halo mask ``.reg`` file

[SPIKE_PARAMETERS]
------------------

SPIKE_MAKE : bool
    create mask for diffraction spikes of bright stars if True
SPIKE_MASKMODEL_PATH : str
    path to diffraction spike geometry (``.reg`` file)
SPIKE_MAG_LIM :
    faint stellar magnitude limit for spike mask
SPIKE_SCALE_FACTOR : float
    factor to scale between magnitude (relative to pivot) and spike mask size
SPIKE_MAG_PIVOT : float
    pivot stellar magnitude
SPIKE_FLAG_VALUE : int (power of two)
    diffraction spike pixel value
SPIKE_REG_FILE : str
    output spike mask ``.reg`` file

[MESSIER_PARAMETERS]
--------------------

MESSIER_MAKE : bool
    create mask around Messier objects if True
MESSIER_CAT_PATH : str
    path to Messier catalogue
MESSIER_SIZE_PLUS : float
   fraction to increase Messier mask
MESSIER_FLAG_VALUE : int (power of 2)
    Messier mask pixel value

[MD_PARAMETERS]
---------------

MD_MAKE : bool
    account for missing data (= zero-valued pixels) if True
MD_THRESH_FLAG : float
    threshold; if relative number of missing data is larger than this threshold,
    image is marked as flagged
MD_THRESH_REMOVE : float
    threshold; if relative number of missing data is larger than this threshold,
    image is marked for removal
MD_REMOVE : bool
    image is removed if marked for removal 

[OTHER]
-------

TEMP_DIRECTORY : str
    path to temporary dictionary
KEEP_INDIVIDUAL_MASK : bool
    keep individual masks in addition to merged mask file
KEEP_REG_FILE : bool
    keep .reg mask file

"""

__all__ = ['mask']
