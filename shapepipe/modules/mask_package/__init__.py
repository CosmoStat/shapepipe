"""MASK MODULE.

This package contains the module for ``mask``.

:Author: Axel Guinot

:Parent module: ``split_exp_runner`` or None

:Input: Single-exposure single-CCD image, weight file, flag file (optional),
        and star catalogue (optional)

:Output: Single-exposure single-CCD flag files

Description
===========

This module creates masks for bright stars, diffraction spikes, deep sky
objects (from the Messier and NGC catalogues), borders, and other artifacts. If
a flag file is given as input, for example from pre-processing, the mask that
is created by this module is joined with the mask from this external flag file.
In this case the config flag ``USE_EXT_FLAG`` needs to be set to ``True``. To
distinguish the newly created output flag file from the input ones, a prefix
can added as specificed by the config entry ``PREFIX``.

An NGC catalogue with positions, sizes, and types is provided with
``shapepipe``,
`source <http://www.klima-luft.de/steinicke/ngcic/rev2000/Explan.htm>`_.

Masked pixels of different mask types are indicated by integers, which
conveniently are powers of two such that they can be combined bit-wise.

To mask bright stars, this module either creates a star catalogue from the
online
`guide star catalogue <https://heasarc.gsfc.nasa.gov/W3Browse/all/gsc.html>`_
database relevant to the the footprint. This is done by calling a CDs
(Centre de Données astronomique de Strasbourg)
`client program <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_.
Note that this requires online access,
which in some cases is not granted on compute nodes of a cluster. In this case,
set the config flag ``USE_EXT_STAR = False``. Alternatively, a star
catalogue can be created before running this module via the script
``create_star_cat``. During the processing of this module, this star catalogue
is read from disk, with ``USE_SET_STAR = True``.

The masking is done with the software ``WeightWatcher`` :cite:`marmo:08`,
which is installed by ``ShapePipe`` by default.

Module-specific config file entries
===================================

USE_EXT_FLAG : bool
    Use external flag file to join with the mask created here;
    if ``True`` flag file needs to be given on input
USE_EXT_STAR : bool
    Read external star catalogue instead of creating one during the
    call of this module;
    if ``True`` star catalogue file needs to be given on input
MASK_CONFIG_PATH : str
    Path to mask config file
HDU : int, optional
    HDU of external flag FITS file; the default value is ``0``
PREFIX : str, optional
    Prefix to be appended to output file name ``flag``;
    helps to distinguish the file patterns of newly created and external
    mask files
CHECK_EXISTING_DIR : str, optional
    If given, search this directory for existing mask files; the
    corresponding images will then not be processed

Mask config file
================

An additional configuration file is used by the mask module, its path is
``MASK_CONFIG_PATH`` in the module config section, see above. The following
describes the config file sections and their entries.

[PROGRAM_PATH]
--------------

WW_PATH : str, optional
    Full path to the WeightWatcher executable (``ww``) on the system ; if not
    set the version controlled WeightWatcher installation in the ShapePipe
    environment will be used
WW_CONFIG_FILE : str
    Path to the WeightWatcher configuration file
CDSCLIENT_PATH : str, optional
    Path to CDS client executable; required if ``USE_EXT_STAR = False``

[BORDER_PARAMETERS]
-------------------

BORDER_MAKE : bool
    Create mask around borders if ``True``
BORDER_WIDTH : int
    Width of border mask in pixels
BORDER_FLAG_VALUE : int
    Border mask pixel value, power of 2

[HALO_PARAMETERS]
-----------------

HALO_MAKE : bool
    Create mask for halos of bright stars if ``True``
HALO_MASKMODEL_PATH : str
    Path to halo mask geometry (``.reg`` file)
HALO_MAG_LIM : float
    Faint stellar magnitude limit for halo mask
HALO_SCALE_FACTOR : float
    Factor to scale between magnitude (relative to pivot) and halo mask size
HALO_MAG_PIVOT : float
    Pivot stellar magnitude
HALO_FLAG_VALUE : int
    Halo mask pixel value, power of 2
HALO_REG_FILE : str
    Output halo mask ``.reg`` file

[SPIKE_PARAMETERS]
------------------

SPIKE_MAKE : bool
    Create mask for diffraction spikes of bright stars if ``True``
SPIKE_MASKMODEL_PATH : str
    Path to diffraction spike geometry (``.reg`` file)
SPIKE_MAG_LIM :
    Faint stellar magnitude limit for spike mask
SPIKE_SCALE_FACTOR : float
    Factor to scale between magnitude (relative to pivot) and spike mask size
SPIKE_MAG_PIVOT : float
    Pivot stellar magnitude
SPIKE_FLAG_VALUE : int
    Diffraction spike pixel value, power of two
SPIKE_REG_FILE : str
    Output spike mask ``.reg`` file

[MESSIER_PARAMETERS]
--------------------

MESSIER_MAKE : bool
    Create mask around Messier objects if ``True``
MESSIER_CAT_PATH : str
    Path to Messier catalogue
MESSIER_SIZE_PLUS : float
    Fraction to increase Messier mask
MESSIER_FLAG_VALUE : int
    Messier mask pixel value, power of 2

[NGC_PARAMETERS]
--------------------

NGC_MAKE : bool
    Create mask around NGC objects if ``True``
NGC_CAT_PATH : str
    Path to NGC catalogue
NGC_SIZE_PLUS : float
    Fraction to increase NGC mask
NGC_FLAG_VALUE : int
    NGC mask pixel value, power of 2

[MD_PARAMETERS]
---------------

MD_MAKE : bool
    Account for missing data (zero-valued pixels) if ``True``
MD_THRESH_FLAG : float
    Threshold; if relative number of missing data is larger than this
    threshold, image is marked as flagged
MD_THRESH_REMOVE : float
    Threshold; if relative number of missing data is larger than this
    threshold, image is marked for removal
MD_REMOVE : bool
    Image is removed if marked for removal

[OTHER]
-------

TEMP_DIRECTORY : str
    Path to temporary dictionary
KEEP_INDIVIDUAL_MASK : bool
    Keep individual masks in addition to merged mask file
KEEP_REG_FILE : bool
    Keep ``.reg`` mask file

"""

__all__ = ['mask']
