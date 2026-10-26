"""MERGE HEADERS PACKAGE.

This package contains the module for ``merge_headers``.

:Author: Axel Guinot

:Parent module: ``find_exposures_runner``

:Input: An ``exp_numbers`` text file listing the tile's exposures; the
    per-exposure header files (``.npy``) are found under ``EXP_BASE_DIR``

:Output: Single SQL file with combined header information

Description
===========

This pipeline module merges the WCS information (image transformation and
distortions, computed during astrometrical calibration) for each CCD of
single-exposure images. The merged information is saved as a single SQL file.

Module-specific config file entries
===================================

EXP_BASE_DIR : str
    Root of the per-exposure work directories; the header ``.npy`` files are
    collected from ``<EXP_BASE_DIR>/<prefix>/<base>/`` for every exposure
    listed in the input ``exp_numbers`` file

"""

__all__ = ["merge_headers"]
