# -*- coding: utf-8 -*-

"""MERGE HEADERS PACKAGE

This package contains the module(s) for ``merge_headers``.

:Author: Axel Guinot

:Parent module: ``split_exp_runner``

:Input: numpy binary files (``.npy``) with single-exposure header information

:Output: single SQL file with combined header information

Description
===========

This pipeline module merges the WCS information (image transformation and
distortions, computed during astrometrical calibration) for each CCD of
single-exposure images. The merged information is saved as a single SQL file.

Module-specific config file entries
===================================

OUTPUT_PATH : str, optional
    overrides the module output dir under ``[FILE]:OUTPUT_DIR``
"""

__all__ = ['merge_headers']
