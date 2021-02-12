# -*- coding: utf-8 -*-

"""MATCH EXTERNAL RUNNER

This module matches an external catalogue to a ShapePipe (SExtractor) catalog

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Xavier Jimenez

:Date: 01/2021

:Package: ShapePipe

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io

import os
import re

import numpy as np

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky


def get_cat(path):

    cat = io.FITSCatalog(path)
    cat.open()

    return cat


def get_data(path, hdu_no):

    cat = get_cat(path)
    data = cat.get_data(hdu_no)
    col_names = cat.get_col_names(hdu_no=hdu_no)
    ext_names = cat.get_ext_name()
    cat.close()

    return data, col_names, ext_names


def get_ra_dec(data, col_ra, col_dec):

    ra = data[col_ra]
    dec = data[col_dec]

    return ra, dec


class MatchCats(object):
    """MatchCat initialisation.

    Parameters
    ----------
    input_file_list : list of strings
        list of input catalogue paths to be pasted.
    output_path : string
        output file path of pasted catalog
    w_log :
        log file
    tolerance : Quantity
        tolerance, with units.arcsec
    col_match : list of strings, optional, default=None
        (internal data) column name(s) to copy into matched output catalog
    hdu_no : int
        (internal) catalog hdu number
    mode : string
        mode, 'CLASSIC' or 'MULTI-EPOCH'
    external_cat_path : string
        external catalog path
    external_col_match : list of strings
        external data column name(s) for matching
    external_col_copy : list of strings
        column name(s) to copy into matched output catalog
    external_hdu_no : int, optional, default=1
        external catalog hdu number
    mark_non_matched : float, optional, None
        if not None, output not only matched but all objects, and mark
        non-matched objects with this value
    """

    def __init__(self, input_file_list, output_path, w_log,
                 tolerance, col_match, hdu_no, mode,
                 external_cat_path, external_col_match,
                 external_col_copy,
                 external_hdu_no=1,
                 mark_non_matched=None):

        self._input_file_list = input_file_list
        self._output_path = output_path
        self._w_log = w_log

        self._tolerance = tolerance

        self._col_match = col_match
        self._hdu_no = hdu_no
        self._mode = mode

        self._external_cat_path = external_cat_path
        self._external_col_match = external_col_match
        self._external_col_copy = external_col_copy
        self._external_hdu_no = external_hdu_no

        self._mark_non_matched = mark_non_matched

    def process(self):

        # Load external and internal data
        external_data, dummy1, dummy2 = get_data(self._external_cat_path, self._external_hdu_no)
        external_ra, external_dec = get_ra_dec(external_data,
                                               self._external_col_match[0],
                                               self._external_col_match[1])
        external_coord = SkyCoord(ra=external_ra, dec=external_dec, unit='deg')

        data, col_names, ext_names = get_data(self._input_file_list[0], self._hdu_no)
        ra, dec = get_ra_dec(data, self._col_match[0], self._col_match[1])
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')

        # Todo: cut duplicates

        # Match objects in external cat to internal cat. idx=indices to external object
        # for each object in internal cat e.g. external_coord[idx[0]] is the match for
        # coord[0].
        idx, d2d, d3d = match_coordinates_sky(coord, external_coord, nthneighbor=1)

        # Find close neighbours, idx_close is True for all close matches
        idx_close = d2d < self._tolerance

        if not any(idx_close):
            self._w_log.info('No match for {} with distance < {} arcsec found, no output created'
                             ''.format(self._input_file_list[0], self._tolerance))
        else:

            # Get indices in internal and external catalogues of pair-wise matches
            w = np.array([(i, ide) for (i, ide) in enumerate(idx) if idx_close[i]])
            id_sub = w[:, 0]
            id_ext_sub = w[:, 1]
            id_all = np.arange(len(idx))

            if self._mark_non_matched:
                # Output all objects
                id_data = id_all
                id_ext = idx
            else:
                # Output only matched objects
                id_data = id_sub
                id_ext = id_ext_sub

            self._w_log.info('{} objects matched out of {}'
                             ''.format(len(id_sub), len(idx)))

            # Copy matched objects from internal catalogue to output data
            matched = {}
            for col in col_names:
                matched[col] = data[col][id_data]

            # Copy columns from external catalogue to output data
            for col in self._external_col_copy:
                matched[col] = external_data[col][id_ext]
                if self._mark_non_matched:
                    for i, i_ext in enumerate(idx):
                        if not idx_close[i]:
                            matched[col][i] = self._mark_non_matched

            # Write FITS file
            out_cat = io.FITSCatalog(self._output_path, SEx_catalog=True,
                                     open_mode=io.BaseCatalog.OpenMode.ReadWrite)
            out_cat.save_as_fits(data=matched, ext_name='MATCHED', sex_cat_path=self._input_file_list[0])

            # Write all extensions if in multi-epoch mode
            if self._mode == 'MULTI-EPOCH':
                hdu_me_list = [i for i, name in enumerate(ext_names)
                               if 'EPOCH' in name]
                for hdu_me in hdu_me_list:
                    data_me, col_names_me, dummy = get_data(self._input_file_list[0], hdu_me)
                    matched_me = {}
                    for col_me in col_names_me:
                        matched_me[col_me] = data_me[col_me][id_data]
                    out_cat.save_as_fits(data=matched_me, ext_name=ext_names[hdu_me])

            # TODO: Compute stats


@module_runner(version='1.0',
               input_module='sextractor_runner',
               file_pattern='tile_sexcat',
               file_ext='.fits',
               depends=['numpy', 'astropy'],
               run_method='parallel')
def match_external_runner(input_file_list, run_dirs, file_number_string,
                          config, w_log):

    # Processing
    tmp = config.getfloat('MATCH_EXTERNAL_RUNNER', 'TOLERANCE')
    tolerance = tmp * u.arcsec

    # Internal data
    col_match = config.getlist('MATCH_EXTERNAL_RUNNER', 'COL_MATCH')
    if config.has_option('MATCH_EXTERNAL_RUNNER', 'HDU'):
        hdu_no = config.getint('MATCH_EXTERNAL_RUNNER', 'HDU')
    else:
        hdu_no = 2

    mode = config.get('MATCH_EXTERNAL_RUNNER', 'MODE')
    valid_modes = ['CLASSIC', 'MULTI-EPOCH']
    if mode not in valid_modes:
        raise ValueError('mode \'{}\' is invalid, must be one of {}'.format(mode, valid_modes))

    # External data
    external_cat_path = config.getexpanded('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_CAT_PATH')
    external_col_match = config.getlist('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_COL_MATCH')

    # TODO: optional or 'none', 'all'
    # Also TODO: change column name if already present in internal cat
    external_col_copy = config.getlist('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_COL_COPY')

    if config.has_option('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_HDU'):
        external_hdu_no = config.getint('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_HDU')
    else:
        external_hdu_no = 1

    # Output
    if config.has_option('MATCH_EXTERNAL_RUNNER', 'OUTPUT_FILE_PATTERN'):
        output_file_pattern = config.get('MATCH_EXTERNAL_RUNNER', 'OUTPUT_FILE_PATTERN')
    else:
        output_file_pattern = 'cat_matched'

    if config.has_option('MATCH_EXTERNAL_RUNNER', 'MARK_NON_MATCHED'):
        mark_non_matched = config.getfloat('MATCH_EXTERNAL_RUNNER', 'MARK_NON_MATCHED')
    else:
        mark_non_matched = None

    file_ext = 'fits'

    output_path = '{}/{}{}.{}'.format(run_dirs['output'],
                                      output_file_pattern,
                                      file_number_string,
                                      file_ext)

    inst = MatchCats(input_file_list, output_path, w_log,
                     tolerance, col_match, hdu_no, mode,
                     external_cat_path, external_col_match,
                     external_col_copy,
                     external_hdu_no=external_hdu_no,
                     mark_non_matched=mark_non_matched)

    inst.process()

    return None, None
