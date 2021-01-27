# -*- coding: utf-8 -*-

"""MATCH EXTERNAL RUNNER

This module pastes different (SExtractor) catalogs of objects with identical IDs.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Axel Guinot

:Date: 10/2020

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

import pandas as pd


def get_data_ra_dec(path, hdu_no, col_ra, col_dec):

    cat = io.FITSCatalog(path)
    cat.open()
    data = cat.get_data(hdu_no)
    ra = data[col_ra]
    dec = data[col_dec]

    return data, ra, dec

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
    hdu_no : list of int
        (internal) catalog hdu number
    external_cat_path : string
        external catalog path
    external_col_match : list of strings
        external data column name(s) for matching
    external_col_match : list of strings, optional, default=None
        column name(s) to copy into matched output catalog
    hdu_no : int, optional, default=1
        external catalog hdu number
    """

    def __init__(self, input_file_list, output_path, w_log,
                 tolerance, col_match, hdu_no,
                 external_cat_path, external_col_match,
                 external_col_copy,
                 external_hdu_no=1):

        self._input_file_list = input_file_list
        self._output_path = output_path
        self._w_log = w_log

        self._tolerance = tolerance

        self._col_match = col_match
        self._hdu_no = hdu_no

        self._external_cat_path = external_cat_path
        self._external_col_match = external_col_match
        self._external_col_copy = external_col_copy
        self._external_hdu_no = external_hdu_no

    def process(self):

        external_data, external_ra, external_dec = get_data_ra_dec(
            self._external_cat_path, self._external_hdu_no,
            self._external_col_match[0], self._external_col_match[1])
        external_coord = SkyCoord(ra=external_ra, dec=external_dec, unit='deg')

        # (loop over inputs)
        for i in range(len(self._input_file_list)):
            data, ra, dec = get_data_ra_dec(self._input_file_list[i], self._hdu_no[i],
                                            self._col_match[0], self._col_match[1])
            coord = SkyCoord(ra=ra, dec=dec, unit='deg')
            # Todo: cut duplicates

            # Match objects in external cat to cat
            idx, d2d, d3d = match_coordinates_sky(coord, external_coord, nthneighbor=1)

            # Find close neighbours
            isdup = d2d < self._tolerance

            if not any(isdup==True):
                self._w_log.info('No match for {} with distance < {} arcsec found, no output created'
                                 ''.format(self._input_file_list[i], self._tolerance))
            else:

                # df_d2d to plot distribution of distances
                #df_to_match = pd.DataFrame(data={'isdup': isdup, 'idx': idx, 'd2d': d2d})
                #df_d2d = df_to_match.copy()
                #df_d2d['RA'], df_d2d['DEC'] = ra_external, dec_external
                #df_to_match.query("isdup == True", inplace=True)
                #df_d2d.query("isdup == True", inplace=True)
                #df_to_match.drop(columns=['isdup'], inplace=True)
                #df_d2d.drop(columns=['isdup'], inplace=True)

                idx_sub = np.array([(i,ide) for (i,ide) in enumerate(idx) if isdup[i] == True])[:,1]

                z_spec_sub = []
                for i in idx_sub:
                    z_spec_sub.append(external_data[self._external_col_copy][i])
                z_spec_sub = np.array(z_spec_sub)

                # Create new, matched cat
                t = Table(data)
                df_to_cut = t.to_pandas()
                df_matched = df_to_cut.copy()
                df_matched['isdup'] = isdup
                df_matched.query("isdup == True", inplace=True)
                df_matched.drop(columns=['isdup'], inplace=True)
                df_matched[self._external_col_copy] = z_spec_sub

                # Write output
                t = Table.from_pandas(df_matched)
                t.write(self._output_path, overwrite=True)

                # TODO: Compute stats


@module_runner(version='1.0',
               input_module='sextractor_runner',
               file_pattern='tile_sexcat',
               file_ext='.fits',
               depends=['numpy', 'astropy', 'pandas'],
               run_method='parallel')
def match_external_runner(input_file_list, run_dirs, file_number_string,
                          config, w_log):


    # Processing
    tmp = config.getfloat('MATCH_EXTERNAL_RUNNER', 'TOLERANCE')
    tolerance = tmp * u.arcsec

    # Internal data
    col_match = config.getlist('MATCH_EXTERNAL_RUNNER', 'COL_MATCH')
    if config.has_option('MATCH_EXTERNAL_RUNNER', 'HDU'):
        tmp = config.getlist('MATCH_EXTERNAL_RUNNER', 'HDU')
        hdu_no = [int(i) for i in tmp]
        if len(hdu_no) != len(input_file_list):
            raise IndexError('Different lengths for input file list ({}) and'
                             'HDU ({})'
                             ''.format(len(input_file_list), len(hdu_no)))
    else:
        hdu_no = [2] * len(input_file_list)

    # External data
    external_cat_path = config.getexpanded('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_CAT_PATH')
    external_col_match = config.getlist('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_COL_MATCH')

    # TODO: optional, list, 'all'
    external_col_copy = config.get('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_COL_COPY')

    if config.has_option('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_HDU'):
        external_hdu_no = config.getint('MATCH_EXTERNAL_RUNNER', 'EXTERNAL_HDU')
    else:
        external_hdu_no = 1

    # Output
    if config.has_option('MATCH_EXTERNAL_RUNNER', 'OUTPUT_FILE_PATTERN'):
        output_file_pattern = config.get('MATCH_EXTERNAL_RUNNER', 'OUTPUT_FILE_PATTERN')
    else:
        output_file_pattern = 'cat_matched'

    file_ext = 'fits'

    output_path = '{}/{}{}.{}'.format(run_dirs['output'],
                                      output_file_pattern,
                                      file_number_string,
                                      file_ext)

    inst = MatchCats(input_file_list, output_path, w_log,
                     tolerance, col_match, hdu_no,
                     external_cat_path, external_col_match,
                     external_col_copy,
                     external_hdu_no=external_hdu_no)

    inst.process()

    return None, None
