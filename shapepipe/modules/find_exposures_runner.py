# -*- coding: utf-8 -*-

"""FIND_EXPOSURES RUNNER

This module runs find_exposures: Identify exposures that are used in selected
tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: January 2019

"""

import re

import astropy.io.fits as fits

import shapepipe.pipeline.file_io as io
from shapepipe.modules.module_decorator import module_runner


class find_exposures():

    """Class for identifying exposures that are used for a given tile.

    Parameters
    ----------
    img_tile_path: string
        path to tile image file
    output_path: string
        output file path
    w_log: log file class
        log file

    Returns
    -------
    None
    """

    def __init__(self, img_tile_path, output_path, w_log):

        self._img_tile_path = img_tile_path
        self._output_path = output_path
        self._w_log = w_log

    def process(self):

        exp_list_uniq = self.get_exposure_list()

        f_out = open(self._output_path, 'w')
        if len(exp_list_uniq) > 0:
            for exp in exp_list_uniq:
                print(exp, file=f_out)
        f_out.close()

    def get_exposure_list(self):
        """Return list of exposure file used for the tile in process, from tiles FITS header

        Parameters
        ----------

        Returns
        -------
        exp_list_uniq: list of string
            list of exposure basenames
        """

        try:
            hdu = fits.open(self._img_tile_path)
            hist = hdu[0].header['HISTORY']

        except Exception:
            self._w_log.info('Error while reading tile image FITS file '
                             '{}, continuing...'.format(self._img_tile_path))

        exp_list = []

        # Get exposure file names
        for h in hist:
            temp = h.split(' ')

            pattern = r'(.*)\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise IndexError('re match \'{}\' failed for filename '
                                 '\'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)

            exp_list.append(exp_name)

        exp_list_uniq = list(set(exp_list))

        self._w_log.info('Found {} exposures used in tile'.format(
                         len(exp_list_uniq)))
        self._w_log.info('{} duplicates were removed'.format(
                         len(exp_list) - len(exp_list_uniq)))

        return exp_list_uniq


@module_runner(version='1.0', file_pattern=['image'], file_ext='.fits',
               depends=['numpy', 'astropy', 'sip_tpv'], numbering_scheme='_0')
def find_exposures_runner(input_file_list, run_dirs, file_number_string,
                          config, w_log):

    input_file_name = input_file_list[0]

    output_path = '{}/exp_numbers{}.txt'.format(run_dirs['output'], file_number_string)

    inst = find_exposures(input_file_name,
                          output_path,
                          w_log)
    inst.process()

    return None, None
