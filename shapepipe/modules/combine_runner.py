# -*- coding: utf-8 -*-
  
"""COMBINE RUNNER

This module combines the output from previous runs.
Symbolic links to previous result files are created
in this module's output directory.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 10/2020

:Package: ShapePipe

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io

import os
import re

import re
import os

class Combine(object):
    """ Combine

    Combine results from previous runs.

    Parameters
    ----------
    input_file_list : list of strings
        input file paths
    output_dir : string
        output directory name
    w_log :
        log file handler

    """

    def __init__(self, input_file_list, output_dir, w_log):

        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._w_log = w_log

    def process(self):

        n_skipped = 0
        n_created = 0
        for item in self._input_file_list:
            input_path = item[0]
            input_base_name = os.path.basename(input_path)
            output_path = '{}/{}'.format(self._output_dir, input_base_name)
            if os.path.exists(output_path):
                n_skipped = n_skipped + 1
            else:
                os.symlink(input_path, output_path)
                n_created = n_created + 1

        self._w_log.info('Number of links created = {}'.format(n_created))
        self._w_log.info('Number of (duplicate) files skipped = {}'.format(n_skipped))
        self._w_log.info('Total number of input files = {}'
                         ''.format(len(self._input_file_list)))


@module_runner(version='1.0',
               run_method='serial')
def combine_runner(input_file_list, run_dirs, file_number_string,
                   config, w_log):

    combine = Combine(input_file_list, run_dirs['output'], w_log)

    combine.process()

    return None, None
