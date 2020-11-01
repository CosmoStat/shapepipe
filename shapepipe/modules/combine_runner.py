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

    """

    def __init__(self, input_file_list, w_log):

        self._input_file_list = input_file_list
        self._w_log = w_log

    def process(self):

        pass

        for input_file in self._input_file_list:
            print(input_file)

        # Create dir
        # Create links
        # update log


@module_runner(version='1.0',
               run_method='serial')
def combine_runner(input_file_list, run_dirs, file_number_string,
                   config, w_log):


    combine = Combine(input_file_list, w_log)

    combine.process()

    return None, None
