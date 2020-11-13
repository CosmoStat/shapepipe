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
from shapepipe.pipeline.file_handler import find_files
from shapepipe.pipeline import file_io as io

import os
import re

class Combine(object):
    """ Combine

    Combine results from previous runs.

    Parameters
    ----------
    input_dir : list of strings
        input directories
    file_pattern : list of strings
        input file patterns
    file_ex : list of strings
        input file extensions
    output_dir : string
        output directory name
    w_log :
        log file handler

    """

    def __init__(self, input_dir, file_pattern, file_ext, output_dir, w_log):

        self._input_dir = input_dir
        self._file_pattern = file_pattern
        self._file_ext = file_ext
        self._output_dir = output_dir
        self._w_log = w_log

    def process(self):

        n_skipped = 0
        n_created = 0
        n_input = 0

        for input_dir in self._input_dir:

            for file_pattern, file_ext in zip(self._file_pattern, self._file_ext):
                file_list = find_files(input_dir, file_pattern, file_ext)
                n_input = n_input + len(file_list)

                # Loop over input types
                for input_path in file_list:

                    input_base_name = os.path.basename(input_path)
                    output_path = '{}/{}'.format(self._output_dir, input_base_name)
                    try:
                        os.symlink(input_path, output_path)
                        n_created = n_created + 1
                    except FileExistsError:
                        self._w_log.info('File {} exists, skipping symlink'.format(output_path))
                        n_skipped = n_skipped + 1

        self._w_log.info('Number of links created = {}'.format(n_created))
        self._w_log.info('Number of (duplicate) files skipped = {}'.format(n_skipped))
        self._w_log.info('Total number of input files = {}'
                         ''.format(n_input))


def my_format(path, name, ext=''):
    """ See FileHandler.format()

    """

    return '{}/{}{}'.format(path, name, ext)


def get_all(module, runs):
    """ See RunLog.get_all()

    """

    all_runs = [run for run in runs if module in
                run.split()[1].split(',')]
    if len(all_runs) == 0:
        raise RuntimeError('No previous run of module \'{}\' '
                           'found'.format(module))

    all_runs = all_runs[::-1]

    return all_runs


@module_runner(version='1.0',
               run_method='serial')
def combine_runner(input_file_list, run_dirs, file_number_string,
                   config, w_log):

    file_pattern = config.getlist('COMBINE_RUNNER', 'FILE_PATTERN')
    file_ext = config.getlist('COMBINE_RUNNER', 'FILE_EXT')

    output_dir = config.getexpanded('FILE', 'OUTPUT_DIR')
    run_log_file = my_format(output_dir,
                             config.get('FILE', 'RUN_LOG_NAME'),
                             '.txt')

    # See RunLog._get_list()
    with open(run_log_file, 'r') as run_log:
            lines = run_log.readlines()
    runs = [line.rstrip() for line in lines]

    dir_list = config.getlist('COMBINE_RUNNER', 'INPUT_DIR')
    # See FileHander._check_input_dir_list()
    input_dir = []
    for dir in dir_list:
        if 'all' in dir.lower():
            module = dir.lower().split(':')[1]
            all_runs = get_all(module, runs)
            input_dir.extend([my_format(my_format(
                                  run.split(' ')[0],
                                  module), 'output')
                                  for run in all_runs])

    print(input_dir)

    print(input_file_list)

    combine = Combine(input_dir, file_pattern, file_ext, run_dirs['output'], w_log)

    #combine.process()

    return None, None
