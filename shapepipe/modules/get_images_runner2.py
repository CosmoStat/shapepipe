# -*- coding: utf-8 -*-

"""GET IMAGES RUNNER

This module copies all images required for processing

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.utilities.canfar import vosHandler

import os
import re
import sys
import glob


class GetImages(object):

    def __init__(self, retrieve, options, image_number_list,
                 input_numbering, input_file_pattern,
                 input_file_ext, w_log, check_existing_dir=None):
        """GetImages initialisation.

        Parameters
        ----------
        retrieve : string
            retrieve/download method
        option : string
            retrieve options
        image_number_list : list of string
            image numbers self._options = options
        input_numbering : string
            numbering scheme, python regexp
        input_file_pattern : list of strings
            file pattern including input number template of input files
        input_file_ext : list of strings
            input file extensions
        w_log:
            log file
        check_existing_dir : string, optional, default=None
            if not None, only retrieve image if not existing at this directory

        Returns
        --------
        self: GetImages class instance
        """

        self._retrieve = retrieve
        self._options = options
        self._image_number_list = image_number_list
        self._input_numbering = input_numbering
        self._input_file_pattern = input_file_pattern
        self._input_file_ext = input_file_ext
        self._w_log = w_log
        self._check_existing_dir = check_existing_dir

    def get_file_list(self, dest_dir, output_file_pattern=None):
        """Return list of file paths to retrieve.

        Parameters
        ----------
        dest_dir: list of stringsr
            input directory or url
        output_file_pattern: list of strings, optional, default=None
            output file base patterns excluding numbering scheme,
            if None use input file patterns

        Returns
        -------
        list_all_files: list of list of strings
            complete file paths, one list for each input file type
        """

        list_all_files = []
        for i in range(len(dest_dir)):
            in_path = dest_dir[i]
            in_pattern = self._input_file_pattern[i]
            in_ext = self._input_file_ext[i]

            list_files_per_type = []
            for number in self._image_number_list:

                if output_file_pattern:
                    # For output:
                    # - replace dots ('.') with dashes ('-') to avoid confusion
                    #   with file extension delimiters
                    # - remove letters in number

                    number_final = re.sub(r'\.', '-', number)
                    number_final = re.sub('[a-zA-Z]', '', number_final)

                    # Keep initial dot in extension
                    x = in_ext[1:]
                    x2 = re.sub(r'\.', '', x)
                    ext_final = in_ext[0] + x2
                    fbase = '{}{}'.format(output_file_pattern[i], number_final)
                else:
                    fbase = re.sub(self._input_numbering, number, in_pattern)
                    ext_final = in_ext

                    # Remove 'p' for LSB images
                    # fbase = re.sub('p', '', fbase)

                fpath = '{}/{}{}'.format(in_path, fbase, ext_final)
                list_files_per_type.append(fpath)
            list_all_files.append(list_files_per_type)

        return list_all_files

    def retrieve(self, all_inputs, all_outputs):
        """Retreve all files.

        Parameters
        ----------
        all_inputs: list of list of strings
            input file paths, one list for each input file type
        all_outputs: list of list of strings
            output file paths, one list for each input file type
        """

        for in_per_type, out_per_type in zip(all_inputs, all_outputs):
            for i in range(len(in_per_type)):
                if self._check_existing_dir:
                    out_base = os.path.basename(in_per_type[i])
                    path = glob.glob('{}/**/{}'
                                     ''.format(self._check_existing_dir,
                                               out_base),
                                     recursive=True)
                    if path:
                        self._w_log.info('{} found, skipping'
                                         ''.format(path[0]))
                        continue
                self._w_log.info('Retrieving {}'.format(in_per_type[i]))
                self.retrieve_one(in_per_type[i], out_per_type[i])

    def retrieve_one(self, in_path, out_path):
        """Retrieve one file.

        Parameters
        ----------
        in_path : string
            input file path
        out_path : string
            output file path
        """

        if self._retrieve == 'vos':
            sys.argv = []
            sys.argv.append('vcp')
            if self._options:
                for opt in self._options.split(' '):
                    sys.argv.append(opt)
            sys.argv.append(in_path)
            sys.argv.append(out_path)

            vcp = vosHandler('vcp')
            vcp()

        elif self._retrieve == 'symlink':
            src = in_path
            dst = out_path
            os.symlink(src, dst)
            if not os.path.exists(src):
                w_log.info('Warning: Source of symlink \'{}\' '
                           'does not exist'
                           ''.format(src))


def read_image_numbers(path):
    """Read image numbers from file.

    Parameters
    ----------
    path : string
        input file path

    Returns
    -------
    image_number_list : list of string
        image numbers
    """

    image_number_list = []
    with open(path) as f:
        for line in f:
            image_number_list.append(line.strip())

    return image_number_list


@module_runner(version='1.0',
               depends=['numpy'],
               run_method='serial')
def get_images_runner2(input_file_list, run_dirs, file_number_string,
                       config, w_log):

    # Input image numbers from all input tile fils
    all_image_numbers = []
    for input_file in input_file_list:
        numbers_from_tile = read_image_numbers(input_file[0])
        all_image_numbers.append(numbers_from_tile)
    flat_list = [item for sublist in all_image_numbers for item in sublist]
    w_log.info('{} image numbers read in total'.format(len(flat_list)))

    # Get unique number list
    image_number_list = list(set(flat_list))
    w_log.info('{} unique exposures numbers'.format(len(image_number_list)))

    # Read config file section

    # Paths
    input_dir = config.getlist('GET_IMAGES_RUNNER2', 'INPUT_PATH')
    nitem = len(input_dir)
    input_file_pattern = config.getlist('GET_IMAGES_RUNNER2', 'INPUT_FILE_PATTERN')
    input_file_ext = config.getlist('GET_IMAGES_RUNNER2', 'INPUT_FILE_EXT')
    output_file_pattern = config.getlist('GET_IMAGES_RUNNER2', 'OUTPUT_FILE_PATTERN')

    input_numbering = config.get('GET_IMAGES_RUNNER2', 'INPUT_NUMBERING')

    if config.has_option('GET_IMAGES_RUNNER2', 'OUTPUT_PATH'):
        output_dir = config.getexpanded('GET_IMAGES_RUNNER2', 'OUTPUT_PATH')
    else:
        output_dir = run_dirs['output']
        output_dir = [output_dir] * nitem

    if any(len(lst) != nitem for lst in [input_file_pattern, input_file_ext,
                                         output_dir, output_file_pattern]):
        raise ValueError('Lists INPUT_DIR, INPUT_FILE_PATTERN, INPUT_FILE_EXT, '
                         'OUTPUT_DIR, OUTPUT_FILE_PATTERN  need to '
                         'have equal length')

    # Method to retrieve images
    retrieve = config.get('GET_IMAGES_RUNNER2', 'RETRIEVE')
    retrieve_ok = ['vos', 'symlink']
    if retrieve not in retrieve_ok:
        raise ValueError('key RETRIEVE={} is invalid, must be in {}'.format(retrieve, retrieve_ok))

    if config.has_option('GET_IMAGES_RUNNER2', 'RETRIEVE_OPTIONS'):
        options = config.getexpanded('GET_IMAGES_RUNNER2', 'RETRIEVE_OPTIONS')
    else:
        options = None

    # Check for already retrieved files
    if config.has_option('GET_IMAGES_RUNNER2', 'CHECK_EXISTING_DIR'):
        check_existing_dir = config.getexpanded('GET_IMAGES_RUNNER2', 'CHECK_EXISTING_DIR')
    else:
        check_existing_dir = None

    inst = GetImages(retrieve, options, image_number_list, input_numbering,
                     input_file_pattern, input_file_ext, w_log,
                     check_existing_dir=check_existing_dir)

    # Assemble input and output file lists
    all_inputs = inst.get_file_list(input_dir)
    all_outputs = inst.get_file_list(output_dir, output_file_pattern=output_file_pattern)
    inst.retrieve(all_inputs, all_outputs)

    return None, None
