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


def in2out_pattern(number):
    """Transform input to output number pattern or image ID

    Parameters
    ----------
    number : string
        input number

    Returns
    -------
    number_final : string
        output number
    """

    # replace dots ('.') with dashes ('-') to avoid confusion
    # with file extension delimiters
    number_final = re.sub(r'\.', '-', number)

    # remove letters in number
    number_final = re.sub('[a-zA-Z]', '', number_final)

    return number_final


class GetImages(object):

    def __init__(self, copy, options, image_number_list,
                 input_numbering, input_file_pattern,
                 input_file_ext, w_log, check_existing_dir=None,
                 n_expected=None):
        """GetImages initiatliation.

        Parameters
        ----------
        copy : string
            copy/download method
        option : string
            copy options
        image_number_list : list of string self._copy = copy
            image numbers self._options = options
        input_numbering : string self._image_number_list = image_number_list
            numbering scheme, python regexp self._input_numbering = input_numbering
        input_file_pattern : list of strings self._input_file_ext = input_file_ext
            file pattern including input number template of input files self._w_log = w_log
        input_file_ext : list of strings
            input file extensions
        w_log :
            log file
        check_existing_dir : string, optional, default=None
            if not None, only retrieve image if not existing at this path (recursively)
        n_expected : int, optional, default=None
            number of expected files per type and ID to download/check for existence
        """

        self._copy = copy
        self._options = options
        self._image_number_list = image_number_list
        self._input_numbering = input_numbering
        self._input_file_pattern = input_file_pattern
        self._input_file_ext = input_file_ext
        self._w_log = w_log
        self._check_existing_dir = check_existing_dir
        self._n_expected = n_expected

    def get_file_list(self, dest_dir, output_file_pattern=None):
        """Return lists of file paths to copy.

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
                    # Transform input to output number patterns

                    number_final = in2out_pattern(number)

                    # Keep initial dot in extension
                    x = in_ext[1:]
                    x2 = re.sub(r'\.', '', x)
                    ext_final = in_ext[0] + x2
                    fbase = '{}{}'.format(output_file_pattern[i], number_final)
                else:
                    fbase = re.sub(self._input_numbering, number, in_pattern)
                    ext_final = in_ext

                if output_file_pattern and output_file_pattern[i] == '*':
                    # copy all input files to output dir, do not append
                    # extension
                    # fpath = '{}/.'.format(in_path)
                    fpath = in_path
                else:
                    fpath = '{}/{}{}'.format(in_path, fbase, ext_final)

                list_files_per_type.append(fpath)
            list_all_files.append(list_files_per_type)

        return list_all_files

    def copy(self, all_inputs, all_outputs):
        """Copy all files.

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
                    if path and len(path) == self._n_expected:
                        self._w_log.info('{} found, skipping'
                                         ''.format(path[0]))
                        continue
                self._w_log.info('Retrieving {}'.format(in_per_type[i]))
                self.copy_one(in_per_type[i], out_per_type[i])

    def copy_one(self, in_path, out_path):

        if self._copy == 'vos':
            sys.argv = []
            sys.argv.append('vcp')
            if self._options:
                for opt in self._options.split(' '):
                    sys.argv.append(opt)
            sys.argv.append(in_path)
            sys.argv.append(out_path)

            self._w_log.info('Command \'{}\''
                             ''.format(' '.join(sys.argv)))

            vcp = vosHandler('vcp')
            vcp()

        elif self._copy == 'symlink':
            src = in_path

            # Get all input file names if INPUT_FILE_PATTERN contains '*'
            all_src = glob.glob(src)
            dst = out_path
            for src in all_src:
                if os.path.isdir(dst):
                    # OUTPUT_FILE_PATTERN is '*', so dst is not regular file
                    # but directory. Append input file name
                    dst_name = '{}/{}'.format(dst, os.path.basename(src))
                else:
                    # dst is regular file
                    dst_name = dst
                os.symlink(src, dst_name)


def read_image_numbers(path):
    """Read image numbers from file.

    Parameters
    ----------
    path: string
        input file path

    Returns
    -------
    image_number_list: list of string
        image numbers
    """

    image_number_list = []
    with open(path) as f:
        for line in f:
            image_number_list.append(line.strip())

    return image_number_list


@module_runner(version='1.0',
               depends=['numpy'],
               run_method='serial',
               numbering_scheme='_0')
def get_images_runner(input_file_list, run_dirs, file_number_string,
                      config, w_log):

    # Input image numbers from all input tile files
    all_image_numbers = []
    for input_file in input_file_list:
        numbers_from_tile = read_image_numbers(input_file[0])
        all_image_numbers.append(numbers_from_tile)
    flat_list = [item for sublist in all_image_numbers for item in sublist]
    w_log.info('Number of images IDs = {}'.format(len(flat_list)))

    # Get unique number list
    image_number_list = list(set(flat_list))
    w_log.info('Number of unique image IDs = {}'.format(len(image_number_list)))

    # Read config file section

    # Copying/download method
    copy = config.get('GET_IMAGES_RUNNER', 'COPY')
    copy_ok = ['vos', 'symlink']
    if copy not in copy_ok:
        raise ValueError('key COPY={} is invalid, must be in {}'.format(copy, copy_ok))

    # Paths
    input_dir = config.getlist('GET_IMAGES_RUNNER', 'INPUT_PATH')
    nitem = len(input_dir)
    input_file_pattern = config.getlist('GET_IMAGES_RUNNER', 'INPUT_FILE_PATTERN')
    input_file_ext = config.getlist('GET_IMAGES_RUNNER', 'INPUT_FILE_EXT')
    output_file_pattern = config.getlist('GET_IMAGES_RUNNER', 'OUTPUT_FILE_PATTERN')

    input_numbering = config.get('GET_IMAGES_RUNNER', 'INPUT_NUMBERING')

    if config.has_option('GET_IMAGES_RUNNER', 'OUTPUT_PATH'):
        output_dir = config.getexpanded('GET_IMAGES_RUNNER', 'OUTPUT_PATH')
    else:
        output_dir = run_dirs['output']

    # Create array to make it compatible with input dir
    output_dir = [output_dir] * nitem

    if any(len(lst) != nitem for lst in [input_dir, input_file_pattern,
                                         input_file_ext, output_file_pattern]):
        raise ValueError('Lists INPUT_PATH ({}), INPUT_FILE_PATTERN ({}), '

                         'INPUT_FILE_EXT ({}), OUTPUT_FILE_PATTERN ({}) '
                         'need to have equal length'
                         ''.format(len(input_dir),
                                   len(input_file_pattern),
                                   len(input_file_ext),
                                   len(output_file_pattern)))

    # Copying/download method
    copy = config.get('GET_IMAGES_RUNNER', 'COPY')
    copy_ok = ['vos', 'symlink']
    if copy not in copy_ok:
        raise ValueError('key COPY={} is invalid, must be in {}'.format(copy, copy_ok))

    if config.has_option('GET_IMAGES_RUNNER', 'COPY_OPTIONS'):
        options = config.getexpanded('GET_IMAGES_RUNNER', 'COPY_OPTIONS')
    else:
        options = None

   # Check for already retrieved files
    if config.has_option('GET_IMAGES_RUNNER', 'CHECK_EXISTING_DIR'):
        check_existing_dir = config.getexpanded('GET_IMAGES_RUNNER', 'CHECK_EXISTING_DIR')
        if config.has_option('GET_IMAGES_RUNNER', 'N_EXPECTED'):
            n_expected = config.getint('GET_IMAGES_RUNNER', 'N_EXPECTED')
        else:
            n_expected = 1
    else:
        check_existing_dir = None
        n_expected = None

    inst = GetImages(copy, options, image_number_list, input_numbering,
                     input_file_pattern, input_file_ext, w_log,
                     check_existing_dir=check_existing_dir, n_expected=n_expected)

    # Assemble input and output file lists
    all_inputs = inst.get_file_list(input_dir)
    all_outputs = inst.get_file_list(output_dir, output_file_pattern=output_file_pattern)
    inst.copy(all_inputs, all_outputs)

    return None, None
