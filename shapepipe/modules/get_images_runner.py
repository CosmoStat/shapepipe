# -*- coding: utf-8 -*-

"""GET IMAGES RUNNER

This module copies all images required for processing

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
import re
import sys


class GetImages(object):

    def __init__(self, copy, options, image_number_list,
                 input_numbering, input_file_pattern,
                 input_file_ext, w_log):
        """GetImages initiatliation.

        Parameters
        ----------
        copy: string
            copy/download method
        option: string
            copy options
        image_number_list: list of string self._copy = copy
            image numbers self._options = options
        input_numbering: string self._image_number_list = image_number_list
            numbering scheme, python regexp self._input_numbering = input_numbering
        input_file_pattern: list of strings self._input_file_ext = input_file_ext
            file pattern including input number template of input files self._w_log = w_log
        input_file_ext: list of strings
            input file extensions
        w_log:
            log file
        """

        self._copy = copy
        self._options = options
        self._image_number_list = image_number_list
        self._input_numbering = input_numbering
        self._input_file_pattern = input_file_pattern
        self._input_file_ext = input_file_ext
        self._w_log = w_log

    def get_file_list(self, dest_dir, output_file_pattern=None):
        """Return list of file paths to copy.

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
                self.copy_one(in_per_type[i], out_per_type[i])

    def copy_one(self, in_path, out_path):

        if self._copy == 'vos':
            sys.argv = []
            sys.argv.append('vcp')
            for opt in self._options.split(' '):
                sys.argv.append(opt)
            sys.argv.append(in_path)
            sys.argv.append(out_path)

            try:
                from vos.commands.vcp import vcp
            except:
                raise ImportError('vos modules not found, re-install ShapePipe with \'install_pipeline --vos\'')

            try:
                vcp()
            except:
                raise ValueError('Error in \'vcp\' command: \'{}\''.format(' '.join(sys.argv)))


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
               depends=['numpy', 'vos'],
               run_method='serial')
def get_images_runner(input_file_list, run_dirs, file_number_string,
                      config, w_log):

    # Input image numbers from all input tile files
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
        output_dir = [output_dir] * nitem

    if any(len(lst) != nitem for lst in [input_file_pattern, input_file_ext,
                                         output_dir, output_file_pattern]):
        raise ValueError('Lists INPUT_DIR, INPUT_FILE_PATTERN, INPUT_FILE_EXT, '
                         'OUTPUT_DIR, OUTPUT_FILE_PATTERN  need to '
                         'have equal length')

    # Copying/download method
    copy = config.get('GET_IMAGES_RUNNER', 'COPY')
    copy_ok = ['vos', 'symlink']
    if copy not in copy_ok:
        raise ValueError('key COPY={} is invalid, must be in {}'.format(copy, copy_ok))

    if config.has_option('GET_IMAGES_RUNNER', 'COPY_OPTIONS'):
        options = config.get('GET_IMAGES_RUNNER', 'COPY_OPTIONS')

    inst = GetImages(copy, options, image_number_list, input_numbering,
                     input_file_pattern, input_file_ext, w_log)

    # Assemble input and output file lists
    all_inputs = inst.get_file_list(input_dir)
    all_outputs = inst.get_file_list(output_dir, output_file_pattern=output_file_pattern)
    inst.copy(all_inputs, all_outputs)

    return None, None
