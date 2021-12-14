# -*- coding: utf-8 -*-

"""GET IMAGES

This module copies all images required for processing

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2019 - 2021

:Package: ShapePipe

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.utilities.canfar import vosHandler

import os
import re
import sys
import glob


# pragma: no cover
def read_image_numbers(path):
    """Read Image Numbers
    Read image numbers from file.

    Parameters
    ----------
    path : str
        input file path

    Returns
    -------
    image_number_list : list of str
        image numbers
    """

    image_number_list = []
    with open(path) as f:
        for line in f:
            image_number_list.append(line.strip())

    return image_number_list


def in2out_pattern(number):
    """In2out Pattern
    Transform input to output number pattern or image ID

    Parameters
    ----------
    number : str
        input number

    Returns
    -------
    number_final : str
        output number
    """

    # replace dots ('.') with dashes ('-') to avoid confusion
    # with file extension delimiters
    number_final = re.sub(r'\.', '-', number)

    # remove letters in number
    number_final = re.sub('[a-zA-Z]', '', number_final)

    return number_final


class GetImages(object):

    def __init__(
        self,
        retrieve_method,
        retrieve_options,
        input_file_list,
        input_numbering,
        input_file_pattern,
        input_file_ext,
        output_file_pattern,
        w_log,
        check_existing_dir=None,
        n_expected=None,
        n_try=3,
    ):
        """Get Images

        Class handling retrieval of input images

        Parameters
        ----------
        retrieve_method : str
            copy/download method
        retrieve_option : str
            retrieve options
        input_file_list : list of str
            input files
        input_numbering : str
            numbering scheme, python regexp
        input_file_pattern : list of str
            file pattern including input number template of input files
        input_file_ext : list of str
            input file extensions
        output_file_pattern : list of str
            output file patterns
        w_log : logging.Logger
            log file
        check_existing_dir : str, optional, default=None
            if not None, only retrieve image if not existing at this
            path (recursively)
        n_expected : int, optional, default=None
            number of expected files per type and ID to download/check for
            existence
        n_try : int, optional, default=3
            number of attempts for VOs download
        """

        self._retrieve_method = retrieve_method
        self._retrieve_options = retrieve_options
        self._input_file_list = input_file_list
        self._input_numbering = input_numbering
        self._input_file_pattern = input_file_pattern
        self._input_file_ext = input_file_ext
        self._output_file_pattern = output_file_pattern
        self._w_log = w_log
        self._check_existing_dir = check_existing_dir
        self._n_expected = n_expected
        self._n_try = n_try

    def process(self, input_dir, output_dir):
        """Process

        Main function to process GetImages

        Parameters
        ----------
        input_dir : str
            input directory
        output_dir : str
            output directory
        """

        # Input image numbers from all input tile files
        all_image_numbers = []
        for input_file in self._input_file_list:
            numbers_from_tile = read_image_numbers(input_file[0])
            all_image_numbers.append(numbers_from_tile)

        # List of unique input images
        flat_list = [item for sublist in all_image_numbers for item in sublist]
        self._w_log.info(f'Number of total image IDs = {len(flat_list)}')

        # Get unique number list
        image_number_list = list(set(flat_list))
        self._w_log.info(
            f'Number of unique image IDs = {len(image_number_list)}'
        )

        # Create array to make it compatible with input dir
        nitem = len(input_dir)

        # Make sure output_dir is list and compatible to input lists
        output_dir = [output_dir] * nitem

        # Check consistency of list lengths
        if any(
            len(lst) != nitem for lst in [
                input_dir,
                self._input_file_pattern,
                self._input_file_ext,
                self._output_file_pattern
            ]
        ):
            raise ValueError(
                f'Lists INPUT_PATH ({len(input_dir)}), '
                + f'INPUT_FILE_PATTERN ({len(self._input_file_pattern)}), '
                + f'INPUT_FILE_EXT ({llen(self._input_file_ext)}), '
                + f'OUTPUT_FILE_PATTERN ({len(self._output_file_pattern)}) '
                + 'need to have equal length'
            )

        # Assemble input and output file lists
        all_inputs = self.get_file_list(
            image_number_list,
            input_dir,
            use_output_file_pattern=False
        )
        all_outputs = self.get_file_list(
            image_number_list,
            output_dir,
            use_output_file_pattern=True
        )

        # Retrieve files
        self.retrieve(all_inputs, all_outputs)

    def get_file_list(self, image_number_list, dest_dir, use_output_file_pattern=False):
        """Get File List
        Return lists of file paths to retrieve.

        Parameters
        ----------
        image_number_list : list of str
            image numbers
        dest_dir : list of str
            input directory or url
        use_output_file_pattern : bool, optional, default=False
            if True, use output file base patterns excluding numbering scheme;
            if False, use input file patterns

        Returns
        -------
        list_all_files : list of list of str
            complete file paths, one list for each input file type
        """

        list_all_files = []
        for idx in range(len(dest_dir)):
            in_path = dest_dir[idx]
            in_pattern = self._input_file_pattern[idx]
            in_ext = self._input_file_ext[idx]

            list_files_per_type = []
            for number in image_number_list:

                if use_output_file_pattern:
                    # Transform input to output number patterns

                    number_final = in2out_pattern(number)

                    # Keep initial dot in extension
                    x = in_ext[1:]
                    x2 = re.sub(r'\.', '', x)
                    ext_final = in_ext[0] + x2
                    fbase = (
                        f'{self._output_file_pattern[idx]}{number_final}'
                    )
                else:
                    fbase = re.sub(self._input_numbering, number, in_pattern)
                    ext_final = in_ext

                if use_output_file_pattern and self._output_file_pattern[idx] == '*':
                    # retrieve all input files to output dir, do not append
                    # extension
                    # fpath = '{}/.'.format(in_path)
                    fpath = in_path
                else:
                    fpath = '{}/{}{}'.format(in_path, fbase, ext_final)

                list_files_per_type.append(fpath)
            list_all_files.append(list_files_per_type)

        return list_all_files

    def retrieve(self, all_inputs, all_outputs):
        """Retrieve
        Retrieve all files.

        Parameters
        ----------
        all_inputs: list of list of str
            input file paths, one list for each input file type
        all_outputs: list of list of str
            output file paths, one list for each input file type
        """

        for in_per_type, out_per_type in zip(all_inputs, all_outputs):
            for idx in range(len(in_per_type)):
                if self._check_existing_dir:
                    out_base = os.path.basename(in_per_type[idx])
                    path = glob.glob('{}/**/{}'
                                     ''.format(self._check_existing_dir,
                                               out_base),
                                     recursive=True)
                    if path and len(path) == self._n_expected:
                        self._w_log.info('{} found, skipping'
                                         ''.format(path[0]))
                        continue
                self.retrieve_one(in_per_type[idx], out_per_type[idx])

    def retrieve_one(self, in_path, out_path):
        """Retrieve One
        Retrieve one file.

        Parameters
        ----------
        in_path : str
            input path
        out_path : str
            output path
        """

        if self._retrieve_method == 'vos':
            sys.argv = []
            sys.argv.append('vcp')
            if self._retrieve_options:
                for opt in self._retrieve_options.split(' '):
                    sys.argv.append(opt)
            sys.argv.append(in_path)
            sys.argv.append(out_path)

            log_cmd = ' '.join(sys.argv)
            vcp = vosHandler('vcp')
            print(log_cmd)

            attempt = 0
            while attempt < self._n_try:
                try:
                    vcp()
                    self._w_log.info(f'Success of command vcp after {attempt}/{self._n_try} attempts')
                    break
                except:
                    attempt += 1
                    self._w_log.info(f'Error with command vcp, attempt {attempt}/{self._n_try}')

            sys.argv = None

        elif self._retrieve_method == 'symlink':
            src = in_path

            # Get all input file names if INPUT_FILE_PATTERN contains '*'
            all_src = glob.glob(src)
            if len(all_src) == 0:
                raise IndexError(
                    f'No input file found corresponding to \'{src}\''
                )

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
