"""GET IMAGES.

This module copies all images required for processing.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import glob
import os
import re
import sys

from shapepipe.modules.module_decorator import module_runner
from shapepipe.utilities.canfar import vosHandler


# pragma: no cover
def read_image_numbers(path):
    """Read Image Numbers.

    Read image numbers from file.

    Parameters
    ----------
    path : str
        Input file path

    Returns
    -------
    list
        Image numbers

    """
    image_number_list = []
    with open(path) as file:
        for line in file:
            image_number_list.append(line.strip())

    return image_number_list


def in2out_pattern(number):
    """Get In2out Pattern.

    Transform input to output number pattern or image ID.

    Parameters
    ----------
    number : str
        Input number

    Returns
    -------
    str
        Output number

    """
    # replace dots ('.') with dashes ('-') to avoid confusion
    # with file extension delimiters
    number_final = re.sub(r'\.', '-', number)

    # remove letters in number
    number_final = re.sub('[a-zA-Z]', '', number_final)
    # make robust for more generalized file names
    number_final = re.sub(r'_', '', number_final)
    return number_final


class GetImages(object):
    """Get Images.

    Class handling retrieval of input images.

    Parameters
    ----------
    retrieve_method : str
        Copy/download method
    retrieve_option : str
        Retrieve options
    input_file_list : list
        Input files
    input_numbering : str
        Numbering scheme, python regexp
    input_file_pattern : list
        File pattern including input number template of input files
    input_file_ext : list
        Input file extensions
    output_file_pattern : list
        Output file patterns
    w_log : logging.Logger
        Log file
    check_existing_dir : str, optional
        If not ``None``, only retrieve image if not existing at this
        path (recursively)
    n_expected : int, optional
        Number of expected files per type and ID to download/check for
        existence
    n_try : int, optional
        Number of attempts for VOs download, default is ``3``

    """

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
        """Process.

        Main function to process GetImages.

        Parameters
        ----------
        input_dir : str
            Input directory
        output_dir : str
            Output directory

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
                + f'INPUT_FILE_EXT ({len(self._input_file_ext)}), '
                + f'OUTPUT_FILE_PATTERN ({len(self._output_file_pattern)}) '
                + 'need to have equal length'
            )

        # Assemble input and output file lists
        all_inputs = self.get_file_list(
            image_number_list,
            input_dir,
            use_output_file_pattern=False
        )
        all_outputs_orig = self.get_file_list(
            image_number_list,
            output_dir,
            use_output_file_pattern=False
        )
        all_outputs_renamed = self.get_file_list(
            image_number_list,
            output_dir,
            use_output_file_pattern=True
        )

        # Retrieve files
        self.retrieve(all_inputs, all_outputs_orig, all_outputs_renamed)

    def get_file_list(
        self,
        image_number_list,
        dest_dir,
        use_output_file_pattern=False,
    ):
        """Get File List.

        Return lists of file paths to retrieve.

        Parameters
        ----------
        image_number_list : list
            Image numbers
        dest_dir : list
            Input directory or url
        use_output_file_pattern : bool, optional
            If ``True``, use output file base patterns excluding numbering
            scheme; if ``False``, use input file patterns; default is ``False``

        Returns
        -------
        list
            Complete file paths, one list for each input file type

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

                if (
                    use_output_file_pattern
                    and self._output_file_pattern[idx] == '*'
                ):
                    # retrieve all input files to output dir, do not append
                    # extension
                    fpath = in_path
                else:
                    fpath = f'{in_path}/{fbase}{ext_final}'

                list_files_per_type.append(fpath)
            list_all_files.append(list_files_per_type)

        return list_all_files

    def retrieve(self, all_inputs, all_outputs_orig, all_outputs_renamed):
        """Retrieve.

        Retrieve all files.

        Parameters
        ----------
        all_inputs: list
            Input file paths, one list for each input file type
        all_outputs: list
            Output file paths, one list for each input file type

        """
        for in_per_type, out_per_type_orig, out_per_type_renamed in zip(
            all_inputs,
            all_outputs_orig,
            all_outputs_renamed
        ):
            for idx in range(len(in_per_type)):
                if self._check_existing_dir:
                    out_base = os.path.basename(out_per_type_orig[idx])
                    path = glob.glob(
                        f'{self._check_existing_dir}/**/{out_base}',
                        recursive=True,
                    )
                    if path:
                        if len(path) == self._n_expected:
                            self._w_log.info(
                                f'{path[0]} found, skipping download'
                            )
                            continue
                        else:
                            self._w_log.info(
                                f'{len(path)} instead of {self._n_expected} '
                                + 'existing files found at'
                                + f' {self._check_existing_dir}'
                                + ', downloading images'
                            )
                    else:
                        self._w_log.info(
                            'No existing images found at'
                            + f' {self._check_existing_dir},'
                            + ' downloading images'
                        )
                self.retrieve_one(
                    in_per_type[idx],
                    out_per_type_orig[idx],
                    out_per_type_renamed[idx]
                )

    def retrieve_one(self, in_path, out_path_orig, out_path_renamed):
        """Retrieve One.

        Retrieve one file.

        Parameters
        ----------
        in_path : str
            Input path
        out_path : str
            Output path

        """
        if self._retrieve_method == 'vos':
            sys.argv = []
            sys.argv.append('vcp')
            if self._retrieve_options:
                for opt in self._retrieve_options.split(' '):
                    sys.argv.append(opt)
            sys.argv.append(in_path)
            sys.argv.append(out_path_orig)

            log_cmd = ' '.join(sys.argv)
            vcp = vosHandler('vcp')
            self._w_log.info(log_cmd)

            # Download file from VOSpace
            attempt = 0
            while attempt < self._n_try:
                try:
                    vcp()
                    self._w_log.info(
                        'Success of command vcp after '
                        + f'{attempt}/{self._n_try} attempts'
                    )
                    break
                except Exception:
                    attempt += 1
                    self._w_log.info(
                        'Error with command vcp, attempt '
                        + f'{attempt}/{self._n_try}'
                    )

            sys.argv = None

            # Create symbolic link to downloaded file with
            # link name in ShapePipe numbering format
            os.symlink(out_path_orig, out_path_renamed)

        elif self._retrieve_method == 'symlink':
            src = in_path

            # Get all input file names if INPUT_FILE_PATTERN contains '*'
            all_src = glob.glob(src)
            if len(all_src) == 0:
                raise IndexError(
                    f'No input file found corresponding to \'{src}\''
                )

            dst = out_path_renamed
            for src in all_src:
                if os.path.isdir(dst):
                    # OUTPUT_FILE_PATTERN is '*', so dst is not regular file
                    # but directory. Append input file name
                    dst_name = f'{dst}/{os.path.basename(src)}'
                else:
                    # dst is regular file
                    dst_name = dst
                os.symlink(src, dst_name)
