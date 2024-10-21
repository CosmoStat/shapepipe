"""MERGE HEADERS.

This module merges the output *header* files of the ``split_exp``
module. It creates a binnary file that contains the WCS of each CCD for each
exposure.

:Author: Axel Guinot

"""

import os
import re

import numpy as np
from sqlitedict import SqliteDict


def merge_headers(input_file_list, output_dir):
    """Merge Headers.

    This function opens the files in the input file list and merges them into
    a SqliteDict file provided they match the appropriate pattern.

    Parameters
    ----------
    input_file_list : list
        List of input files
    output_dir : str
        Output path

    Raises
    ------
    TypeError
        For invalid ``output_dir`` type

    """
    if not isinstance(output_dir, str):
        raise TypeError(
            "Output directory for merge headers must be a string "
            + f"not {type(output_dir)}."
        )

    # Open SqliteDict file
    final_file = SqliteDict(f"{output_dir}/log_exp_headers.sqlite")
    # Set matching pattern
    pattern = "headers-"

    for file_path in input_file_list:
        # Extract file path
        file_path_scalar = file_path[0]
        # Check file name against pattern
        file_name = os.path.split(file_path_scalar)[1]
        file_base_name = os.path.splitext(file_name)[0]
        pattern_split = re.split(pattern, file_base_name)
        if len(pattern_split) < 2:
            raise IndexError(
                f'Regex "{pattern}" not found in base name "{file_base_name}".'
            )
        key = pattern_split[1]
        # Load Numpy binary file
        final_file[key] = np.load(file_path_scalar, allow_pickle=True)

    # Close file
    final_file.commit()
    final_file.close()
