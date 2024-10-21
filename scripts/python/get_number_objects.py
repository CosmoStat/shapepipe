#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script get_number_objects.py

Get number of objects in a (last-run SExtractor) catalogue.

:Author: Martin Kilbinger <martin.kilblinger@cea.fr>

"""

import sys
import copy
import glob

from optparse import OptionParser
from astropy.io import fits

from shapepipe.pipeline.run_log import get_last_dir, get_all_dirs
from shapepipe.utilities import cfis


class param:
    """General class to store (default) variables"""

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)


def params_default():
    """Set default parameter values.

    Returns
    -------
    class param
        parameter values

    """
    p_def = param(
        input_path=".",
        input_name_base="final_cat",
        hdu_num=1,
    )

    return p_def


def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class param
        parameter values

    Returns
    -------
    list
        command line options
        command line str

    """
    usage = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # IO
    parser.add_option(
        "-i",
        "--input_path",
        dest="input_path",
        type="string",
        default=p_def.input_path,
        help=f"input path, default='{p_def.input_path}'",
    )
    parser.add_option(
        "-n",
        "--input_name_base",
        dest="input_name_base",
        type="string",
        default=p_def.input_name_base,
        help=f"input name base, default='{p_def.input_name_base}'",
    )
    parser.add_option(
        "-l",
        "--list_tile_ID_path",
        dest="tile_ID_list_path",
        type="string",
        default=None,
        help=f"tile ID list, default: Use all data in input files",
    )

    # Control
    parser.add_option(
        "-p",
        "--param_path",
        dest="param_path",
        type="string",
        default=None,
        help="parameter file path, default=None",
    )

    parser.add_option(
        "",
        "--hdu_num",
        dest="hdu_num",
        type="int",
        default=p_def.hdu_num,
        help=f"input HDU number, default='{p_def.hdu_num}'",
    )

    parser.add_option(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="verbose output",
    )

    options, args = parser.parse_args()

    return options, args


def check_options(options):
    """Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Returns
    -------
    bool
        Result of option check. False if invalid option value.

    """
    return True


def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.

    Parameters
    ----------
    p_def:  class param
        parameter values
    optiosn: tuple
        command line options

    Returns
    -------
    class param
        updated paramter values

    """
    param = copy.copy(p_def)

    # Update keys in param according to options values
    for key in vars(param):
        if key in vars(options):
            setattr(param, key, getattr(options, key))

    # Add remaining keys from options to param
    for key in vars(options):
        if not key in vars(param):
            setattr(param, key, getattr(options, key))

    # Do extra stuff if necessary

    return param


def main(argv=None):

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    # Save command line arguments to log file
    f_log = cfis.log_command(argv, close_no_return=False)

    pattern = "sexcat"
    run_log_file = "output/log_run_sp.txt"

    # For v1
    # module = 'sextractor_runner_run_1'

    # For v2
    module = "sextractor_runner"

    # Get all sextractor output directories
    all_dir = get_all_dirs(run_log_file, module)
    paths = []

    # Find tile runs
    for path in all_dir:
        if "run_sp_tile_Sx" in path:
            paths.append(path)
    paths = sorted(paths)

    # Get latest run
    last_dir = paths[-1]

<<<<<<< HEAD
    # Get all output SExtractor catalogues
    file_list = glob.glob(f'{last_dir}/{pattern}*.fits')
=======
    file_list = glob.glob(f"{last_dir}/{pattern}*.fits")
>>>>>>> origin/v1.4
    if len(file_list) == 0:
        raise ValueError(f"No files {last_dir}/{pattern}*.fits found")

    # Add up number of objects over all catalogues
    n_obj = 0
    hdu_no = -1
    for fpath in file_list:
        hdu_list = fits.open(fpath)
        header = hdu_list[-1].header
        n_obj += int(header["NAXIS2"])

    # Compute average
    n_obj = int(n_obj / len(file_list))

    print(n_obj)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
