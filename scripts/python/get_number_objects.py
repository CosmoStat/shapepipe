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
        input_path='.',
        input_name_base='final_cat',
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
    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # IO
    parser.add_option(
        '-i',
        '--input_path',
        dest='input_path',
        type='string',
        default=p_def.input_path,
        help=f'input path, default=\'{p_def.input_path}\''
    )
    parser.add_option(
        '-n',
        '--input_name_base',
        dest='input_name_base',
        type='string',
        default=p_def.input_name_base,
        help=f'input name base, default=\'{p_def.input_name_base}\''
    )
    parser.add_option(
        '-l',
        '--list_tile_ID_path',
        dest='tile_ID_list_path',
        type='string',
        default=None,
        help=f'tile ID list, default: Use all data in input files'
    )

    # Control
    parser.add_option(
        '-p',
        '--param_path',
        dest='param_path',
        type='string',
        default=None,
        help='parameter file path, default=None'
    )

    parser.add_option(
        '',
        '--hdu_num',
        dest='hdu_num',
        type='int',
        default=p_def.hdu_num,
        help=f'input HDU number, default=\'{p_def.hdu_num}\''
    )

    parser.add_option(
        '-v',
        '--verbose',
        dest='verbose',
        action='store_true',
        help='verbose output'
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


def read_param_file(path, verbose=False):
    """Return parameter list read from file

    Parameters
    ----------
    path: str
        input file name
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    list of str
        parameter names

    """
    param_list = []

    if path:

        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue 
                entry = line.rstrip()
                if not entry or entry == '':
                    continue
                param_list.append(entry)

    if verbose:
        if len(param_list) > 0: 
            print(f'Copying {len(param_list)} columns', end='')
        else:
            print('Copying all columns', end='')
        print(' into merged catalogue')

    # Check for multiples
    multiples = []
    for param in param_list:
        if param_list.count(param) > 1:
            multiples.append(param)

    if len(multiples) > 0:
        print('The following parameters are more than one times '
              'in the parameter file: ', end='')
        for m in multiples:
            print(m, end=' ')
        print()
        raise ValueError('Multiple identical keys found')

    return param_list
                            

def get_data(path, hdu_num, param_list):
    """Return data of selected columns from FITS file.

    Parameters
    ----------
    path: str
        input file name
    hdu_num: int
        HDU number
    param_list: list of str
        parameters to be extracted. If none, copy
        all columns

    Returns
    -------
    numpy array
        data columns

    """
    hdu_list = fits.open(path)
    hdu = hdu_list[hdu_num]

    if param_list:
        cols = []
        for p in param_list:
            cols.append(hdu.columns[p]) 
        coldefs = fits.ColDefs(cols)
        hdu_new = fits.BinTableHDU.from_columns(coldefs)
        d = hdu_new.data
    else:
        d = hdu.data

    return d


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

    pattern = 'sexcat'
    run_log_file = 'output/log_run_sp.txt'

    # For v1
    #module = 'sextractor_runner_run_1'

    # For v2
    module = "sextractor_runner"
    all_dir = get_all_dirs(run_log_file, module)
    paths = []
    for path in all_dir:
        if "run_sp_tile_Sx" in path:
            paths.append(path)
    paths = sorted(paths)
    last_dir = paths[-1]

    file_list = glob.glob(f'{last_dir}/{pattern}*.fits')
    if len(file_list) == 0:
        raise ValueError(f'No files {last_dir}/{pattern}*.fits found')

    n_obj = 0
    hdu_no = -1
    for fpath in file_list:
        hdu_list = fits.open(fpath)
        header = hdu_list[-1].header
        n_obj += int(header['NAXIS2'])

    n_obj = int(n_obj / len(file_list))

    print(n_obj)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
