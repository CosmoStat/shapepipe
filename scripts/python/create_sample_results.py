#!/usr/bin/env python

"""Script create_sample_results.py

Create directory with links to results for a given (sub-)sample.

:Author: Martin Kilbinger

:Date: 07/2020
"""

import re
import os
import sys
import glob
import copy
import io
from contextlib import redirect_stdout
from optparse import OptionParser
import cfis
from shapepipe.utilities.file_system import mkdir


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class cfis.param
        parameter values
    """

    p_def = cfis.param(
    )

    return p_def


def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class cfis.param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # I/O
    parser.add_option('', '--input_IDs', dest='input_IDs', type='string',
         help='input tile ID file specifying sample')
    parser.add_option('-i', '--input_dir', dest='input_dir', type='string',
         help='input directory name')
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string',
         help='output directory name')

    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')

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
    erg: bool
        Result of option check. False if invalid option value.
    """

    if options.input_IDs is None:
        print('No input ID file list given (option \'--input_IDs\')')
        return False
    if options.input_dir is None:
        print('No input directory name given (option \'--input_dir\')')
        return False
    if options.output_dir is None:
        print('No output directory name given (option \'--output_dir\')')
        return False

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
    param: class param
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


def read_ID_list(input_ID_path, verbose=False):
    """Return input ID list from file.

    Parameters
    ----------
    input_ID_path: string
        input ID file path
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    input_IDs: list of strings
        input ID list
    """

    if verbose:
        print('Reading input ID list...')

    input_IDs = []
    with open(input_ID_path) as f:
        for line in f:
            input_IDs.append(line.rstrip())

    if verbose:
        print('{} IDs found in input file'.format(len(input_IDs)))

    return input_IDs


def create_links(input_dir, output_dir, input_IDs, result_base_names, verbose=False):
    """Create symbolic links to result files corresponding to (sub-)sample.

    Parameters
    ----------
    input_dir: string
        input directory
    output_dir: string
        output directory
    input_IDs: list of strings
        tile ID list
    results_base_names: list of strings
        file base names
    verbose: bool, optional, default=False
        verbose output if True
    """

    if verbose:
        print('Creating links...')

    n_total = {}
    n_created = 0
    n_existed = 0
    for ID in input_IDs:
        n_total[ID] = 0
        for base in result_base_names:
            name = '{}_{}.tgz'.format(base, ID)
            src = '{}/{}'.format(os.path.abspath(input_dir), name)
            link_name = '{}/{}'.format(output_dir, name)

            #if verbose:
                #print('Creating link {} <- {}'.format(src, link_name))

            if not os.path.exists(src):
                #raise IOError('Source file \'{}\' does not exist'.format(src))
                print('Source file \'{}\' does not exist, skipping'.format(src))
            elif not os.path.exists(link_name):
                os.symlink(src, link_name)
                n_created = n_created + 1
            else:
                n_existed = n_existed + 1

            n_total[ID] = n_total[ID] + 1

    n_expected = len(input_IDs) * len(result_base_names)
    if verbose:
        print('{:5d} links created'.format(n_created))
        print('{:5d} links existed already'.format(n_existed))
        print('{:5d}/{} links available now'.format(n_created + n_existed, n_expected))
        n_tot = sum(n_total.values())
        print('{:5d} as cross-check'.format(n_tot))


def main(argv=None):

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    # Save calling command
    cfis.log_command(argv)
    if param.verbose:
        cfis.log_command(argv, name='sys.stdout')


    ### Start main program ###

    if param.verbose:
        print('Start of program {}'.format(os.path.basename(argv[0])))

    input_IDs = read_ID_list(param.input_IDs, verbose=param.verbose)

    result_base_names = ['psfex', 'psfex_interp_exp', 'setools_mask', 'setools_stat', 'setools_plot',
                         'pipeline_flag', 'final_cat', 'logs']

    if os.path.isdir(param.output_dir):
        if param.verbose:
            print('Directory {} already exists, continuing...'.format(param.output_dir))
    else:
        mkdir(param.output_dir)

    create_links(param.input_dir, param.output_dir, input_IDs, result_base_names, verbose=param.verbose)


    ### End main program

    if param.verbose:
        print('End of program {}'.format(os.path.basename(argv[0])))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
