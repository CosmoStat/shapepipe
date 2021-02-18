#!/usr/bin/env python

"""Script canfar_avail_results.py

Check whether results files are available on vos.

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

from shapepipe.utilities.canfar import vosHandler

import cfis


# Global VLS definition
vls = vosHandler('vls')


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
        input_IDs  = '.',
        input_vos = 'cosmostat/kilbinger/results',
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
    parser.add_option('-i', '--input_IDs', dest='input_IDs', type='string', default=p_def.input_IDs,
         help='input tile ID file(s) or directory path, default=\'{}\''.format(p_def.input_IDs))
    parser.add_option('', '--input_vos', dest='input_vos', type='string', default=p_def.input_vos,
         help='input path on vos, default=\'{}\''.format(p_def.input_vos))
    parser.add_option('-o', '--output_not_avail', dest='output_not_avail', type='string',
         help='output file for not-available IDs, default no output')

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


def read_input_files(input_path, verbose=False):
    """Return list of ID files.

    Parameters
    ----------
    input_path: string
        list of files or directory name
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    ID_files: list of strings
        file names with tile IDs
    """

    if os.path.isdir(input_path):
        input_files = glob.glob('{}/*'.format(input_path))
    else:
        input_files =  cfis.my_string_split(input_path, stop=True, sep=' ')

    ID_files = []
    for f in input_files:
        if os.path.isdir(f):
            if verbose:
                print('Skipping directory \'{}\''.format(f))
        else:
            ID_files.append(f)

    if verbose:
        print('{} input files found'.format(len(ID_files)))

    return ID_files


def check_results(ID_files, input_vos, result_base_names, n_complete, verbose=False):
    """Count the number of result files uploaded to vos for each input ID file.

    Parameters
    ----------
    ID_files: list of strings
        file name with tile IDs
    input_vos: string
        vos input directory
    result_base_names: list of strings
        result file base names
    n_complete: int
        number of files for complete result set
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    n_found: dictionary
        number of files found for each input file and each ID
    n_IDs: dictionary
        number of ID for each input file
    IDs_not_avail: list
        IDs that are not available on vos
    """

    cmd = 'vls'
    vos_dir = 'vos:cfis/{}'.format(input_vos)

    sys.argv = []
    sys.argv.append(cmd)
    sys.argv.append(vos_dir)
    #sys.argv.append('> .vls.tmp')
    #if certfile:
        #sys.argv.append('--certfile={}'.format(certfile))

    if verbose:
        print('Getting vos directory content from vls...')
    f = io.StringIO()

    try:
        with redirect_stdout(f):
            vls()
    except:
        print('Error during vls command')
        raise
 
    vls_out = f.getvalue()

    n_found = {}
    n_IDs = {}
    IDs_not_avail = []

    # Loop over all input files
    for ID_list in ID_files:
        with open(ID_list) as f:
            if verbose:
                print('Checking ID list file {}...'.format(ID_list))
            n_found[ID_list] = {}
            n_IDs[ID_list] = 0

            # Loop over all lines = IDs in file
            for line in f:
                ID = line.rstrip()
                n_found[ID_list][ID] = 0

                # Count how many result files are available
                for base in result_base_names:
                    name = '{}_{}.tgz'.format(base, ID)
                    if name in vls_out:
                        n_found[ID_list][ID] = n_found[ID_list][ID] + 1
                n_IDs[ID_list] = n_IDs[ID_list] + 1

                # If not complete set found, add to not-avail list
                if n_found[ID_list][ID] != n_complete:
                    IDs_not_avail.append(ID)

    return n_found, n_IDs, IDs_not_avail


def output_summary(n_found, n_IDs, n_complete):
    """Create output with summary of result availability.

    Parameters
    ----------
    n_found: dictionary
        number of files found for each input file and each ID
    n_IDs: dictionary
        number of ID for each input file
    n_complete: int
        number of files of a complete set of results
    """

    for ID_list in n_found.keys():
        nf = sum(value == n_complete for value in n_found[ID_list].values())
        print('{}: {}/{} ({:.1f}%) complete'.format(os.path.basename(ID_list), nf, n_IDs[ID_list],  nf/n_IDs[ID_list]*100))


def output_IDs(ID_list, output):
    """Write IDs to file

    Parameters
    ----------
    ID_list: list of string
        IDs
    output: string
        output file name
    """

    f = open(output, 'w')
    for ID in ID_list:
        print(ID, file=f)
    f.close()


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

    ID_files = read_input_files(param.input_IDs, verbose=param.verbose)

    result_base_names = ['psfex', 'psfex_interp_exp', 'setools_mask', 'setools_stat', 'setools_plot',
                         'final_cat', 'pipeline_flag', 'logs']

    n_complete = len(result_base_names)

    n_found, n_IDs, IDs_not_avail = check_results(ID_files, param.input_vos, result_base_names, n_complete, verbose=param.verbose)

    output_summary(n_found, n_IDs, n_complete)

    if param.output_not_avail:
        output_IDs(IDs_not_avail, param.output_not_avail)


    ### End main program

    if param.verbose:
        print('End of program {}'.format(os.path.basename(argv[0])))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
