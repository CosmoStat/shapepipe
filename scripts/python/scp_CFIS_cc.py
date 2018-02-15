#!/usr/bin/env python

"""
:Name:        scp_CFIS_cc.py

:Description:
    Secure-copies CFIS images to and from cc@in2p3.

:Author:      Martin Kilbinger <martin.kilbinger@cea.fr>

:Date:        2017
"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import copy
import glob
import re

import numpy as np
import pylab as plt

from astropy.io import ascii
from astropy.table import Table, Column

from optparse import OptionParser
from optparse import OptionGroup

import stuff
import cfis



def exclude(f, exclude_list):
    """Return True if f is on exclude_list

    Parameters
    ----------
    f: string
        file name
    exclude_list: list of strings
        list of files

    Returns
    -------
    is_in_exclude: bool
        True (False) if f is in list
    """

    return f in exclude_list


def get_copy_list(param, exclude_list, include_list=None, verbose=False):
    """Read directory and return list of files to copy accounting for files to exclude.
    """

    # Input file candidates to copy
    if include_list == None:
        # All files in directory
        files = glob.glob('*')
    else:
        # Read files from ascii file
        files = CFIS.read_list(include_list)

    pattern = CFIS.get_file_pattern(param.pattern, param.band, param.type)
    
    dst_list = []
    n_exc = 0
    for f in files:

        # Test if file matches pattern
        m = re.findall(pattern, f)
        if len(m) != 0:

            # Test if file is not on exclude list
            if CFIS.exclude(m[0], exclude_list) == False:
                dst_list.append(m[0])
            elif verbose == True:
                n_exc += 1

    if verbose == True:
        print('Excluding {} files'.format(n_exc))

    return dst_list


def scp(to_copy, t, to_cc=True, dry_run=False, verbose=False):
    """Secure-copies files to cc.

    Parameters
    ----------
    to_copy: list of strings
        list of file names to copy
    t: string
        type
    to_cc: bool, optional, default=True
        If True (False), copy to (from) cc
    dry_run: bool
        do not copy if True (default: False)
    verbose: bool
        verbose mode if true (defaul: False)
    Returns
    -------
    None
    """

    if dry_run == True:
        sdry = ' (dry run)'
    else:
        sdry = ''
        print('Password for cc: Kidt9uslYon')

    if t == 'raw':
        subdir = 'pitcairn'
    else:
        subdir = 'tiles'

    cc_dir = 'cc:/sps/euclid/Users/mkilbing/astro/data/CFIS/{}'.format(subdir)
    if to_cc:
        arg = ' '.join(to_copy)
        dest = cc_dir
    else:
        to_copy.insert(0, ' ')
        ' {}'.format(cc_dir).join(to_copy)
        dest = '.'
    stuff.run_cmd('scp {} {}'.format(arg, dest), verbose=verbose, run=not dry_run)

    if verbose == True:
        print('scp-ed {} files{}'.format(len(to_copy), sdry))



def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class stuff.param
        parameter values
    """

    p_def = stuff.param(
        band         = 'r',
        type         = 'tiles',
        pattern      = '',
    )

    return p_def



def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class stuff.param
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

    parser.add_option('-i', '--include_list', dest='include_list', type='string', default=None,
        help='Include files from this file, default=None, use all files in current directory')
    parser.add_option('-x', '--exclude_list', dest='exclude_list', type='string', default=None,
        help='Exclude files from this file, default=none (no excluded files)')

    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='type', type='string', default=p_def.type,
        help='data type, one of \'tiles\' (default)| \'cat\'|\'weight\'|\'raw\'')
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern,
        help='file pattern to match, e.g.~\'^21\d{5}p\', default=none (=all match)')

    parser.add_option('', '--from_cc', dest='from_cc', action='store_true', default=False,
        help='Copy from cc (default: to cc)')
    parser.add_option('-n', '--dry-run', dest='dry_run', action='store_true', default=False,
        help='dry run, only print commands')
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

    see_help = 'See option \'-h\' for help.'

    return True



def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class stuff.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class stuff.param
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
    param.to_cc = not param.from_cc

    return param



def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)


    # Save calling command
    stuff.log_command(argv)
    if options.verbose:
        stuff.log_command(argv, name='sys.stderr')


    if options.verbose is True:
        print('Start program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    if exclude_list != None:
        exclude_list = CFIS.read_list(param.exclude_list)
    else:
        exclude_list = []

    to_copy = get_copy_list(param, exclude_list, verbose=param.verbose)

    scp(to_copy, param.type, to_cc=param.to_cc, dry_run=param.dry_run, verbose=param.verbose)


    ### End main program

    if options.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

