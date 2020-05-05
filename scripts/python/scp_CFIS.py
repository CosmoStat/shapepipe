#!/usr/bin/env python

"""
:Name:        scp_CFIS.py

:Description:
    Secure-copy CFIS images.

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

import cfis



def get_copy_list(param, exclude_list, include_list=None, verbose=False):
    """Read directory and return list of files to copy accounting for files to exclude.
    """

    # Input file candidates to copy
    if include_list == None:
        # All files in directory
        files = glob.glob('*')
    else:
        # Read files from ascii file
        files = cfis.read_list(include_list)

    pattern = cfis.get_file_pattern(param.pattern, param.band, param.type)
    
    dst_list = []
    n_exc = 0
    for f in files:

        # Test if file matches pattern
        m = re.findall(pattern, f)
        if len(m) != 0:

            # Test if file is not on exclude list
            if cfis.exclude(m[0], exclude_list) == False:
                dst_list.append(m[0])
            elif verbose == True:
                n_exc += 1

    if verbose == True:
        print('Excluding {} files'.format(n_exc))

    if len(dst_list) == 0:
        raise cfis.CfisError('No matching files found')

    return dst_list


def scp(to_copy, t, remote, scp_cmd, to=True, dry_run=False, verbose=False):
    """Secure-copies files.

    Parameters
    ----------
    to_copy: list of strings
        list of file names to copy
    t: string
        type
    remote: string
        remote host name
    to: bool, optional, default=True
        If True (False), copy to (from) remote host
    dry_run: bool
        do not copy if True (default: False)
    verbose: bool
        verbose mode if true (defaul: False)
    Returns
    -------
    None
    """

    n_files = len(to_copy)

    if dry_run == True:
        sdry = ' (dry run)'
    else:
        sdry = ''

    if t == 'exposure':
        subdir = 'pitcairn'
    elif t in ('exposure_weight', 'exposure_weight.fz'):
        subdir = 'weights'
    elif t in ('exposure_flag', 'exposure_flag.fz'):
        subdir = 'flags'
    else:
        subdir = 'tiles'

    remote_dir = '{}:astro/data/CFIS/{}/'.format(remote, subdir)
    if to:
        arg = ' '.join(to_copy)
        dest = remote_dir
    else:
        to_copy.insert(0, ' ')
        arg  = ' {}'.format(remote_dir).join(to_copy)
        dest = '.'

    cfis.run_cmd('{} {} {}'.format(scp_cmd, arg, dest), verbose=verbose, run=not dry_run)
    #arg_list = arg.split(' ')
    #for tc in arg_list:
        #if tc != '':
            #cfis.run_cmd('{} {} {}'.format(scp_cmd, tc, dest), verbose=verbose, run=not dry_run)

    if verbose == True:
        print('scp-ed {} files{}'.format(n_files, sdry))



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
        band         = 'r',
        type         = 'tile',
        pattern      = '',
        remote       = 'cc2',
        scp_cmd      = 'scp',
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

    parser.add_option('-i', '--include_list', dest='include_list', type='string', default=None,
        help='Include files from this file, default=None, use all files in current directory')
    parser.add_option('-x', '--exclude_list', dest='exclude_list', type='string', default=None,
        help='Exclude files from this file, default=none (no excluded files)')

    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='type', type='string', default=p_def.type,
        help='data type, one in tiles (default)|cat|weight[.fz]|exposure|exposure_weight[.fz]|exposure_flag[.fz]]')
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern,
        help='file pattern to match, e.g.~\'^21\d{5}p\', default=none (=all match)')

    parser.add_option('-s', '--scp_cmd', dest='scp_cmd', type='string', default=p_def.scp_cmd,
        help='scp command, default=\'{}\''.format(p_def.scp_cmd))
    parser.add_option('-r', '--remote', dest='remote', type='string', default=p_def.remote,
        help='remote name, default=\'{}\''.format(p_def.remote))
    parser.add_option('', '--from', dest='fromR', action='store_true', default=False,
        help='Copy from remote (default: to remote)')
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
    p_def:  class cfis.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class cfis.param
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
    param.to = not param.fromR

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
    cfis.log_command(argv)
    if options.verbose:
        cfis.log_command(argv, name='sys.stderr')


    if options.verbose is True:
        print('Start program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    if param.exclude_list != None:
        exclude_list = cfis.read_list(param.exclude_list)
    else:
        exclude_list = []

    to_copy = get_copy_list(param, exclude_list, include_list=param.include_list, verbose=param.verbose)

    scp(to_copy, param.type, param.remote, param.scp_cmd, to=param.to, dry_run=param.dry_run, verbose=param.verbose)


    ### End main program

    if options.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

