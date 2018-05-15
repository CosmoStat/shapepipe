#!/usr/bin/env python

"""
:Name:          cfis_download_images.py

:Description:
    Downloads CFIS images from the CADC server, using
    vos tools.

:Author:        Martin Kilbinger <martin.kilbinger@cea.fr>

:Date:           2017
"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import copy
import re

import numpy as np
import pylab as plt

from astropy.io import ascii
from astropy.table import Table, Column

from optparse import OptionParser
from optparse import OptionGroup

import cfis
import stuff



def filter_file_list(in_file_list, pattern, band, type, verbose=False):
    """Return filtered version of file list according to pattern, band, and file type.

    Parameters
    ----------
    in_file_list: list of strings
        input file list
    pattern: string
        base pattern to match
    band: string
        band, one in 'u', 'r'
    type: string
        image type, one in 'tile', 'weight', 'weight.fz', 'exposure', 'exposure_flag',
	'exposure_flag.fz', 'exposure_weight', or 'exposure_weight.fz'
    verbose: bool, optional
        verbose mode if True; default=False

    Returns
    -------
    dst_list: list of strings
        final destination file list
    """

    f = open(in_file_list, 'rU')
    lines = f.readlines()
    f.close()

    dst_list = []

    pattern = cfis.get_file_pattern(pattern, band, type)

    if verbose == True:
        print('Filtering for pattern \'{}\''.format(pattern))

    for file in lines:

        m = re.findall(pattern, file)
        if len(m) != 0:
            dst_list.append(m[0])

    return dst_list


def remove_exclude(dst_list, exclude_list, verbose=False):
    """Remove files from dst_list that are on exclude_list and
       returned cleaned list.

    Parameters
    ----------
    dst_list: list of strings
        destination file list
    exclude_list: list of strings
        list of files to be excluded
    verbose: bool
        verbose mode if True; default=False

    Returns
    -------
    out_list: list of strings
        file list stripped of excluded files
    """

    if len(exclude_list) == 0:
        return dst_list

    out_list = []
    nexc     = 0
    for f in dst_list:
        if cfis.exclude(f, exclude_list) == False:
            out_list.append(f)
        else:
            nexc += 1

    if verbose == True:
        print('{} files excluded'.format(nexc))

    return out_list



def download(dst_list, out_dir, t, dry_run=False, verbose=False):
    """Download files using vos command.
    Parameters
    ----------
    dst_list: list of string
        file list to download
    out_dir: string
        output destination directory
    t: string
        file type, one in 'tile', 'weight', 'weight.fz', 'exposure',
	'exposure_flag', 'exposure_flag.fz', 'exposure_weight', or 'exposure_weight.fz'
    dry_run: bool, optional
        If True do not download but perform dry run; default=False
    verbose: bool, optional
        Verbose mode if True; default=False
    """

    logname = 'logfile.txt'
    if os.path.isfile(logname):
        os.remove(logname)

    if t == 'exposure':
        subdir = 'pitcairn'
    elif t in ('exposure_flag', 'exposure_flag.fz'):
	subdir = 'flags'
    elif t in ('exposure_weight', 'exposure_weight.fz'):
	subdir = 'weights'
    else:
        subdir = 'tile'

    sdry = ''
    if dry_run == True:
        sdry = ' (dry run)'

    n_dl = 0
    n_ex = 0
    for i, dst in enumerate(dst_list):

        dest = '{}/{}'.format(out_dir, dst)
        if not os.path.isfile(dest):
            print('Downloading {}{}'.format(dst_list[i], sdry))
            stuff.run_cmd('vos.sh vcp vos:cfis/{}/{} {}'.format(subdir, dst_list[i], out_dir),
                            verbose=False, run=not dry_run)
            n_dl += 1
        else:
            print('File {} exists'.format(dest))
	    n_ex += 1

    if verbose == True:
        print('Downloaded {} files, skipped {} files{}'.format(n_dl, n_ex, sdry))



def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class param
        parameter values
    """

    p_def = stuff.param(
        in_file_list = 'ls.txt',
        out_dir      = '.',
	    band         = 'r',
        type         = 'tile',
        pattern      = '',
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
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    parser.add_option('-i', '--in_file_list', dest='in_file_list', default=p_def.in_file_list, type='string',
         help='input list with files to download, default = {0}'.format(p_def.in_file_list))
    parser.add_option('-o', '--output_dir', dest='out_dir', default=p_def.out_dir, type='string',
         help='output directory name, default = {0}'.format(p_def.out_dir))

    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='type', type='string', default=p_def.type,
        help='data type, one of \'tile\' (default)| \'cat\'|\'weight\'|\'weight.fz\'|\'exposure\'|'
	     '\'exposure_flag\'|\'exposure_flag.fz\'|\'exposure_weight|\'exposure_weight.fz\'')
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern,
        help='file pattern to match, e.g.~\'^21\d{5}p\', default=none (=all match)')
    parser.add_option('-x', '--exclude_list', dest='exclude_list', type='string', default=None,
        help='Exclude files from this file, default=None')


    parser.add_option('-n', '--dry-run', dest='dry_run', action='store_true', default=False,
        help='dry run, only print commands')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='verbose')

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

    if options.band != 'r' and options.band != 'u':
        stuff.error('Invalid band \'{}\''.format(options.band))

    see_help = 'See option \'-h\' for help.'

    return True



def update_param(p_def, options):
    """Return default parameter, updated according to options.
    
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

    return param



def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

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

    dst_list     = filter_file_list(param.in_file_list, param.pattern, param.band, param.type, verbose=param.verbose)
    if param.exclude_list is not None:
        exclude_list = cfis.read_list(param.exclude_list)
    else:
        exclude_list = []

    dst_list = remove_exclude(dst_list, exclude_list, verbose=param.verbose)

    stuff.mkdir_p(param.out_dir)

    download(dst_list, param.out_dir, param.type, dry_run=param.dry_run, verbose=param.verbose)


    ### End main program

    if options.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

