#!/usr/bin/env python

"""
:Name:          cfis_download_images.py

:Description:
    Downloads CFIS images from the CADC server, using
    vos tools.

:Author:        Martin Kilbinger <martin.kilbinger@cea.fr>

:Date:          2017

:Package:       ShapePipe
"""

import sys
import os
import copy
import re

import numpy as np

from astropy.io import ascii
from astropy.table import Table, Column

from optparse import OptionParser
from optparse import OptionGroup

from shapepipe.utilities.canfar import vosHandler

import cfis


# Global VCP definition
vcp = vosHandler('vcp')


def filter_file_list(in_file_list, pattern, band, type, columns, in_number_only=False, verbose=False):
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
        image type, one in 'tile', 'weight', 'weight.fz', 'cat', 'exposure', 'exposure_flag',
	    'exposure_flag.fz', 'exposure_weight', or 'exposure_weight.fz'
    columns: list of int
        column numbers: name [size]
    in_number_only: bool, optional, default=False
        if True, input file names are image number only
    verbose: bool, optional
        verbose mode if True; default=False

    Returns
    -------
    dst_list: list of strings
        final destination file list
    size_list: list of int
        list of file sizes for files in *dst_list*
    """

    f = open(in_file_list, 'r')
    lines = f.readlines()
    f.close()

    dst_list = []
    size_list = []

    pattern = cfis.get_file_pattern(pattern, band, type)

    if verbose == True:
        print('Filtering for pattern \'{}\''.format(pattern))

    for line in lines:

        l = line.split()
        file_base = l[columns[0]]

        if in_number_only:
            file_name = cfis.get_file_pattern('{}p'.format(file_base), band, type, want_re=False)
        else:
            file_name = file_base

        m = re.findall(pattern, file_name)
        if len(m) != 0:
            dst_list.append(m[0])

            if len(columns) == 2:
                size_list.append(int(l[columns[1]]))

    return dst_list, size_list


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



def download(dst_list, size_list, out_dir, t, dry_run=False, verbose=False, quick=False,
             certfile=None, out_list=None):
    """Download files using vos command.
    Parameters
    ----------
    dst_list: list of string
        file list to download
    size_list: list of int
        list of file sizes for files in *dst_list*
    out_dir: string
        output destination directory
    t: string
        file type, one in 'tile', 'weight', 'weight.fz', 'exposure',
     'exposure_flag', 'exposure_flag.fz', 'exposure_weight', or 'exposure_weight.fz'
    dry_run: bool, optional
        if True do not download but perform dry run; default=False
    verbose: bool, optional, default=False
        verbose mode if True
    quick: bool, optional, default=False
        quick copy mode if True
    certfile: string, optional, default=None
        certificate file for vos identification
    out_list: string, optional, default=None
        output file name for files to download (if dry_run is True)
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
        subdir = 'tiles_DR2'

    sdry = ''
    if dry_run == True:
        sdry = ' (dry run)'
        if out_list:
            f_list = open(out_list, 'w')

    n_dl = 0
    n_ex = 0
    for i, dst in enumerate(dst_list):

        dest = '{}/{}'.format(out_dir, dst)

        do_download = True

        if os.path.isfile(dest):
            if len(size_list) > 0:
                size = os.path.getsize(dest)
                if size != size_list[i]:
                    do_download = True
                    print('File {} incomplete.'.format(dest), end=' ')
                else:
                    do_download = False
                    n_ex += 1
                    print('File {} exists and is complete'.format(dest))
            else:
                # No information about size given: do not overwrite
                # existing file
                print('File {} exists, information whether complete is missing, no download'.format(dest))
                do_download = False

        else:
            do_download = True

        if do_download:

            print('Downloading {}/{}{}'.format(subdir, dst_list[i], sdry))

            if out_list:
                if len(size_list) > 0:
                    print(size_list[i], end=' ', file=f_list)
                print(dst_list[i], file=f_list)

            cmd = 'vcp'
            src = 'vos:cfis/{}/{}'.format(subdir, dst_list[i])

            sys.argv = []
            sys.argv.append(cmd)
            if quick == True:
                sys.argv.append('--quick')
            if certfile:
                sys.argv.append('--certfile={}'.format(certfile))
            sys.argv.append(src)
            sys.argv.append(out_dir)

            vcp()

            n_dl += 1

    if verbose == True:
        print('Downloaded {} files, skipped {} files{}'.format(n_dl, n_ex, sdry))

    if out_list:
        f_list.close()



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

    p_def = cfis.param(
        in_file_list = 'ls.txt',
        out_dir      = '.',
	    band         = 'r',
        type         = 'tile',
        pattern      = '',
        scolumns     = '0',
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

    # Input files
    parser.add_option('-i', '--in_file_list', dest='in_file_list', default=p_def.in_file_list, type='string',
         help='input list with files to download, default = {0}'.format(p_def.in_file_list))
    parser.add_option('-x', '--exclude_list', dest='exclude_list', type='string', default=None,
        help='Exclude files from this file, default=None')

    # Output files
    parser.add_option('-o', '--output_dir', dest='out_dir', default=p_def.out_dir, type='string',
         help='output directory name, default = {0}'.format(p_def.out_dir))
    parser.add_option('-l', '--out_list', dest='out_list', type='string', default=None,
        help='Output file for files to download (only with -n)')

    # Other options
    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='type', type='string', default=p_def.type,
        help='data type, one of \'tile\' (default)| \'cat\'|\'weight\'|\'weight.fz\'|\'exposure\'|'
	     '\'exposure_flag\'|\'exposure_flag.fz\'|\'exposure_weight|\'exposure_weight.fz\'')
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern,
        help='file pattern to match, e.g.~\'^21\d{5}p\', default=none (=all match)')
    parser.add_option('-c', '--in_file_columns', dest='scolumns', type='string', default=p_def.scolumns,
        help='input file columns: name [size], default={}'.format(p_def.scolumns))
    parser.add_option('', '--in_number_only', dest='in_number_only', action='store_true',
        help='input file names are image number only')

    # vcp options
    parser.add_option('', '--certfile', dest='certfile', type='string',
        help='certificate file')
    parser.add_option('-q', '--quick', dest='quick', action='store_true', default=False, help='quick copy mode')

    # Misc options
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
        print('Invalid band \'{}\''.format(options.band))
        return False

    if options.out_list and not options.dry_run:
        print('Ouput list (-l) only valid if dry run (-n)')
        return False

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

    tmp = cfis.my_string_split(param.scolumns, stop=True)
    param.columns = [int(c) for c in tmp]

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
    cfis.log_command(argv)
    if options.verbose:
        cfis.log_command(argv, name='sys.stdout')


    if options.verbose is True:
        print('Start program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    dst_list, size_list = filter_file_list(param.in_file_list, param.pattern, param.band, param.type, \
                                           param.columns, in_number_only=param.in_number_only, \
                                           verbose=param.verbose)

    if param.exclude_list is not None:
        exclude_list = cfis.read_list(param.exclude_list)
    else:
        exclude_list = []

    dst_list = remove_exclude(dst_list, exclude_list, verbose=param.verbose)

    cfis.mkdir_p(param.out_dir)

    download(dst_list, size_list, param.out_dir, param.type, dry_run=param.dry_run,
             verbose=param.verbose, quick=param.quick, certfile=param.certfile, out_list=param.out_list)


    ### End main program

    if options.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))
