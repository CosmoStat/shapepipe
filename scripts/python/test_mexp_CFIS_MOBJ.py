#!/usr/bin/env python

"""Script cfis_select_tileobj_expPSF.py

  Print histograms of object ID's in object multi-exposure files.
  Expected output per file:
   (name, 0)
   0 count_0
   1 count_1
   ...
  For a large area selection, the highest counts should occur for n=3.
  For smaller areas, the number of smaller counts is more frequent.

:Authors: Martin Kilbinger

:Date: 11/12/2018
"""

# Compability with python2.x for x>6
from __future__ import print_function

import glob
import os
import sys
import copy

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup
import numpy as np
from astropy.io import fits

from generic import stuff


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

    base_dir = '{}/data'.format(os.environ['HOME'])

    p_def = stuff.param(
        input_dir = '{}/tiles'.format(base_dir),
        pattern   = 'CFIS_MOBJ',
        key       = 'ID',
        hdu       = 1,
        mode      = 'count_unique',
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
    parser = OptionParser(usage=usage, formatter=stuff.IndentedHelpFormatterWithNL())

    # Input
    parser.add_option('-i', '--input_dir', dest='input_dir', type='string', default=p_def.input_dir,
         help='input directory, default=\'{}\''.format(p_def.input_dir))
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern,
         help='input file pattern, default=\'{}\''.format(p_def.pattern))
    parser.add_option('-k', '--key', dest='key', type='string', default=p_def.key,
         help='FITS column key, default=\'{}\''.format(p_def.key))
    parser.add_option('-m', '--mode', dest='mode', type='string', default=p_def.mode,
         help='output mode, default=\'{}\''.format(p_def.mode))
    parser.add_option('', '--hdu', dest='hdu', type='int', default=p_def.hdu,
         help='FITS hdu, default=\'{}\''.format(p_def.hdu))

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
    if param.verbose:
        stuff.log_command(argv, name='sys.stderr')

    if param.verbose is True:
        print('Start of program {}'.format(os.path.basename(argv[0])))

    ### Start main program ###

    obj_list = glob.glob('{}/{}*'.format(param.input_dir, param.pattern))

    n_zero_max = 3

    for obj in obj_list:
        f = fits.open(obj)
        u, i = np.unique(f[param.hdu].data[param.key], return_counts=True)
        print(obj, len(f[param.hdu].data))

        if param.mode == 'count_unique':
            j = 0
            n_zero = 0
            while True:
                n = len(u[i == j])
                print(j, n)
                if n == 0:
                    n_zero += 1
                if n_zero == n_zero_max:
                    break
                j = j + 1

        elif param.mode == 'hist':
            for uu, ii in zip(u, i):
                print(uu, ii)

            

    ### End main program ###

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))




