#!/usr/bin/env python

"""Script create_weight_format.py

Unzips and removes first (empty) HDU of CFIS tile weight
such that they can be read by SExtractor,

:Authors: Martin Kilbinger, Axel Guinot

:Date: 24/03/2020

:Package: ShapePipe
"""

# Compability with python2.x for x>6
from __future__ import print_function

import sys
import os
import copy
from optparse import OptionParser
from astropy.io import fits

import cfis
from shapepipe.pipeline import file_io as io


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class tuff.param
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

    parser.add_option('-i', '--input', dest='input', type='string',
         help='input (fzipped) weight file name')
    parser.add_option('-o', '--output', dest='output', type='string',
         help='output (unzipped weight) file name, default=input name without \'.fz\' extension')
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

    if not options.input:
        print('input (fzipped weight image) file name not given (option \'-i\')')
        return False

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

    if param.output is None:
        param.output = os.path.splitext(param.input)[0]

    return param


def transform(input_fname, output_fname, verbose=False):
    """Transform input (fzipped weight) image.

    Parameters
    ----------
    input_fname: string
        input file name
    output_fname: string
        output file name
    verbose: bool, optional, default=False
        verbose output if True

    Return
    ------
    None
    """

    data = fits.getdata(input_fname, 1)
    header = fits.getheader(input_fname, 1)
    new_weight = io.FITSCatalog(output_fname, open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    new_weight.save_as_fits(data=data, image=True, image_header=header) 


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
    if param.verbose:
        cfis.log_command(argv, name='sys.stderr')

    transform(param.input, param.output, verbose=param.verbose)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

