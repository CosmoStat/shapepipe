#!/usr/bin/env python

"""Script stats_global.py

Output global stats from a ShapePipe run.

:Author: Martin Kilbinger

:Date: 07/2020
"""

import re
import os
import sys
import copy
import io
import glob

import numpy as np
import matplotlib.pylab as plt

from optparse import OptionParser

import cfis


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
        input_dir = '.',
        pattern = 'star_stat-',
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
    parser.add_option('-i', '--input_dir', dest='input_dir', type='string', \
         default=p_def.input_dir, \
         help='input directory, default=\'{}\''.format(p_def.input_dir))
    parser.add_option('-p', '--pattern', dest='pattern', type='string', \
         default=p_def.pattern, \
         help='input file pattern, default=\'{}\''.format(p_def.pattern))

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


def get_stats_file_list(param):

    # Full paths
    paths = glob.glob('{}/output/*/setools_runner/output/stat/{}*'.
                      format(param.input_dir, param.pattern))

    # File names w/o directories
    names = []
    for f in paths:
        names.append(os.path.basename(f))

    # Get unique files
    # TODO: Check whether multiple files are identical
    names_unq = []
    paths_unq = []
    for i in range(len(names)):
        if names[i] not in names_unq:
            names_unq.append(names[i])
            paths_unq.append(paths[i])

    if param.verbose:
        print('{} files in total, {} unique files found'.
              format(len(paths), len(paths_unq)))

    return paths_unq


def gather_values(paths, verbose=False):
    """Return values found in files

    Parameters
    ----------
    paths: list of string
        file names
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    values: dictionary with lists of float
        values and their keys
    """

    values = {}
    for path in paths:
        with open(path) as f:
            lines = f.readlines()
        for line in lines: 
            m = re.search("#", line)
            if m:
                continue
            m = re.search('(.*) = (\S*)', line)
            if m:
                key = m[1]
                val = m[2]
                if not key in values:
                    values[key] = []
                values[key].append(float(val))

    if verbose:
        print('{} keys created'.format(len(values)))
        for key in values:
            print('#{{values[{}]}} = {}'.format(key, len(values[key])))

    return values


def compute_histograms(values, verbose=False):
    """Compute histograms from dictionary of value lists.

    Parameters
    ----------
    values: dictionary with lists of float
        values and their keys
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    hists: dictionary of histograms
        histograms for all keys
    """

    hists = {}
    for key in values:

        hists[key] = np.histogram(values[key], bins=50)

    return hists


def plot_histograms(hists, verbose=False):
    """Create histogram plots.

    Parameters
    ----------
    hists: dictionary of histograms
        histograms for all keys
    verbose: bool, optional, default=False
        verbose output if True
    """

    fig, (ax) = plt.subplots()

    i = 0
    for key in hists:
        bins = hists[key][1][:-1]
        freq = hists[key][0]

        if verbose:
            print(key, bins, freq)

        ax = plt.gca()
        width = bins[1] = bins[0]
        ax.bar(bins, freq, width=width)

        xmin = min(bins)
        xmax = max(bins)
        plt.xlim(xmin, xmax)
        plt.xlabel(key)
        plt.ylabel('frequency')

        plot_name = 'hist_{}.png'.format(i)
        plt.savefig(plot_name)

        i = i + 1


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

    files = get_stats_file_list(param)

    values = gather_values(files, verbose=param.verbose)

    hists = compute_histograms(values, verbose=param.verbose)

    plot_histograms(hists, verbose=param.verbose)

    ### End main program

    if param.verbose:
        print('End of program {}'.format(os.path.basename(argv[0])))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

