#!/usr/bin/env python

"""Script stats_global.py

Output global stats from a ShapePipe run.

:Author: Martin Kilbinger

:Date: 07/2020

:Package: ShapePipe
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
        output_dir = '.',
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
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', \
         default=p_def.output_dir, \
         help='output directory, default=\'{}\''.format(p_def.output_dir))
    parser.add_option('-c', '--config', dest='config', type='string', \
         help='configuration file, default=none')
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
    options: tuple
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
    """Return list of files with statistics information

    Parameters
    ----------
    param: class param
        parameter values

    Returns
    -------
    paths_unq: list of string
        unique file name list
    """

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
                val = float(m[2])
                if not np.isnan(val):
                    if key not in values:
                        values[key] = []
                    values[key].append(val)
                else:
                    if verbose:
                        print('NaN found in file \'{}\', key \'{}\''
                              .format(path, key))

    if verbose:
        print('{} keys created'.format(len(values)))
        for key in values:
            print('#{{values[{}]}} = {}'.format(key, len(values[key])))

    return values


def compute_histograms(values, config=None, verbose=False):
    """Compute histograms from dictionary of value lists.

    Parameters
    ----------
    values: dictionary with lists of float
        values and their keys
    config: class config
        configuration values
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    hists: dictionary of histograms
        histograms for all keys
    """

    if config and config.has_option('ALL', 'nbins'):
        nbins_global = config.getint('ALL', 'nbins')
    else:
        nbins_global = 50

    print('nbins_global = {}'.format(nbins_global))

    hists = {}
    i = 0
    for key in values:

        si = str(i)

        if config and config.has_option(si, 'nbins'):
            nbins = config.getint(si, 'nbins')
        else:
            nbins = nbins_global

        try:
            hists[key] = np.histogram(values[key], bins=nbins)
        except ValueError as err:
            print('Skipping histogram #{}: {}'.format(i, err))

        i = i + 1

    return hists


def plot_histograms(hists, config=None, output_dir='.', verbose=False):
    """Create histogram plots.

    Parameters
    ----------
    hists : dictionary of histograms
        histograms for all keys
    config : class config
        configuration values
    output_dir : string, optional, default='.'
        output directory
    verbose : bool, optional, default=False
        verbose output if True
    """

    fig, (ax) = plt.subplots()
    plt.tight_layout()

    if config and config.has_option('ALL', 'fontsize'):
        fontsize = config.getint('ALL', 'fontsize')
        plt.rcParams.update({'font.size': fontsize})

    xlim_fac = 0.05

    i = 0
    for key in hists:

        si = str(i)

        if config and \
            config.has_option(si, 'plot') and \
            config.getboolean(si, 'plot') == False:
            if verbose:
                print('Skipping histogram #{} {}'.format(i, key))

        else:

            bins = hists[key][1][:-1]
            freq = hists[key][0]

            fig = plt.figure()
            ax = plt.subplot(111)
            width = bins[1] - bins[0]
            ax.bar(bins, freq, width=width)

            # Overwrite limits if found in config file.
            # If not stretch bin boundaries
            if config and config.has_option(si, 'xmin'):
                xmin = config.getfloat(si, 'xmin')
            else:
                xmin = min(bins)
            if config and config.has_option(si, 'xmax'):
                xmax = config.getfloat(si, 'xmax')
            else:
                xmax = max(bins)

            dx = xmax - xmin

            # If limits from bin boundaries: extend by small
            # amount
            if config and config.has_option(si, 'xmin'):
                xxmin = xmin
            else:
                xxmin = xmin - dx * xlim_fac
            if config and config.has_option(si, 'xmax'):
                xxmax = xmax
            else:
                xxmax = xmax + dx * xlim_fac

            plt.xlim(xxmin, xxmax)

            if config and config.has_option(si, 'xlabel'):
                xlabel = config.get(si, 'xlabel')
            else:
                xlabel = key
            plt.xlabel(xlabel)

            ymin = min(freq)
            ymax = max(freq)
            plt.ylim(ymin, ymax)
            ylabel = 'frequency'
            plt.ylabel(ylabel)

            if config and config.has_option(si, 'title'):
                plt.title(config.get(si, 'title'))

            if config and config.has_option(si, 'fname'):
                file_base = config.getexpanded(si, 'fname')
            else:
                file_base = 'hist_{}'.format(i)

            if verbose:
                print('Creating files \'{}.*\''.format(file_base))

            plt.savefig('{}/{}.png'.format(output_dir, file_base), bbox_inches='tight')
            np.savetxt('{}/{}.txt'.format(output_dir, file_base), np.transpose([bins, freq]),
                       fmt='%10g', header='[{}] [{}]'.format(xlabel, ylabel)),

        i = i + 1



def get_config(config_path, verbose=False):
    """Return configuration file values.

    Parameters
    ----------
    config_path : string
        configuration file path
    verbose : bool, optional, default=False
        verbose output if True

    Returns
    -------
    conf : CustomParser
        configuration values, None if config_path does not exist
    """

    if config_path is None:
        return None

    from shapepipe.pipeline.config import CustomParser

    if verbose:
        print('Reading configuration file \'{}\''.format(config_path))

    if not os.path.exists(config_path):
        raise OSError('Configuration file \'{}\' does not exist'.format(config_path))

    conf = CustomParser()
    conf.read(config_path)

    return conf


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

    config = get_config(param.config, verbose=param.verbose)

    hists = compute_histograms(values, config=config, verbose=param.verbose)

    if os.path.isfile(param.output_dir):
        raise OSError('Output path \'{}\' is a regular file'
                      ''.format(param.output_dir))
    if not os.path.isdir(param.output_dir):
        os.mkdir(param.output_dir)

    plot_histograms(hists, config=config, output_dir=param.output_dir,
                    verbose=param.verbose)

    ### End main program

    if param.verbose:
        print('End of program {}'.format(os.path.basename(argv[0])))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

