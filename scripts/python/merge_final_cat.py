#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script merge_final_cat.py

Merge all final catalogues, created by ShapePipe module 'make_catalogue_runner',
into a joined numpy binary file.

:Authors: Axel Guinot, Martin Kilbinger

:Date: 2020

:Package: ShapePipe
"""

from astropy.io import fits
import numpy as np
import os
import sys
import re
import copy

from optparse import OptionParser

from tqdm import tqdm


class param:
    """General class to store (default) variables
    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)


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

    p_def = param(
        input_path  = '.',
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

    parser.add_option('-i', '--input_path', dest='input_path', type='string',
                      default=p_def.input_path,
                      help='input path, default=\'{}\''.format(p_def.input_path))
    parser.add_option('-p', '--param_path', dest='param_path', type='string',
                      default=None,
                      help='parameter file path, default=None')

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


def read_param_file(path, verbose=False):
    """Return parameter list read from file

    Parameters
    ----------
    path: string
        input file name
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    param_list: list of strings
        parameter names
    """

    param_list = []

    if path:

        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue 
                entry = line.rstrip()
                if not entry or entry == '':
                    continue
                param_list.append(entry)

    if verbose:
        if len(param_list) > 0: 
            print('Copying {} columns'.format(len(param_list)), end='')
        else:
            print('Copying all columns', end='')
        print(' into final catalogue')

    # Check for multiples
    multiples = []
    for param in param_list:
        if param_list.count(param) > 1:
            multiples.append(param)

    if len(multiples) > 0:
        print('The following parameters are more than one times '
              'in the parameter file: ', end='')
        for m in multiples:
            print(m, end=' ')
        print()
        raise ValueError('Multiple identical keys found')

    return param_list
                            

def get_data(path, hdu_num, param_list):
    """Return data of selected columns from FITS file.

    Parameters
    ----------
    path: string
        input file name
    hdu_num: int
        HDU number
    param_list: list of strings
        parameters to be extracted. If none, copy
        all columns

    Returns
    -------
    data: numpy array
        data columns
    """

    hdu_list = fits.open(path)
    hdu = hdu_list[hdu_num]

    if param_list:
        cols = []
        for p in param_list:
            cols.append(hdu.columns[p]) 
        coldefs = fits.ColDefs(cols)
        hdu_new = fits.BinTableHDU.from_columns(coldefs)
        d = hdu_new.data
    else:
        d = hdu.data

    return d


def main(argv=None):

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    path = param.input_path

    param.param_list = read_param_file(param.param_path, verbose=param.verbose)

    l = os.listdir(path=path)
    lpath = []
    for this_l in l:
        lpath.append(os.path.join(path, this_l))

    if param.verbose:
        print('{} final catalog files found'.format(len(lpath)))

    # Determine number of columns and keys
    d_tmp = get_data(lpath[0], 1, param.param_list)

    d = np.zeros(d_tmp.shape, dtype=d_tmp.dtype)
    for key in d_tmp.dtype.names:
        d[key] = d_tmp[key]

    #new_dt = np.dtype(d_tmp.dtype.descr + [('TILE_ID', '>i4')])
    #d = np.zeros(d_tmp.shape, dtype=new_dt)

    #if 'TILE_ID' in d_tmp.dtype.names:
        #d['TILE_ID'].fill(int(''.join(re.findall('\d+', l[0]))))

    # Read all final catalogues and merge
    count = 1
    #for i in tqdm(lpath[1:], total=len(lpath)-1):
    for i in lpath[1:]:
        if ('final_cat' not in i) | ('.npy' in i):
            continue

        try:
            d_tmp = get_data(i, 1, param.param_list)

            dd = np.zeros(d_tmp.shape, dtype=d_tmp.dtype)
            for key in d_tmp.dtype.names:
                dd[key] = d_tmp[key]


            #if 'TILE_ID' in d_tmp.dtype.names:
                #dd['TILE_ID'].fill(int(''.join(re.findall('\d+', i))))

            d = np.concatenate((d, dd))
            count = count + 1
            print('File \'{}\' copied ({}/{})'.format(i, count, len(lpath)))
        except:
            print('Error while copying file \'{}\''.format(i))
            #raise

    # Save merged catalogue as numpy binary file
    print('Saving final np cat')
    np.save('final_cat.npy', d)
    print('Done')

    if param.verbose:
        print('{} catalog files merged with success'.format(count))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
