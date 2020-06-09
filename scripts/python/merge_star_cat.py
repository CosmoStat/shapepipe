#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script merge_psf_validation.py

Merge single-exposure single-CCD validation PSF files.
for testing and validating the PSF model.

:Authors: Axel Guinot, Martin Kilbinger

:Date: 2019, 2020

:Package: ShapePipe
"""


from shapepipe.pipeline.execute import execute

import os
import sys
import re
import copy

import numpy as np
from optparse import OptionParser
from astropy.io import fits

import cfis
from shapepipe.pipeline import file_io as sc


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
        input_file_pattern = 'validation_psf',
        output_path = './psf_cat_full.fits'
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

    parser.add_option('-i', '--input_dir', dest='input_dir', type='string',
         help='input directory')
    parser.add_option('-p', '--input_file_pattern', dest='input_file_pattern', type='string',
         default=p_def.input_file_pattern, help='input file pattern, default=\'{}\''.format(p_def.input_file_pattern))
    parser.add_option('-o', '--output_path', dest='output_path', type='string', default=p_def.output_path,
         help='output file name, default=\'{}\''.format(p_def.output_path))
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

    if not options.input_dir:
        print('input directory not given (option \'-i\')')
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

    return param


def create_full_cat(input_dir, input_file_pattern, output_path, verbose=False):
    """Create full PSF catalogue from input PSF files.

    Parameters
    ----------
    input_dir: string
        input directory
    input_file_pattern: string
        input file pattern
    output_path: string
        output file path
    verbose: bool, optional, default=False

    Returns
    -------
    None
    """

    # Get input file names
    fnames = os.listdir(input_dir)
    starcat_names = ['{}/{}'.format(input_dir, name) for name in fnames if input_file_pattern in name]
    if len(starcat_names) == 0:
        raise IndexError('Input file list is empty')

    if verbose:
        print('Found {} files'.format(len(starcat_names)))

    # Initialise columns to be saved
    x, y, ra, dec = [], [], [], []
    g1_psf, g2_psf, size_psf = [], [], []
    g1, g2, size = [], [], []
    flag_psf, flag_star = [], []
    mag, snr, psfex_acc = [], [], []
    ccd_nb = []

    if verbose:
        print('Reading input FITS PSF catalogues...')
    # Read input files and add to output columns
    for name in starcat_names:
        starcat_j = fits.open(name)

        # positions
        x += list(starcat_j[2].data['X'])
        y += list(starcat_j[2].data['Y'])
        ra += list(starcat_j[2].data['RA'])
        dec += list(starcat_j[2].data['DEC'])

        # shapes (convert sigmas to R^2)
        g1_psf += list(starcat_j[2].data['E1_PSF_HSM'])
        g2_psf += list(starcat_j[2].data['E2_PSF_HSM'])
        size_psf += list(starcat_j[2].data['SIGMA_PSF_HSM']**2)
        g1 += list(starcat_j[2].data['E1_STAR_HSM'])
        g2 += list(starcat_j[2].data['E2_STAR_HSM'])
        size += list(starcat_j[2].data['SIGMA_STAR_HSM']**2)

        # flags
        flag_psf += list(starcat_j[2].data['FLAG_PSF_HSM'])
        flag_star += list(starcat_j[2].data['FLAG_STAR_HSM'])

        # misc
        mag += list(starcat_j[2].data['MAG'])
        snr += list(starcat_j[2].data['SNR'])
        psfex_acc += list(starcat_j[2].data['ACCEPTED'])

        # CCD number
        ccd_nb += [re.split(r"\-([0-9]*)\-([0-9]+)\.", name)[-2]]*len(starcat_j[2].data['RA'])

    # Prepare output FITS catalogue
    output = sc.FITSCatalog(output_path,
                            open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                            SEx_catalog=True)

    # Collect columns
    # convert back to sigma for consistency
    data = {'X': x, 'Y': y, 'RA': ra, 'DEC': dec,
            'E1_PSF_HSM': g1_psf, 'E2_PSF_HSM': g2_psf, 'SIGMA_PSF_HSM': np.sqrt(size_psf),
            'E1_STAR_HSM': g1, 'E2_STAR_HSM': g2, 'SIGMA_STAR_HSM': np.sqrt(size),
            'FLAG_PSF_HSM': flag_psf, 'FLAG_STAR_HSM': flag_star,
            'MAG': mag, 'SNR': snr, 'ACCEPTED': psfex_acc,
            'CCD_NB': ccd_nb}

    # Write file
    if verbose:
        print('Writing full PSF catalog file {}...'.format(output_path))
    output.save_as_fits(data, sex_cat_path=starcat_names[0])


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
        cfis.log_command(argv, name='sys.stdout')

    create_full_cat(param.input_dir, param.input_file_pattern, param.output_path, verbose=param.verbose)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
