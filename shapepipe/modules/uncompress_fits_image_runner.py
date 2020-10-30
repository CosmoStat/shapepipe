# -*- coding: utf-8 -*-

"""UNCOMPRESS_FITS_IMAGE_RUNNER

This module uncompress fits images and save them on a single hdu fits.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner

from astropy.io import fits


@module_runner(version='1.0',
               file_pattern=['image'],
               file_ext=['.fits'],
               numbering_scheme='_0')
def uncompress_fits_image_runner(input_file_list, run_dirs, file_number_string,
                                 config, w_log):

    output_pattern_list = config.getlist('UNCOMPRESS_FITS_IMAGE_RUNNER', 'OUTPUT_PATTERN')

    if config.has_option('UNCOMPRESS_FITS_IMAGE_RUNNER', 'HDU_DATA'):
        data_hdu = config.getint("UNCOMPRESS_FITS_IMAGE_RUNNER", "HDU_DATA")
    else:
        data_hdu = 0

    if len(input_file_list) != len(output_pattern_list):
        raise ValueError('Lists INPUT_PATH ({}) and OUTPUT_PATTERN ({}) '
                         'need to be of equal length.'
                         ''.format(len(input_file_list),
                                   len(output_pattern_list)))

    # Read data from input list files
    for i in range(len(input_file_list)):
        data = fits.getdata(input_file_list[i], data_hdu)
        header = fits.getheader(input_file_list[i], data_hdu)

        # Create and write new FITS file with that HDU only
        hdu = fits.PrimaryHDU(data, header)
        hdul = fits.HDUList([hdu])
        hdul.writeto('{}/{}{}.fits'.format(run_dirs['output'], output_pattern_list[i], file_number_string))

    return None, None
