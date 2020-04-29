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
def uncompress_fits_image_runner(input_file_list, output_dir, file_number_string,
                                 config, w_log):

    output_pattern = config.get("UNCOMPRESS_FITS_IMAGE_RUNNER", "OUTPUT_PATTERN")

    try:
        data_hdu = config.getint("UNCOMPRESS_FITS_IMAGE_RUNNER", "HDU_DATA")
    except:
        data_hdu = 0

    data = fits.getdata(input_file_list[0], data_hdu)
    header = fits.getheader(input_file_list[0], data_hdu)

    hdu = fits.PrimaryHDU(data, header)
    hdul = fits.HDUList([hdu])

    hdul.writeto(output_dir + '/' + output_pattern + file_number_string + '.fits')

    return None, None
