# -*- coding: utf-8 -*-

"""VIGNET MAKER RUNNER

This file contains methods to create postage stamps from images.

:Author: Axel Guinot

"""

import numpy as np
from astropy.wcs import WCS
import re
from sf_tools.image.stamp import FetchStamps
import shapepipe.pipeline.file_io as io
from shapepipe.modules.module_decorator import module_runner


def get_pos(galcat_path, pos_params, pos_type):
    """Get positions

    Get the positions of the given parameters from SExtractor catalog.

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalog
    pos_params : list
        List of string containing the SExtractor's parameters to use as
        positions

    Returns
    -------
    numpy.ndarray
        Array of the positions

    """

    f = io.FITSCatalog(galcat_path, SEx_catalog=True)
    f.open()

    pos = np.array([f.get_data()[pos_params[1]],
                    f.get_data()[pos_params[0]]]).T

    f.close()

    return pos


def convert_pos(pos, image_path):
    """Convert position

    Convert positions from world coordinates to pixel coordinates.

    Parameters
    ----------
    pos : numpy.ndarray
        Positions in world coordinates.
    image_path : str
        Path to the image from where the stamp are created.

    Return
    ------
    numpy.ndarray
        New positions in pixel coordinates.

    """

    f = io.FITSCatalog(image_path)
    f.open()
    h = f.get_header(0)
    f.close()

    w = WCS(h)

    pos_tmp = np.copy(pos)
    pos_tmp[:, [0, 1]] = pos_tmp[:, [1, 0]]

    new_pos = w.all_world2pix(pos_tmp, 1)

    new_pos[:, [0, 1]] = new_pos[:, [1, 0]]

    return new_pos


def get_stamp(img_path, galcat_path, pos, rad):
    """Get stamp

    Extract stamp at given positions on the image

    Parameters
    ----------
    img_path : str
        Path to the image
    galcat_path : str
        Path to the SExtractor catalog
    pos_params : list
        List of string containing the SExtractor's parameters to use as
        positions
    rad : int
        Radius of the stamp, must be odd

    Returns
    -------
    numpy.array
        Array containing the vignets

    """

    img_file = io.FITSCatalog(img_path)
    img_file.open()
    img = img_file.get_data(0)
    img_file.close()

    fs = FetchStamps(img, int(rad))
    fs.get_pixels(np.round(pos).astype(int))

    vign = fs.scan()

    return vign


def get_original_vignet(galcat_path):
    """Get original vignet

    Get the vignets from the SExtractor catalog

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalog

    Returns
    -------
    numpy.ndarray
        Array containing the vignets

    """

    f = io.FITSCatalog(galcat_path, SEx_catalog=True)
    f.open()

    vign = f.get_data()['VIGNET']

    f.close()

    return vign


def make_mask(galcat_path, mask_value):
    """Make mask

    Change the value of the SExtractor mask on vignet

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalog
    mask_value : float
        New value of the mask

    Returns
    -------
    numpy.array
        Array of the vignets with the new mask value

    """

    vign = get_original_vignet(galcat_path)

    vign[np.where(vign == -1e30)] = mask_value

    return vign


def save_vignet(vign, sexcat_path, output_dir, suffix):
    """Save vignet

    Save the vignet into a SExtractor format catalog

    Parameters
    ----------
    sexcat_path : str
        Path to the original SExtractor catalog
    output_dir : str
        Path to the output directory

    """

    s = re.split(r"\-([0-9]*)\-([0-9]+)\.", sexcat_path)
    num = '-{0}-{1}'.format(s[1], s[2])

    output_name = output_dir + '/' + suffix + '_vignet{}.fits'.format(num)
    f = io.FITSCatalog(output_name, SEx_catalog=True,
                       open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    f.save_as_fits(vign, names=['VIGNET'], sex_cat_path=sexcat_path)


@module_runner(input_module='setools_runner',
               file_pattern=['galaxy_selection', 'image'],
               file_ext=['.fits', '.fits'], depends=['numpy', 'astropy', 'sf_tools'])
def vignetmaker_runner(input_file_list, output_dir, file_number_string,
                       config, w_log):

    galcat_path = input_file_list[0]

    do_masking = config.getboolean("VIGNETMAKER_RUNNER", "MASKING")
    if do_masking:
        mask_value = config.getfloat("VIGNETMAKER_RUNNER", "MASK_VALUE")
        vign = make_mask(galcat_path, mask_value)
        save_vignet(vign, galcat_path, output_dir, 'cat')

    else:
        stamp_size = config.getint("VIGNETMAKER_RUNNER", "STAMP_SIZE") - 1
        if stamp_size % 2 != 0:
            raise ValueError("The STAMP_SIZE must be odd")
        rad = int(stamp_size/2)

        suffix = config.getlist("VIGNETMAKER_RUNNER", "SUFFIX")
        if len(suffix) != len(input_file_list[1:]):
            raise ValueError("You must provide a suffix for each image from "
                             "which you extract stamps.")

        pos_type = config.get("VIGNETMAKER_RUNNER", "COORD")
        pos_params = config.getlist("VIGNETMAKER_RUNNER", "POSITION_PARAMS")

        pos = get_pos(galcat_path, pos_params, pos_type)

        for n, img in enumerate(input_file_list[1:]):
            image_path = img

            if pos_type == 'PIX':
                pass
            elif pos_type == 'SPHE':
                pos = convert_pos(pos, image_path)
            else:
                raise ValueError('Coordinates type must be in : PIX (pixel), '
                                 'SPHE (spherical).')

            vign = get_stamp(image_path, galcat_path, pos-1, rad)

            save_vignet(vign, galcat_path, output_dir, suffix[n])

    return None, None
