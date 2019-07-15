# -*- coding: utf-8 -*-

"""SWARP RUNNER

This module run swarp.

:Author: Axel Guinot

"""

import re
import os

import numpy as np
from astropy.wcs import WCS
import sip_tpv as stp

from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io


def get_history(coadd_path, image_path_list):
    """ Get history

    Write in the coadd header the single exposures used and how many CCDs from
    them.

    Parameters
    ----------
    coadd_path : str
        Path to the coadd image.
    image_path_list : list
        List of the single exposures path to check

    """
    coadd_file = io.FITSCatalog(coadd_path, hdu_no=0,
                                open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    coadd_file.open()
    wcs_coadd = WCS(coadd_file.get_header())
    corner_coadd = wcs_coadd.calc_footprint().T

    for img_path in image_path_list:
        if (img_path == '\n') or (img_path == ''):
            continue
        img_path = img_path.replace('\n', '')
        img_path = img_path.replace(' ', '')

        f_tmp = io.FITSCatalog(img_path)
        f_tmp.open()
        n_hdu = len(f_tmp.get_ext_name())
        ccd_inter = 0
        for ext in range(1, n_hdu):
            h_tmp = f_tmp._cat_data[ext].header
            stp.pv_to_sip(h_tmp)
            wcs_tmp = WCS(h_tmp)
            corner_tmp = wcs_tmp.calc_footprint().T
            if (np.min(corner_coadd[0]) > np.max(corner_tmp[0]) or
                    np.max(corner_coadd[0]) < np.min(corner_tmp[0])):
                continue
            if (np.min(corner_coadd[1]) > np.max(corner_tmp[1]) or
                    np.max(corner_coadd[1]) < np.min(corner_tmp[1])):
                continue

            ccd_inter += 1
        f_tmp.close()

        if ccd_inter != 0:
            coadd_file.add_header_card("HISTORY",
                                       "From file {} {} extension(s) used"
                                       "".format(os.path.split(img_path)[1],
                                                 ccd_inter))
    coadd_file.close()


@module_runner(version='1.0',
               file_pattern=['tile'],
               file_ext=['.txt'],
               executes=['swarp'], depends=['numpy', 'astropy', 'sip_tpv'])
def swarp_runner(input_file_list, run_dirs, file_number_string,
                 config, w_log):

    num = '-' + re.split('-', file_number_string)[1]

    exec_path = config.getexpanded("SWARP_RUNNER", "EXEC_PATH")
    dot_swarp = config.getexpanded("SWARP_RUNNER", "DOT_SWARP_FILE")
    image_prefix = config.get("SWARP_RUNNER", "IMAGE_PREFIX")
    weight_prefix = config.get("SWARP_RUNNER", "WEIGHT_PREFIX")

    if config.has_option('SWARP_RUNNER', 'SUFFIX'):
        suffix = config.get('SWARP_RUNNER', 'SUFFIX')
        if (suffix.lower() != 'none') and (suffix != ''):
            suffix = suffix + '_'
        else:
            suffix = ''
    else:
        suffix = ''

    output_image_name = suffix + 'image{0}.fits'.format(num)
    output_weight_name = suffix + 'weight{0}.fits'.format(num)
    output_image_path = '{0}/{1}'.format(run_dirs['output'], output_image_name)
    output_weight_path = '{0}/{1}'.format(run_dirs['output'],
                                          output_weight_name)

    # Get center position
    tmp = os.path.split(os.path.splitext(input_file_list[0])[0])[1]
    tmp = re.split('_|-', tmp)
    ra, dec = tmp[1], tmp[2]

    # Get weight list
    image_file = open(input_file_list[0])
    image_list = image_file.readlines()
    image_file.close()
    weight_list = []
    for img_path in image_list:
        tmp = os.path.split(img_path)
        new_name = tmp[1].replace(image_prefix,
                                  weight_prefix).replace('\n', '')
        weight_list.append('/'.join([tmp[0], new_name]))

    command_line = '{} @{} -c {}' \
                   ' -WEIGHT_IMAGE {}' \
                   ' -IMAGEOUT_NAME {} -WEIGHTOUT_NAME {}' \
                   ' -RESAMPLE_SUFFIX .resamp{}.fits ' \
                   ' -CENTER_TYPE MANUAL -CENTER {},{} ' \
                   ''.format(exec_path, input_file_list[0], dot_swarp,
                             ','.join(weight_list), output_image_path,
                             output_weight_path, num, ra, dec)

    stderr, stdout = execute(command_line)

    check_error = re.findall('error', stdout.lower())
    check_error2 = re.findall('all done', stdout.lower())

    if check_error == []:
        stderr2 = ''
    else:
        stderr2 = stdout
    if check_error2 == []:
        stderr2 = stdout

    get_history(output_image_path, image_list)

    return stdout, stderr2
