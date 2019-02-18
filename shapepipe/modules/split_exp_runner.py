# -*- coding: utf-8 -*-

"""SPLIT EXP RUNNER

This module split the different CCD's hdu of a single exposure into separate files.

:Author: Axel Guinot

"""


import numpy as np
import sip_tpv as stp
from astropy.wcs import WCS

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io


def create_hdus(exp_path, output_dir, output_name, output_sufix, n_hdu=40,
                transf_coord=True, transf_int=False, save_header=False):
    """
    """

    exp_file = io.FITSCatalog(exp_path, hdu_no=1)
    exp_file.open()

    header_file = np.zeros(n_hdu, dtype='O')

    for i in range(1, n_hdu+1):

        h = exp_file._cat_data[i].header
        if transf_coord:
            stp.pv_to_sip(h)

        d = exp_file.get_data(i)
        if transf_int:
            d = d.astype(np.int16)

        file_name = (output_dir + '/' + output_sufix + output_name +
                     '-' + str(i-1) + '.fits')
        new_file = io.FITSCatalog(file_name,
                                  open_mode=io.BaseCatalog.OpenMode.ReadWrite)
        new_file.save_as_fits(data=d, image=True, image_header=h)

        w = WCS(h)

        header_file[i-1] = w

    if save_header:
        file_name = output_dir + '/' + 'headers' + output_name + '.npy'
        np.save(file_name, header_file)


@module_runner(version='1.0', file_pattern=['image', 'weight', 'flag'],
               file_ext=['.fz', '.fz', '.fz'],
               depends=['numpy', 'astropy', 'sip_tpv'])
def split_exp_runner(input_file_list, output_dir, file_number_string,
                     config, w_log):

    file_suffix = config.getlist("SPLIT_EXP_RUNNER", "OUTPUT_SUFFIX")
    n_hdu = config.getint("SPLIT_EXP_RUNNER", "N_HDU")

    for exp_path, exp_suffix in zip(input_file_list, file_suffix):

        transf_int = 'flag' in exp_suffix
        transf_coord = 'image' in exp_suffix
        save_header = 'image' in exp_suffix

        create_hdus(exp_path, output_dir, file_number_string, exp_suffix,
                    n_hdu, transf_coord, transf_int, save_header)

    return None, None
