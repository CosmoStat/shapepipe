# -*- coding: utf-8 -*-

"""SEXTRACTOR RUNNER

This module run SExtractor.

:Author: Axel Guinot

"""

import re
from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io

import numpy as np
from sqlitedict import SqliteDict
from astropy.io import fits


def get_zero_point(image_path, key):
    """Get mag zero point

    This function read the magnitude zero-point from the header image.

    Parameters
    ----------
    image_path: str
        Path to the input image
    key: str
        Key to use in the header for the zero-point

    Returns
    -------
    zp: float
        Value associated to the key provided

    """

    h = fits.getheader(image_path)

    zp = h[key]

    try:
        zp = float(zp)
    except:
        raise ValueError('{} is not a float value for zero-point'.format(zp))

    return zp


def make_post_process(cat_path, f_wcs_path, pos_params, ccd_size):
    """Make post process

    This function will add one hdu by epoch to the SExtractor catalog. Only works for tiles.
    The columns will be: NUMBER same as SExtractor NUMBER
                         EXP_NAME name of the single exposure for this epoch
                         CCD_N extansion where the object is

    Parameters
    ----------
    cat_path: str
        Path to the outputed SExtractor catalog
    f_wcs_path: str
        Path to the log file containing wcs for all single exp CCDs
    pos_params: list
        World coordinates to use to match the objects.
    ccd_size: list
        Size of a ccd [nx, ny]

    """

    cat = io.FITSCatalog(cat_path, SEx_catalog=True, open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    cat.open()

    f_wcs = SqliteDict(f_wcs_path)
    n_hdu = len(f_wcs[list(f_wcs.keys())[0]])

    hist = []
    for i in cat.get_data(1)[0][0]:
        if re.split('HISTORY', i)[0] == '':
            hist.append(i)

    exp_list = []
    pattern = r'([0-9]*)\.(.*)'
    for i in hist:
        m = re.search(pattern, i)
        exp_list.append(m.group(1))

    obj_id = np.copy(cat.get_data()['NUMBER'])

    n_epoch = np.zeros(len(obj_id), dtype='int32')
    for i, exp in enumerate(exp_list):
        pos_tmp = np.ones(len(obj_id), dtype='int32') * -1
        for j in range(n_hdu):
            w = f_wcs[exp][j]
            pix_tmp = w.all_world2pix(cat.get_data()[pos_params[0]], cat.get_data()[pos_params[1]], 0)
            ind = ((pix_tmp[0] > int(ccd_size[0])) & (pix_tmp[0] < int(ccd_size[1])) &
                   (pix_tmp[1] > int(ccd_size[2])) & (pix_tmp[1] < int(ccd_size[3])))
            pos_tmp[ind] = j
            n_epoch[ind] += 1
        exp_name = np.array([exp_list[i] for n in range(len(obj_id))])
        a = np.array([(obj_id[ii], exp_name[ii], pos_tmp[ii]) for ii in range(len(exp_name))],
                     dtype=[('NUMBER', obj_id.dtype), ('EXP_NAME', exp_name.dtype), ('CCD_N', pos_tmp.dtype)])
        cat.save_as_fits(data=a, ext_name='EPOCH_{}'.format(i))
        cat.open()

    cat.add_col('N_EPOCH', n_epoch)

    cat.close()


@module_runner(input_module=['split_exp_runner', 'mask_runner_exp'], version='1.0.1',
               file_pattern=['image', 'weight', 'flag'],
               file_ext=['.fits', '.fits', '.fits'],
               executes=['sex'], depends=['numpy', 'sqlitedict', 'astropy'])
def sextractor_runner_exp(input_file_list, run_dirs, file_number_string,
                          config, w_log):

    num = file_number_string

    exec_path = config.getexpanded("SEXTRACTOR_RUNNER_EXP", "EXEC_PATH")
    dot_sex = config.getexpanded("SEXTRACTOR_RUNNER_EXP", "DOT_SEX_FILE")
    dot_param = config.getexpanded("SEXTRACTOR_RUNNER_EXP", "DOT_PARAM_FILE")

    weight_file = config.getboolean("SEXTRACTOR_RUNNER_EXP", "WEIGHT_IMAGE")
    flag_file = config.getboolean("SEXTRACTOR_RUNNER_EXP", "FLAG_IMAGE")
    psf_file = config.getboolean("SEXTRACTOR_RUNNER_EXP", "PSF_FILE")

    zp_from_header = config.getboolean("SEXTRACTOR_RUNNER_EXP", "ZP_FROM_HEADER")
    if zp_from_header:
        zp_key = config.get("SEXTRACTOR_RUNNER_EXP", "ZP_KEY")
        zp_value = get_zero_point(input_file_list[0], zp_key)

    if config.has_option('SEXTRACTOR_RUNNER_EXP', "CHECKIMAGE"):
        check_image = config.getlist("SEXTRACTOR_RUNNER_EXP", "CHECKIMAGE")
    else:
        check_image = ['']

    if config.has_option('SEXTRACTOR_RUNNER_EXP', 'SUFFIX'):
        suffix = config.get('SEXTRACTOR_RUNNER_EXP', 'SUFFIX')
        if (suffix.lower() != 'none') & (suffix != ''):
            suffix = suffix + '_'
        else:
            suffix = ''
    else:
        suffix = ''

    output_file_name = suffix + 'sexcat{0}.fits'.format(num)
    output_file_path = '{0}/{1}'.format(run_dirs['output'], output_file_name)

    command_line = ('{0} {1} -c {2} -PARAMETERS_NAME {3} -CATALOG_NAME {4}'
                    ''.format(exec_path, input_file_list[0], dot_sex,
                              dot_param, output_file_path))

    if zp_from_header:
        command_line += ' -MAG_ZEROPOINT {0}'.format(zp_value)

    extra = 1
    if weight_file:
        command_line += ' -WEIGHT_IMAGE {0}'.format(input_file_list[extra])
        extra += 1
    if flag_file:
        command_line += ' -FLAG_IMAGE {0}'.format(input_file_list[extra])
        extra += 1
    if psf_file:
        command_line += ' -PSF_NAME {0}'.format(input_file_list[extra])
        extra += 1
    if extra != len(input_file_list):
        raise ValueError("Incoherence between input files and keys related "
                         "to extra files: Found {} extra files, but input "
                         "file list lenght is {}"
                         .format(extra, len(input_file_list)))

    if (len(check_image) == 1) & (check_image[0] == ''):
        check_type = ['NONE']
        check_name = ['none']
    else:
        check_type = []
        check_name = []
        for i in check_image:
            check_type.append(i.upper())
            check_name.append(run_dirs['output'] + '/' + suffix+i.lower()+num+'.fits')

    command_line += (' -CHECKIMAGE_TYPE {0} -CHECKIMAGE_NAME {1}'
                     ''.format(','.join(check_type), ','.join(check_name)))

    w_log.info('Calling command \'{}\''.format(command_line))

    stderr, stdout = execute(command_line)

    check_error = re.findall('error', stdout.lower())
    check_error2 = re.findall('all done', stdout.lower())

    if check_error == []:
        stderr2 = ''
    else:
        stderr2 = stdout
    if check_error2 == []:
        stderr2 = stdout

    if config.getboolean("SEXTRACTOR_RUNNER_EXP", "MAKE_POST_PROCESS"):
        f_wcs_path = config.getexpanded("SEXTRACTOR_RUNNER_EXP", "LOG_WCS")
        pos_params = config.getlist("SEXTRACTOR_RUNNER_EXP", "WORLD_POSITION")
        ccd_size = config.getlist("SEXTRACTOR_RUNNER_EXP", "CCD_SIZE")
        make_post_process(output_file_path, f_wcs_path, pos_params, ccd_size)

    return stdout, stderr2
