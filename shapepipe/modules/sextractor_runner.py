# -*- coding: utf-8 -*-

"""SEXTRACTOR RUNNER

This module run SExtractor.

:Author: Axel Guinot

"""

import re
from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner

from shapepipe.modules.SExtractor_runner import sextractor_script as ss


@module_runner(
    input_module='mask_runner',
    version='1.0.1',
    file_pattern=['image', 'weight', 'flag'],
    file_ext=['.fits', '.fits', '.fits'],
    executes=['sex'],
    depends=['numpy'],
)
def sextractor_runner(input_file_list, run_dirs, file_number_string,
                      config, module_config_sec, w_log):

    num = file_number_string

    exec_path = config.getexpanded(module_config_sec, "EXEC_PATH")
    dot_sex = config.getexpanded(module_config_sec, "DOT_SEX_FILE")
    dot_param = config.getexpanded(module_config_sec, "DOT_PARAM_FILE")
    dot_conv = config.getexpanded(module_config_sec, "DOT_CONV_FILE")

    weight_file = config.getboolean(module_config_sec, "WEIGHT_IMAGE")
    flag_file = config.getboolean(module_config_sec, "FLAG_IMAGE")
    psf_file = config.getboolean(module_config_sec, "PSF_FILE")
    detection_image = config.getboolean(module_config_sec, "DETECTION_IMAGE")
    detection_weight = config.getboolean(
        module_config_sec, "DETECTION_WEIGHT"
    )

    zp_from_header = config.getboolean(module_config_sec, "ZP_FROM_HEADER")
    if zp_from_header:
        zp_key = config.get(module_config_sec, "ZP_KEY")
    else:
        zp_key = None

    bkg_from_header = config.getboolean(module_config_sec, "BKG_FROM_HEADER")
    if bkg_from_header:
        bkg_key = config.get(module_config_sec, "BKG_KEY")
    else:
        bkg_key = None

    if config.has_option(module_config_sec, "CHECKIMAGE"):
        check_image = config.getlist(module_config_sec, "CHECKIMAGE")
    else:
        check_image = ['']

    if config.has_option(module_config_sec, 'SUFFIX'):
        suffix = config.get(module_config_sec, 'SUFFIX')
    else:
        suffix = None

    SE_caller = ss.sextractor_caller(
        input_file_list,
        run_dirs['output'],
        num,
        dot_sex,
        dot_param,
        dot_conv,
        weight_file,
        flag_file,
        psf_file,
        detection_image,
        detection_weight,
        zp_from_header,
        bkg_from_header,
        zero_point_key=zp_key,
        background_key=bkg_key,
        check_image=check_image,
        output_suffix=suffix,
    )

    command_line = SE_caller.make_command_line(exec_path)
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

    if config.getboolean(module_config_sec, "MAKE_POST_PROCESS"):
        f_wcs_path = config.getexpanded(module_config_sec, "LOG_WCS")
        pos_params = config.getlist(module_config_sec, "WORLD_POSITION")
        ccd_size = config.getlist(module_config_sec, "CCD_SIZE")
        ss.make_post_process(
            SE_caller.path_output_file,
            f_wcs_path,
            pos_params, ccd_size
        )

    return stdout, stderr2
