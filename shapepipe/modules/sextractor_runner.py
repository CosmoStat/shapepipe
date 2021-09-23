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

    exec_path = config.getexpanded("SEXTRACTOR_RUNNER", "EXEC_PATH")
    dot_sex = config.getexpanded("SEXTRACTOR_RUNNER", "DOT_SEX_FILE")
    dot_param = config.getexpanded("SEXTRACTOR_RUNNER", "DOT_PARAM_FILE")
    dot_conv = config.getexpanded("SEXTRACTOR_RUNNER", "DOT_CONV_FILE")

    weight_file = config.getboolean("SEXTRACTOR_RUNNER", "WEIGHT_IMAGE")
    flag_file = config.getboolean("SEXTRACTOR_RUNNER", "FLAG_IMAGE")
    psf_file = config.getboolean("SEXTRACTOR_RUNNER", "PSF_FILE")
    detection_image = config.getboolean("SEXTRACTOR_RUNNER", "DETECTION_IMAGE")
    detection_weight = config.getboolean(
        "SEXTRACTOR_RUNNER", "DETECTION_WEIGHT"
    )

    zp_from_header = config.getboolean("SEXTRACTOR_RUNNER", "ZP_FROM_HEADER")
    if zp_from_header:
        zp_key = config.get("SEXTRACTOR_RUNNER", "ZP_KEY")
    else:
        zp_key = None

    bkg_from_header = config.getboolean("SEXTRACTOR_RUNNER", "BKG_FROM_HEADER")
    if bkg_from_header:
        bkg_key = config.get("SEXTRACTOR_RUNNER", "BKG_KEY")
    else:
        bkg_key = None

    if config.has_option('SEXTRACTOR_RUNNER', "CHECKIMAGE"):
        check_image = config.getlist("SEXTRACTOR_RUNNER", "CHECKIMAGE")
    else:
        check_image = ['']

    if config.has_option('SEXTRACTOR_RUNNER', 'SUFFIX'):
        suffix = config.get('SEXTRACTOR_RUNNER', 'SUFFIX')
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

    if config.getboolean("SEXTRACTOR_RUNNER", "MAKE_POST_PROCESS"):
        f_wcs_path = config.getexpanded("SEXTRACTOR_RUNNER", "LOG_WCS")
        pos_params = config.getlist("SEXTRACTOR_RUNNER", "WORLD_POSITION")
        ccd_size = config.getlist("SEXTRACTOR_RUNNER", "CCD_SIZE")
        ss.make_post_process(
            SE_caller.path_output_file,
            f_wcs_path,
            pos_params, ccd_size
        )

    return stdout, stderr2
