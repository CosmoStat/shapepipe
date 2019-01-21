# -*- coding: utf-8 -*-

"""SEXTRACTOR RUNNER

This module run SExtractor.

:Author: Axel Guinot

"""
import re

from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module='mask_runner', version='1.0',
               file_pattern=['image', 'weight', 'flag'], file_ext=['.fits','.fits','.fits'],
               executes='sex')
def sextractor_runner(input_file_list, output_dir, file_number_string,
                   config, w_log):

    num = file_number_string
    output_file_name = 'sexcat{0}.fits'.format(num)
    output_file_path = '{0}/{1}'.format(output_dir, output_file_name)

    exec_path = config.getexpanded("SEXTRACTOR", "EXEC_PATH")
    dot_sex = config.getexpanded("SEXTRACTOR", "DOT_SEX_FILE")
    dot_param = config.getexpanded("SEXTRACTOR", "DOT_PARAM_FILE")

    weight_file = config.getboolean("SEXTRACTOR", "WEIGHT_IMAGE")
    flag_file = config.getboolean("SEXTRACTOR", "FLAG_IMAGE")
    psf_file = config.getboolean("SEXTRACTOR", "PSF_FILE")

    check_image = config.getlist("SEXTRACTOR", "CHECKIMAGE")

    command_line = '{0} {1} -c {2} -PARAMETERS_NAME {3} -CATALOG_NAME {4}'.format(exec_path, input_file_list[0], dot_sex, dot_param, output_file_path)

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
        raise ValueError("Incoherence between input files and keys related to extra files.")

    if (len(check_image) == 1) & (check_image[0] == ''):
        check_type = ['NONE']
        check_name = ['none']
    else:
        check_type = []
        check_name = []
        for i in check_image:
            check_type.append(i.upper())
            check_name.append(output_dir + '/' +i.lower()+num+'.fits')
    
    command_line += ' -CHECKIMAGE_TYPE {0} -CHECKIMAGE_NAME {1}'.format(','.join(check_type), ','.join(check_name))

    stderr, stdout = execute(command_line)

    check_error = re.findall('error', stdout.lower())
    check_error2 = re.findall('all done', stdout.lower())

    if check_error == []:
        stderr2 = ''
    else:
        stderr2 = stdout
    if check_error2 == []:
        stderr2 = stdout

    return stdout, stderr2
