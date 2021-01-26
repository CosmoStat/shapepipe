# -*- coding: utf-8 -*-

"""PSFEX RUNNER

This module run PSFEx.

:Author: Axel Guinot

"""

import re
import os
from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module='setools_runner', version='1.0',
               file_pattern=['star_selection'], file_ext=['.fits'],
               executes='psfex')
def psfex_runner(input_file_list, run_dirs, file_number_string,
                 config, w_log):

    exec_path = config.getexpanded("PSFEX_RUNNER", "EXEC_PATH")
    dot_sex = config.getexpanded("PSFEX_RUNNER", "DOT_PSFEX_FILE")
    outcat_name = '{0}/psfex_cat{1}.cat'.format(run_dirs['output'],
                                                file_number_string)

    command_line = ('{0} {1} -c {2} -PSF_DIR {3} -OUTCAT_NAME {4}'
                    ''.format(exec_path, input_file_list[0], dot_sex,
                              run_dirs['output'], outcat_name))

    if config.has_option('PSFEX_RUNNER', "CHECKIMAGE"):
        check_image = config.getlist("PSFEX_RUNNER", "CHECKIMAGE")
    else:
        check_image = ['']
    if (len(check_image) == 1) & (check_image[0] == ''):
        check_type = ['NONE']
        check_name = ['none']
    else:
        suffix = re.split(file_number_string, os.path.splitext(
                          os.path.split(input_file_list[0])[-1])[0])[0]
        check_type = []
        check_name = []
        for i in check_image:
            check_type.append(i.upper())
            check_name.append(run_dirs['output'] + '/' + suffix + '_' +
                              i.lower() + file_number_string+'.fits')

    command_line += (' -CHECKIMAGE_TYPE {0} -CHECKIMAGE_NAME {1}'
                     ''.format(','.join(check_type), ','.join(check_name)))

    w_log.info('Running command \'{}\''
               ''.format(command_line))

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
