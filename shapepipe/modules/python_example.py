# -*- coding: utf-8 -*-

"""PYTHON MODULE EXAMPLE

This module defines methods for an example Python module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import time
from numpy.random import randint
from shapepipe.modules.module_decorator import module_runner


class Dummy(object):

    def __init__(self, sleep_time=None):

        if not isinstance(sleep_time, type(None)):
            self.sleep_time = sleep_time

        else:
            self.sleep_time = randint(1, 10)

    def _wait(self):

        time.sleep(self.sleep_time)

    def _read_file(self, file_name):

        with open(file_name) as data_file:
            content = data_file.read().replace('\n', '')

        return content

    def read_files(self, file_name1, file_name2):

        self._wait()
        content1 = self._read_file(file_name1)
        content2 = self._read_file(file_name2)

        self.content = '{} and {}'.format(content1, content2)

    def write_file(self, file_name, message):

        new_content = message + str(self.content)

        text_file = open(file_name, 'w')
        text_file.write(new_content)
        text_file.close()


@module_runner(version='1.0', file_pattern=['numbers', 'letters'],
               file_ext='.txt', depends='numpy', run_method='parallel')
def python_example(input_file_list, run_dirs, file_number_string,
                   config, w_log):

    output_file_name = ('{}/pyex_output{}.cat'.format(run_dirs['output'],
                        file_number_string))
    message = config.get('PYTHON_EXAMPLE', 'MESSAGE')

    inst = Dummy()
    inst.read_files(*input_file_list)
    inst.write_file(output_file_name, message)

    return inst.content, None
