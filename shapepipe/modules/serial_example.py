# -*- coding: utf-8 -*-

"""SERIAL MODULE EXAMPLE

This module defines methods for an example serial module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner


class Dummy(object):

    def __init__(self):

        pass

    def _read_file(self, file_name):

        with open(file_name) as data_file:
            content = data_file.read().replace('\n', '')

        return content

    def read_files(self, input_file_list):

        self.content = ''

        for file_list in input_file_list:
            for file_name in file_list:
                self.content += self._read_file(file_name)

    def write_file(self, file_name):

        text_file = open(file_name, 'w')
        text_file.write(self.content)
        text_file.close()


@module_runner(input_module='python_example', version='1.0',
               file_pattern=['numbers', 'letters', 'pyex_output'],
               file_ext=['.txt', '.txt', '.cat'], depends='numpy',
               run_method='serial')
def serial_example(input_file_list, run_dirs, file_number_string,
                   config, w_log):

    output_file_name = ('{}/serial_output{}.cat'.format(run_dirs['output'],
                        file_number_string))

    inst = Dummy()
    inst.read_files(input_file_list)
    inst.write_file(output_file_name)

    return inst.content, None
