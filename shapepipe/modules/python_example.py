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
               file_ext='.txt', depends='numpy')
def python_example(worker_dict, filehd, config, w_log):

    input_file_names = worker_dict['process']
    output_dir = filehd.output_dir
    output_file_name = ('{}/{}.cat'.format(output_dir,
                        worker_dict['job_name']))
    message = config.get('PYTHON_EXAMPLE', 'MESSAGE')

    inst = Dummy()
    inst.read_files(*input_file_names)
    inst.write_file(output_file_name, message)

    return inst.content, None
