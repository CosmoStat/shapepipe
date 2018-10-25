# -*- coding: utf-8 -*-

"""PYTHON MODULE EXAMPLE

This module defines methods for an example Python module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import time
import numpy as np
from .module_decorator import module_runner


class Dummy(object):

    def __init__(self, sleep_time=None):

        if not isinstance(sleep_time, type(None)):
            self.sleep_time = sleep_time

        else:
            self.sleep_time = np.random.randint(1, 10)

    def _wait(self):

        time.sleep(self.sleep_time)

    def read_file(self, file_name):

        self._wait()
        self.content = np.genfromtxt(file_name)

    def write_file(self, file_name, message):

        new_content = message + str(self.content)

        text_file = open(file_name, 'w')
        text_file.write(new_content)
        text_file.close()


@module_runner(n_inputs=1, ext='.txt')
def python_example(worker_dict, filehd, config, w_log):

    input_file_name = worker_dict['process']
    output_dir = filehd.output_dir
    output_file_name = ('{}/{}.cat'.format(output_dir,
                        worker_dict['job_name']))
    message = config.get('PYTHON_EXAMPLE', 'MESSAGE')

    inst = Dummy()
    inst.read_file(input_file_name)
    inst.write_file(output_file_name, message)

    return inst.content, None
