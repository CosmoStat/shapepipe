# -*- coding: utf-8 -*-

"""RUN LOG HANDLING

This module defines methods for creating a run log.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import numpy as np


class RunLog(object):

    def __init__(self, run_log_file, current_run):

        self.run_log_file = run_log_file
        self.current_run = current_run
        self._write()
        self._get_list()

    def _write(self):

        with open(self.run_log_file, 'a') as run_log:
            run_log.write('{}\n'.format(self.current_run))

    def _get_list(self):

        with open(self.run_log_file, 'r') as run_log:
            lines = run_log.readlines()

        self._runs = [line.rstrip() for line in lines]

    def get_last(self):

        return self._runs[self._runs.index(self.current_run) - 1]

    def get_run(self, search_string):

        runs = [run for run in self._runs if search_string in run]

        if len(runs) < 1:
            raise RuntimeError('No runs found matching search string.')

        elif len(runs) > 1:
            raise RuntimeError('More than one run found matching search '
                               'string')

        return runs[0]
