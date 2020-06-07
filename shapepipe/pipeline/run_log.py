# -*- coding: utf-8 -*-

"""RUN LOG HANDLING

This module defines methods for creating a run log.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import numpy as np


class RunLog(object):
    """ Run Log Class

    This class manages the run log for ShapePipe.

    Parameters
    ----------
    run_log_file : str
        Name of the run log file
    module_list : list
        List of modules to be run
    current_run : str
        Name of the current run

    """

    def __init__(self, run_log_file, module_list, current_run):

        self.run_log_file = run_log_file
        self._module_list = ','.join(module_list)
        self.current_run = current_run
        self._write()
        self._get_list()

    def _write(self):
        """ Write

        Write current run to the run log.

        """

        with open(self.run_log_file, 'a') as run_log:
            run_log.write('{} {}\n'.format(self.current_run,
                          self._module_list))

    def _get_list(self):
        """ Get List

        Get a list of all runs in the run log.

        """

        with open(self.run_log_file, 'r') as run_log:
            lines = run_log.readlines()

        self._runs = [line.rstrip() for line in lines]

    def get_last(self, module):
        """ Get Last

        Get the last run of the pipeline for a given module.

        Parameters
        ----------
        module : str
            Module name

        Returns
        -------
        str
            The last run for a given module

        """

        all_runs = [run for run in self._runs if module in
                    run.split()[1].split(',')]
        if len(all_runs) == 0:
            raise RuntimeError('Last-run module \'{}\' not '
                               'found'.format(module))
        last_run = all_runs[::-1][0]

        return last_run.split(' ')[0]

    def get_run(self, search_string):
        """ Get Run

        Get a specific run that matches the input search string.

        Parameters
        ----------
        search_string : str
            Pattern to match to runs

        Returns
        -------
        str
            The run matching the search string

        Raises
        ------
        RuntimeError
            If no runs found that match the search string
        RuntimeError
            If more than one run is found matches the search string

        """

        runs = [run for run in self._runs if search_string in run]

        if len(runs) < 1:
            raise RuntimeError('No runs found matching search string.')

        elif len(runs) > 1:
            raise RuntimeError('More than one run found matching search '
                               'string')

        return runs[0].split(' ')[0]
