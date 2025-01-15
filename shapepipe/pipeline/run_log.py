"""RUN LOG HANDLING.

This module defines methods for creating a run log.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import numpy as np

from shapepipe.pipeline.shared import split_module_run


class RunLog(object):
    """Run Log Class.

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
        self._runs = get_list(run_log_file)

    def _write(self):
        """Write.

        Write current run to the run log.

        """
        with open(self.run_log_file, 'a') as run_log:
            run_log.write(f'{self.current_run} {self._module_list}\n')

    def get_run(self, search_string):
        """Get Run.

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
            raise RuntimeError(
                f'No runs found matching search string \'{search_string}\'.'
            )

        elif len(runs) > 1:
            raise RuntimeError(
                'More than one run found matching search string '
                + f'\'{search_string}\''
            )

        return runs[0].split(' ')[0]


def get_list(run_log_file):
    """Get List.

    Get a list of all runs in the run log.

    Parameters
    ----------
    run_log_file : str
        Run log file name

    Returns
    -------
    list
        Run log file entries

    """
    with open(run_log_file, 'r') as run_log:
        lines = run_log.readlines()

    runs = [line.rstrip() for line in lines]

    return runs


def get_all(runs, module):
    """Get All.

    Get all previous pipeline runs of a given model.

    Parameters
    ----------
    runs : list
        Log file entries consisting of directory and module(s)
    module : str
        Module name

    Raises
    ------
    RuntimeError
        If no previous runs are found

    Returns
    -------
    list
        All run paths for a given module

    """
    module_base, _ = split_module_run(module)

    all_runs = [
        run for run in runs
        if module_base in run.split()[1].split(',')
    ]
    if len(all_runs) == 0:
        raise RuntimeError(
            f'No previous run of module \'{module_base}\' found'
        )

    all_runs = all_runs[::-1]

    return all_runs


def get_last(runs, module):
    """Get Last.

    Get the last run of the pipeline for a given module.

    Parameters
    ----------
    runs : list
        Log file entries consisting of directory and module(s)
    module : str
        Module name

    Returns
    -------
    str
        The last run for a given module

    """
    all_runs = get_all(runs, module)
    last_run = all_runs[0]

    return last_run.split(' ')[0]


def get_all_dirs(run_log_file, module):
    """Get All Dirs.

    Return directory paths corresponding to all runs of given module.

    Parameters
    ----------
    run_log_file : str
        Run log file name
    module : str
        Module name

    Returns
    -------
    list
        Directory names of all module runs

    """
    runs = get_list(run_log_file)
    all_runs = get_all(runs, module)

    all_dirs = []
    for run in all_runs:
        dir_name = run.split(" ")[0]
        all_dirs.append(f"{dir_name}/{module}/output")
 
    return all_dirs


def get_last_dir(run_log_file, module):
    """Get Last Dir.

    Return directory path corresponding to last run of given module.

    Parameters
    ----------
    run_log_file : str
        Run log file name
    module : str
        Module name

    Returns
    -------
    str
        Directory name of last module run

    """
    all_dirs = get_all_dirs(run_log_file, module)
    last_dir = all_dirs[0]

    return last_dir



