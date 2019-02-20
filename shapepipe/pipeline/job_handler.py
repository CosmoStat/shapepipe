# -*- coding: utf-8 -*-

"""JOB HANDLER

This module defines a class for handling pipeline jobs.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from modopt.interface.errors import warn
from joblib import Parallel, delayed, cpu_count
from configparser import ConfigParser
from logging import Logger
from shapepipe.pipeline.worker_handler import WorkerHandler

import psutil
import gc


def memory_usage_psutil():
    gc.collect()
    return psutil.Process().memory_info().rss / (2 ** 20)


class JobHandler(object):
    """ Job Handler

    This class handles the submition of jobs to workers distributed among a
    specified number of CPUs.

    Parameters
    ----------
    module : str
        Module name
    filehd : FileHandler
        File handler instance
    config : CustomParser
        Configuaration parser instance
    log : logging.Logger
        Logging instance
    batch_size : int, optional
        Number of jobs to submitted simultaneously, the default is 1
    timeout : int, optional
        Timeout limit for a given job in seconds, the default is None
    verbose : bool, optional
        Verbose setting, default is True

    """

    def __init__(self, module, filehd, config, log, batch_size=1,
                 backend='loky', timeout=None, verbose=True):

        self.filehd = filehd
        self.log = log
        self.config = config
        self._module = module
        self._module_runner = self.filehd.module_runners[self._module]
        self.error_count = 0
        self._verbose = verbose

        # Set the bacth size
        if self.config.has_option('JOB', 'SMP_BATCH_SIZE'):
            self.batch_size = self.config.getint('JOB', 'SMP_BATCH_SIZE')
        else:
            self.batch_size = batch_size

        # Set the backend
        if self.config.has_option('JOB', 'SMP_BACKEND'):
            self.backend = self.config.get('JOB', 'SMP_BACKEND')
        else:
            self.backend = backend

        # Set the job timeout limit
        if self.config.has_option('JOB', 'TIMEOUT'):

            time_str = self.config.get('JOB', 'TIMEOUT')

            if ':' in time_str:
                self.timeout = self.hms2sec(time_str)
            else:
                self.timeout = int(time_str)

        else:
            self.timeout = timeout

        # Set up module in file handler
        self.filehd.set_up_module(self._module)

        # Set the total number of jobs
        self._n_jobs = len(self.filehd.process_list)

        # Set the job names
        self.job_names = [self._set_job_name(job) for job in
                          range(self._n_jobs)]

        self._log_job_parameters()

    @property
    def config(self):
        """ Config

        This method defines the configuation parser instance

        Raises
        ------
        TypeError
            For incorrect input type

        """

        return self._config

    @config.setter
    def config(self, value):

        if not isinstance(value, ConfigParser):
            raise TypeError('config must be an instane of '
                            'configparser.ConfigParser')

        self._config = value

    @property
    def log(self):
        """ Log

        This method defines the logging instance

        Raises
        ------
        TypeError
            For incorrect input type

        """

        return self._log

    @log.setter
    def log(self, value):

        if not isinstance(value, Logger):
            raise TypeError('log must be an instance of logging.Logger.')

        self._log = value

    @property
    def batch_size(self):
        """ Batch Size

        This method defines the job batch size.

        Raises
        ------
        ValueError
            For invalid batch size value

        """

        return self._batch_size

    @batch_size.setter
    def batch_size(self, value):

        if not isinstance(value, int) or (value < 1):
            raise ValueError('Batch size must be an integer >= 1.')

        if value > cpu_count():
            warn('Batch size exeeds the number of available CPUs.')

        self._batch_size = value

    @property
    def timeout(self):
        """ Timeout Limit

        This method defines the timeout limit for all jobs.

        Raises
        ------
        TypeError
            For incorrect input type
        ValueError
            For invalid timeout limit value

        """

        return self._timeout

    @timeout.setter
    def timeout(self, value):

        if not isinstance(value, (type(None), int)):
            raise TypeError('Timeout must be None or an integer.')

        if isinstance(value, int) and (value < 1):
            raise ValueError('Timeout limit must be >= 1.')

        self._timeout = value

    def finish_up(self):
        """ Finish Up

        Finish up JobHandler session.

        """

        self._check_for_errors()
        self._check_missed_processes()
        self.log.info('All processes complete')
        self.log.info('')

        if self._verbose:
            print('All processes complete')
            print('Memory Usage:', memory_usage_psutil())
            print('')

    def submit_smp_jobs(self):
        """ Submit Jobs

        This method submits the jobs and checks the results.

        """

        self._distribute_smp_jobs()
        self.finish_up()

    @staticmethod
    def hms2sec(time_str):
        """ HMS to Seconds

        Convert a string from hours, minutes and seconds to seconds.

        Parameters
        ----------
        time_str :  str
            Time string

        Returns
        -------
        int
            Time in seconds

        Notes
        -----
        Time strings should take the form 'HH:MM:SS'.

        """

        h, m, s = time_str.split(':')

        return int(h) * 3600 + int(m) * 60 + int(s)

    @staticmethod
    def _set_job_name(num):
        """ Set Job Name

        This method creates a job name for a given process number.

        Parameters
        ----------
        num : int
            Process number

        Returns
        -------
        str
            Job name

        """

        return 'process_{}'.format(num)

    def _log_job_parameters(self):
        """ Log Job Parameters

        This method logs the class instance parameters.

        """

        text = 'Starting job handler with:'
        cpu_info = ' - Number of available CPUs: {}'.format(cpu_count())
        job_info = ' - Total number of jobs: {}'.format(self._n_jobs)
        batch_info = ' - Batch size: {}'.format(self.batch_size)
        time_info = ' - Timeout Limit: {}s'.format(self.timeout)
        module_info = ' - Module: {}'.format(self._module)

        if self._verbose:
            print('Starting job handler with:')
            print(cpu_info)
            print(job_info)
            print(batch_info)
            print(time_info)
            print(module_info)

        # Log process properties
        self.log.info(text)
        self.log.info(cpu_info)
        self.log.info(job_info)
        self.log.info(batch_info)
        self.log.info(time_info)
        self.log.info(module_info)

    def _distribute_smp_jobs(self):
        """ Distribute SMP Jobs

        This method distributes the jobs to the workers using SMP.

        """

        result = (Parallel(n_jobs=self.batch_size, backend=self.backend)
                  (delayed(WorkerHandler(verbose=self._verbose).worker)
                   (job_name, process,
                    self.filehd.get_worker_log_name(self._module, job_name,
                                                    process[0]),
                    self.filehd.output_dir, self.config, self.timeout,
                    self._module_runner)
                   for job_name, process in
                   zip(self.job_names, self.filehd.process_list.items())))

        self.worker_dicts = result

    def _check_for_errors(self):
        """ Check for Errors

        This method checks the worker dictionaries for errors and exceptions.

        """

        # Check worker dictionaries for errors
        self._check_exception_status()
        self._check_stderr_status()

    def _check_exception_status(self):
        """ Check Exception Status

        This method checks the worker dictionaries for exceptions raised by
        Python and logs the instances.

        """

        for worker_dict in self.worker_dicts:
            if worker_dict['exception']:
                self.log.info('ERROR: {} recorded in: {}'.format(
                              worker_dict['exception'], worker_dict['log']))
                self.error_count += 1

    def _check_stderr_status(self):
        """ Check STDERR Status

        This method checks the worker dictionaries for errors raised by
        stderr and logs the instances.

        """

        for worker_dict in self.worker_dicts:
            if worker_dict['stderr']:
                self.log.info('ERROR: stderr recorded in: {}'.format(
                              worker_dict['log']))
                self.error_count += 1

    def _check_missed_processes(self):
        """ Check Missed Processes

        This method checks the file handler for processes that were not
        submitted.

        """

        missed_txt = (' - The following processes were not submitted to '
                      'workers:')

        if self.filehd.missed:

            self.log.info(missed_txt)
            self.log.info(' - {}'.format(self.filehd.missed))

            if self._verbose:
                print(missed_txt)
                print(' - {}'.format(self.filehd.missed))
