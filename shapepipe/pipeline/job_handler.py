# -*- coding: utf-8 -*-

"""JOB HANDLER

This module defines a class for handling pipeline jobs.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from configparser import ConfigParser
from gc import collect
from joblib import Parallel, delayed, cpu_count
from logging import Logger
from modopt.interface.errors import warn
from shapepipe.pipeline.worker_handler import WorkerHandler


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
    job_type : str, optional
        Job type, the default is 'parallel'
    parallel_mode : str, optional
        Parallisation mode, default is 'smp'
    batch_size : int, optional
        Number of jobs to submitted simultaneously, the default is None
    backend : str, optional
        Joblib backend, the default is None (which corresponds to 'loky')
    timeout : int, optional
        Timeout limit for a given job in seconds, the default is None
    verbose : bool, optional
        Verbose setting, default is True

    """

    def __init__(self, module, filehd, config, log, job_type='parallel',
                 parallel_mode='smp', batch_size=None, backend=None,
                 timeout=None, verbose=True):

        self.filehd = filehd
        self.log = log
        self.job_type = job_type
        self.parallel_mode = parallel_mode
        self.config = config
        self.batch_size = batch_size
        self.backend = backend
        self.timeout = timeout
        self._module = module
        self._module_runner = self.filehd.module_runners[self._module]
        self.error_count = 0
        self._verbose = verbose

        # Set up module in file handler
        self.filehd.set_up_module(self._module)

        # Set the total number of jobs
        self._n_procs = len(self.filehd.process_list)

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
    def job_type(self):
        """ Job Type

        This method defines the job type

        Raises
        ------
        TypeError
            For incorrect input type

        """

        return self._job_type

    @job_type.setter
    def job_type(self, value):

        if value not in ('serial', 'parallel'):
            raise TypeError('{} is not a valid job type.'.format(value))

        self._job_type = value

    @property
    def parallel_mode(self):
        """ Parallel Mode

        This method defines the mode of parallelisation.

        Raises
        ------
        TypeError
            For incorrect input type

        """

        return self._parallel_mode

    @parallel_mode.setter
    def parallel_mode(self, value):

        if value not in ('smp', 'mpi'):
            raise TypeError('{} is not a valid parallel mode.'.format(value))

        self._parallel_mode = value

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

        if (isinstance(value, type(None)) and
                self.config.has_option('JOB', 'SMP_BATCH_SIZE')):
            value = self.config.getint('JOB', 'SMP_BATCH_SIZE')

        elif isinstance(value, type(None)):
            value = 1

        if not isinstance(value, int) or (value < 1):
            raise ValueError('Batch size must be an integer >= 1.')

        if value > cpu_count():
            warn('Batch size exeeds the number of available CPUs.')

        self._batch_size = value

    @property
    def backend(self):
        """ Backend

        This method defines the joblib backend. The default is 'loky'.

        Raises
        ------
        ValueError
            For invalid backend value

        """

        return self._backend

    @backend.setter
    def backend(self, value):

        if (isinstance(value, type(None)) and
                self.config.has_option('JOB', 'SMP_BACKEND')):
            value = self.config.get('JOB', 'SMP_BACKEND').lower()
        elif isinstance(value, type(None)):
            value = 'loky'

        if value not in ('loky', 'multiprocessing', 'threading'):
            raise ValueError('{} is not a valid joblib backend.'.format(value))

        self._backend = value

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

        if (isinstance(value, type(None)) and
                self.config.has_option('JOB', 'TIMEOUT')):
            value = self.config.get('JOB', 'TIMEOUT')
            value = self.hms2sec(value) if ':' in value else int(value)

        if not isinstance(value, (type(None), int)):
            raise TypeError('Timeout must be None or an integer.')

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
            print('')

        collect()

        self.clean_up()

    def submit_jobs(self):
        """ Submit Jobs

        Submit jobs in serial or parallel.

        """

        if self.job_type == 'serial':
            self.submit_serial_job()
        else:
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

    def _log_job_parameters(self):
        """ Log Job Parameters

        This method logs the class instance parameters.

        """

        text = 'Starting job handler with:'
        module_info = ' - Module: {}'.format(self._module)
        cpu_info = ' - Number of available CPUs: {}'.format(cpu_count())
        proc_info = ' - Total number of processes: {}'.format(self._n_procs)
        job_type = ' - Job Type: {}'.format(self.job_type)
        batch_info = ' - Batch size: {}'.format(self.batch_size)
        time_info = ' - Timeout Limit: {}s'.format(self.timeout)

        if self._verbose:
            print('Starting job handler with:')
            print(module_info)
            print(cpu_info)
            print(proc_info)
            print(job_type)
            if self.job_type == 'parallel' and self.parallel_mode == 'smp':
                print(batch_info)
            print(time_info)

        # Log process properties
        self.log.info(text)
        self.log.info(module_info)
        self.log.info(cpu_info)
        self.log.info(proc_info)
        self.log.info(job_type)
        if self.job_type == 'parallel' and self.parallel_mode == 'smp':
            self.log.info(batch_info)
        self.log.info(time_info)

    def _distribute_smp_jobs(self):
        """ Distribute SMP Jobs

        This method distributes the jobs to the workers using SMP.

        """

        result = (Parallel(n_jobs=self.batch_size, backend=self.backend)
                  (delayed(WorkerHandler(verbose=self._verbose).worker)
                   (process[1:], process[0],
                    self.filehd.get_worker_log_name(self._module,
                                                    process[0]),
                    self.filehd.module_run_dirs, self.config, self.timeout,
                    self._module_runner)
                   for process in self.filehd.process_list))

        self.worker_dicts = result

    def submit_serial_job(self):
        """ Submit Serial Job

        Submit a single serial job with access to all processes.

        """

        wh = WorkerHandler(verbose=self._verbose)
        process = self.filehd.process_list

        result = wh.worker(process, '',
                           self.filehd.get_worker_log_name(self._module,
                                                           '_serial'),
                           self.filehd.module_run_dirs, self.config,
                           self.timeout, self._module_runner)

        self.worker_dicts = [result]

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

    def clean_up(self):
        """ Finish

        Finish job handler instance.

        """

        self.filehd.remove_process_mmap()
