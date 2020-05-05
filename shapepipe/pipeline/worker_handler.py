# -*- coding: utf-8 -*-

"""WORKER HANDLER

This module defines a class for handling pipeline wokers.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import platform
from os import getpid
from threading import active_count
from modopt.interface.errors import catch_error, warn
from modopt.interface.log import set_up_log, close_log
from shapepipe.pipeline.timeout import with_timeout


class WorkerHandler(object):
    """ Worker Handler

    This class defines the worker to process a given job.

    """

    def __init__(self, verbose=True):

        self.worker_dict = {}
        self._stdout = None
        self._stderr = None
        self._verbose = verbose

    def worker(self, process, job_name, w_log_name, run_dirs, config,
               timeout, module_runner):
        """ Worker

        This method defines a worker.

        Parameters
        ----------
        process : np.ndarray
            File(s) to be processed
        w_log_name : str
            Worker log name
        module_runner : function
            Module runner
        run_dirs : dict
            Run directories
        config : CustomParser
            Configuaration parser instance
        timeout : int
            Timeout limit in seconds

        Returns
        -------
        dict
            Worker dictionary

        """

        self._w_log_name = w_log_name
        self._run_dirs = run_dirs
        self._config = config
        self._module_runner = module_runner
        self._prepare_worker(process, job_name, timeout,
                             module_runner.__name__)
        self._create_worker_log()
        self._run_worker()
        close_log(self.w_log, verbose=False)

        return self.worker_dict

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

        return 'process{}'.format(num)

    def _prepare_worker(self, process, job_name, timeout, module):
        """ Prepare Worker

        This method defines a worker instance dictionary.

        Parameters
        ----------
        Process : str
            File to be processed
        config : CustomParser
            Configuaration parser instance
        timeout : int
            Timeout limit in seconds
        module : str
            Module runner name

        """

        self.worker_dict['pid'] = getpid()
        self.worker_dict['threads'] = active_count()
        self.worker_dict['node'] = platform.node()
        self.worker_dict['system'] = platform.system()
        self.worker_dict['machine'] = platform.machine()
        self.worker_dict['exception'] = False
        self.worker_dict['stderr'] = False
        self.worker_dict['process'] = list(process)
        self.worker_dict['file_number_string'] = job_name
        self.worker_dict['job_name'] = self._set_job_name(job_name)
        self.worker_dict['timeout'] = timeout
        self.worker_dict['module'] = module

    def _create_worker_log(self):
        """ Create Worker Log

        This method prepares a logging instance for the worker and logs the
        worker parameters.

        """

        process_size = len(str(self.worker_dict['process']))

        if self._verbose:

            print(' - {} PID: {} '.format(
                  self.worker_dict['job_name'],
                  self.worker_dict['pid']), end='')

            if (process_size <
                    self._config.getint('WORKER', 'PROCESS_PRINT_LIMIT')):
                print('processing {} {}'.format(
                      self.worker_dict['file_number_string'],
                      self.worker_dict['process']))
            else:
                print()

        self.w_log = set_up_log(self._w_log_name, verbose=False)
        self.worker_dict['log'] = self.w_log.name
        self.w_log.info('Worker process running with:')
        self.w_log.info(' - Job Name: {}'.format(
                        self.worker_dict['job_name']))
        self.w_log.info(' - PID: {}'.format(self.worker_dict['pid']))
        self.w_log.info(' - Threads: {}'.format(self.worker_dict['threads']))
        self.w_log.info(' - Node: {}'.format(self.worker_dict['node']))
        self.w_log.info(' - System: {}'.format(self.worker_dict['system']))
        self.w_log.info(' - Machine: {}'.format(self.worker_dict['machine']))
        self.w_log.info(' - Timeout Limit: {}'.format(
                        self.worker_dict['timeout']))
        self.w_log.info(' - Process: {}'.format(self.worker_dict['process']))

    def _run_worker(self):
        """ Run Worker

        This method runs the worker with a given timeout limit and catches the
        corresponding errors.

        """

        try:
            with_timeout(self.worker_dict['timeout'],
                         self.w_log.name)(self._worker_execution)()

        except Exception as err:
            catch_error(err, self.w_log)
            self.worker_dict['exception'] = type(err).__name__

    def _worker_execution(self):
        """ Worker Execution

        This method executes a worker job and logs the results.

        """

        self._run_module()
        self._log_stdout()

    def _run_module(self):
        """ Run Module

        This method runs a module script.

        Raises
        ------
        RuntimeError
            For non-existent module runner

        """

        self.w_log.info(' - Running module: {}'.format(
                        self.worker_dict['module']))

        file_number_string = self.worker_dict['file_number_string']
        input_file_list = self.worker_dict['process']

        self._stdout, self._stderr = self._module_runner(input_file_list,
                                                         self._run_dirs,
                                                         file_number_string,
                                                         self._config,
                                                         self.w_log)

    def _log_stdout(self):
        """ Log STDOUT

        This method logs the stdout and stderr output of the job.

        """

        self.w_log.info('Process produced the following output: {}'.format(
                        self._stdout))

        if self._stderr:
            self.w_log.info('Process produced the following error(s): {}'
                            ''.format(self._stderr))
            self.worker_dict['stderr'] = True
