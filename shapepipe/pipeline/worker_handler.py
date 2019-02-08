# -*- coding: utf-8 -*-

"""WORKER HANDLER

This module defines a class for handling pipeline wokers.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from os import getpid
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

    def worker(self, job_name, process, filehd, config, timeout, module):
        """ Worker

        This method defines a worker.

        Parameters
        ----------
        job_name : str
            Job name
        Process : str
            File to be processed
        config : CustomParser
            Configuaration parser instance
        timeout : int
            Timeout limit in seconds
        module : str
            Module runner name

        Returns
        -------
        dict
            Worker dictionary

        """

        self._filehd = filehd
        self._config = config
        self._prepare_worker(job_name, process, timeout, module)
        self._create_worker_log()
        self._run_worker()
        close_log(self.w_log, verbose=False)

        return self.worker_dict

    def _prepare_worker(self, job_name, process, timeout, module):
        """ Prepare Worker

        This method defines a worker instance dictionary.

        Parameters
        ----------
        job_name : str
            Job name
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
        self.worker_dict['exception'] = False
        self.worker_dict['stderr'] = False
        self.worker_dict['process'] = process
        self.worker_dict['job_name'] = job_name
        self.worker_dict['timeout'] = timeout
        self.worker_dict['module'] = module

    def _create_worker_log(self):
        """ Create Worker Log

        This method prepares a logging instance for the worker and logs the
        worker parameters.

        """

        if self._verbose:
            print(' - {} PID: {} processing {}'.format(
                  self.worker_dict['job_name'], self.worker_dict['pid'],
                  self.worker_dict['process']))

        self.w_log = set_up_log(self._filehd.get_worker_log_name(
                                self.worker_dict['module'],
                                self.worker_dict['job_name'],
                                self.worker_dict['process'][0]),
                                verbose=False)
        self.worker_dict['log'] = self.w_log.name
        self.w_log.info('Worker process running with:')
        self.w_log.info(' - Job Name: {}'.format(
                        self.worker_dict['job_name']))
        self.w_log.info(' - PID: {}'.format(self.worker_dict['pid']))
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

        moduler_runner = (self._filehd.module_runners[
                          self.worker_dict['module']])

        file_number_string = self.worker_dict['process'][0]
        input_file_list = self.worker_dict['process'][1]
        output_dir = self._filehd.output_dir

        self._stdout, self._stderr = moduler_runner(input_file_list,
                                                    output_dir,
                                                    file_number_string,
                                                    self._config, self.w_log)

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
