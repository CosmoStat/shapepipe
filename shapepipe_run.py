#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""SHAPEPIPE SCRIPT

This script runs the shape measurement pipeline.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from datetime import datetime
from modopt.interface.errors import catch_error
from modopt.interface.log import set_up_log, close_log
from shapepipe.info import shapepipe_logo, line, __installs__
from shapepipe.pipeline.args import create_arg_parser
from shapepipe.pipeline.config import create_config_parser
from shapepipe.pipeline.dependency_handler import DependencyHandler
from shapepipe.pipeline.file_handler import FileHandler
from shapepipe.pipeline.job_handler import JobHandler

from mpi4py import MPI
from shapepipe.pipeline.mpi_run import mpi_run


class ShapePipe():
    """ ShapePipe

    ShapePipe runner class.

    """

    def __init__(self):

        self._args = create_arg_parser()
        self._config = create_config_parser(self._args.config)
        self._set_run_name()
        self._modules = self._config.getlist('EXECUTION', 'MODULE')
        self._filehd = FileHandler(self._run_name, self._modules, self._config)
        self._verbose = self._config.getboolean('DEFAULT', 'VERBOSE')
        self._error_count = 0
        self._prep_run()

    def _set_run_name(self):
        """ Set Run Name

        Set the name of the current pipeline run.

        """

        self._run_name = self._config.get('DEFAULT', 'RUN_NAME')

        if self._config.getboolean('DEFAULT', 'RUN_DATETIME'):
            self._run_name += datetime.now().strftime('_%Y-%m-%d_%H-%M-%S')

    def _create_pipeline_log(self):
        """ Create Pipeline Log

        Create a general logging instance for the pipeline run.

        """

        self.log = set_up_log(self._filehd.log_name, verbose=False)

        start_text = 'Starting ShapePipe Run: {}'.format(self._run_name)

        self.log.info(shapepipe_logo())
        self.log.info(start_text)
        self.log.info('')

        if self._verbose:
            print(shapepipe_logo())
            print(start_text)
            print('')

    def _close_pipeline_log(self):
        """ Close Pipeline Log

        Close general logging instance for the pipeline run.

        """

        final_error_count = ('A total of {} errors were recorded.'.format(
                             self._error_count))
        end_text = 'Finishing ShapePipe Run'

        self.log.info(final_error_count)
        self.log.info(end_text)
        self.log.info(line())
        close_log(self.log, verbose=False)

        if self._verbose:
            print(final_error_count)
            print(end_text)
            print(line())

    def _get_module_depends(self, property):
        """ Get Module Dependencies

        List the Python packages and executables needed to run the modules.

        Returns
        -------
        tuple
            List of python dependencies, list of system executables

        """

        prop_list = []

        module_runners = self._filehd.module_runners

        for module in module_runners.keys():

            if self._config.has_option(module.upper(), property.upper()):
                prop_list += self._config.getlist(module.upper(),
                                                  property.upper())
            else:
                prop_list += getattr(module_runners[module], property)

            if self._filehd.get_add_module_property(module, property):
                prop_list += self._filehd.get_add_module_property(module,
                                                                  property)

        return prop_list

    def _check_dependencies(self):
        """ Check Dependencies

        Check that all pipeline dependencies have been installed.

        """

        module_dep = self._get_module_depends('depends') + __installs__
        module_exe = self._get_module_depends('executes')

        dh = DependencyHandler(module_dep, module_exe)

        dep_text = 'Checking Python Dependencies:'
        exe_text = 'Checking System Executables:'

        self.log.info(dep_text)
        if self._verbose:
            print(dep_text)

        for dep in dh.check_dependencies():

            self.log.info(dep)

            if self._verbose:
                print(dep)

        self.log.info('')
        if self._verbose:
            print('')

        self.log.info(exe_text)
        if self._verbose:
            print(exe_text)

        for exe in dh.check_executables():

            self.log.info(exe)

            if self._verbose:
                print(exe)

        self.log.info('')
        if self._verbose:
            print('')

    def _check_module_versions(self):
        """ Check Module Version

        Check versions of the modules.

        """

        ver_text = 'Checking Module Versions:'

        self.log.info(ver_text)
        if self._verbose:
            print(ver_text)

        for module in self._modules:

            module_txt = (' - {} {}'.format(
                          module,
                          self._filehd.module_runners[module].version))

            self.log.info(module_txt)
            if self._verbose:
                print(module_txt)

        self.log.info('')
        if self._verbose:
            print('')

    def _prep_run(self):
        """ Run

        Run the pipeline.

        """

        # Make output directories for the pipeline run
        self._filehd.create_global_run_dirs()

        # Make a log for the pipeline run
        self._create_pipeline_log()

        # Check the pipeline dependencies
        self._check_dependencies()

        # Check the versions of these modules
        self._check_module_versions()

        # for module in self._modules:
        #
        #     # Create a job handler for the current module
        #     jh = JobHandler(module, filehd=self._filehd,
        #                     config=self._config,
        #                     log=self.log, verbose=self._verbose)
        #
        #     # Submit the job handler jobs
        #     jh.submit_jobs()
        #
        #     # Update error count
        #     self._error_count += jh.error_count
        #
        # # Finish and close the pipeline log
        # self._close_pipeline_log()


def split(container, count):

    return [container[_i::count] for _i in range(count)]


def print_thing(thing, rank):

    print('{}: {}'.format(rank, thing))


def run(pipe, comm):

    if comm.rank == 0:
        modules = pipe._modules

    else:
        modules = None

    modules = comm.bcast(modules, root=0)

    for module in modules:

        if comm.rank == 0:
            filehd = pipe._filehd
            config = pipe._config
            verbose = pipe._verbose
            jh = JobHandler(module, filehd=filehd, config=config,
                            log=pipe.log, verbose=verbose)
            timeout = jh.timeout
            job_names = jh._job_names
            process_list = list(jh.filehd.process_list.items())
            jobs = split(list(zip(job_names, process_list)), comm.size)
        else:
            jh = None
            timeout = None
            jobs = None
            filehd = None
            config = None
            verbose = None

        filehd = comm.bcast(filehd, root=0)
        config = comm.bcast(config, root=0)
        verbose = comm.bcast(verbose, root=0)
        timeout = comm.bcast(timeout, root=0)
        jobs = comm.scatter(jobs, root=0)

        results = []
        for job in jobs:
            res = mpi_run(job, filehd, config, timeout, module, verbose)
            results.append(res)

        results = comm.gather(results, root=0)

        if comm.rank == 0:
            jh._worker_dicts = [_i for temp in results for _i in temp]
            jh._check_for_errors()
            jh._check_missed_processes()
            jh.log.info('All processes complete')
            jh.log.info('')

            if jh._verbose:
                print('All processes complete')
                print('')

            pipe._error_count += jh.error_count

    if comm.rank == 0:
        pipe._close_pipeline_log()

    # if comm.rank == 0:
    #
    #     for module in pipe._modules[:1]:
    #
    #         # Create a job handler for the current module
    #         jh = JobHandler(module, filehd=pipe._filehd, config=pipe._config,
    #                         log=pipe.log, verbose=pipe._verbose)
    #
    # filehd = comm.bcast(pipe._filehd, root=0)
    # config = comm.bcast(pipe._config, root=0)
    # verbose = comm.bcast(pipe._verbose, root=0)

        #     exit()
        #
        #     # Submit the job handler jobs
        #     jh.submit_jobs()
        #
        #     # Update error count
        #     pipe._error_count += jh.error_count
        #
        # pipe._close_pipeline_log()


def main(args=None):

    comm = MPI.COMM_WORLD

    try:

        if comm.rank == 0:
            pipe = ShapePipe()
        else:
            pipe = None

        run(pipe, comm)

    except Exception as err:
        if comm.rank == 0:
            catch_error(err, pipe.log)
            return 1


if __name__ == "__main__":
    main()
