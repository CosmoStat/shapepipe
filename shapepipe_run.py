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

    def _get_module_depends(self):
        """ Get Module Dependencies

        List the Python packages and executables needed to run the modules.

        Returns
        -------
        tuple
            List of python dependencies, list of system executables

        """

        depends = []
        executes = []

        module_runners = self._filehd.module_runners

        for module in module_runners.keys():
            depends += module_runners[module].depends
            executes += module_runners[module].executes

        return depends, executes

    def _check_dependencies(self):
        """ Check Dependencies

        Check that all pipeline dependencies have been installed.

        """

        module_dep, module_exe = self._get_module_depends()
        dependencies = __installs__ + module_dep

        dh = DependencyHandler(dependencies, module_exe)

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

    def run(self):
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

        for module in self._modules:

            # Create a job handler for the current module
            jh = JobHandler(module, filehd=self._filehd,
                            config=self._config,
                            log=self.log, verbose=self._verbose)

            # Submit the job handler jobs
            jh.submit_jobs()

            # Update error count
            self._error_count += jh.error_count

        # Finish and close the pipeline log
        self._close_pipeline_log()


def main(args=None):

    try:
        pipe = ShapePipe()
        pipe.run()

    except Exception as err:
        catch_error(err, pipe.log)
        return 1


if __name__ == "__main__":
    main()
