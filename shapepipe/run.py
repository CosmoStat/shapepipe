"""SHAPEPIPE RUN.

This module sets up a given run of the shape measurement pipeline.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import sys
from datetime import datetime

from joblib import cpu_count
from modopt.interface.errors import catch_error
from modopt.interface.log import close_log, set_up_log

from shapepipe.info import __installs__, line, shapepipe_logo
from shapepipe.pipeline.args import create_arg_parser
from shapepipe.pipeline.config import create_config_parser
from shapepipe.pipeline.dependency_handler import DependencyHandler
from shapepipe.pipeline.file_handler import FileHandler
from shapepipe.pipeline.job_handler import JobHandler
from shapepipe.pipeline.mpi_run import split_mpi_jobs, submit_mpi_jobs

try:
    from mpi4py import MPI
except ImportError:  # pragma: no cover
    import_mpi = False
else:
    import_mpi = True


class ShapePipe():
    """ShapePipe.

    ShapePipe runner class.

    """

    def __init__(self):

        self.log = None

    def set_up(self):
        """Set Up.

        Set up ShapePipe properties.

        """
        self._args = create_arg_parser()
        self.config = create_config_parser(self._args.config)
        self._set_run_name()
        self.modules = self.config.getlist('EXECUTION', 'MODULE')
        self.mode = self.config.get('EXECUTION', 'MODE').lower()
        self.exclusive=self._args.exclusive
        self.verbose = self.config.getboolean('DEFAULT', 'VERBOSE')
        self.filehd = FileHandler(
            self._run_name,
            self.modules,
            self.config,
            exclusive=self._args.exclusive,
            verbose=self.verbose,
        )
        self.error_count = 0
        self._prep_run()

    def _set_run_name(self):
        """Set Run Name.

        Set the name of the current pipeline run.

        """
        self._run_name = self.config.get('DEFAULT', 'RUN_NAME')

        if self.config.getboolean('DEFAULT', 'RUN_DATETIME'):
            self._run_name += datetime.now().strftime('_%Y-%m-%d_%H-%M-%S')

    def _create_pipeline_log(self):
        """Create Pipeline Log.

        Create a general logging instance for the pipeline run.

        """
        self.log = set_up_log(self.filehd.log_name, verbose=False)

        start_text = f'Starting ShapePipe Run: {self._run_name}'

        self.log.info(shapepipe_logo())
        self.log.info(start_text)
        self.log.info('')

        if self.verbose:
            print(shapepipe_logo(colour=True))
            print(start_text)
            print('')

        # Temporary fix to give file handler access to the log. This should
        # be improved at some point.
        self.filehd.log = self.log

    def close_pipeline_log(self):
        """Close Pipeline Log.

        Close general logging instance for the pipeline run.

        Raises
        ------
        RunTimeError
            if error occurs during pipeline run

        """
        if self.error_count == 1:
            plur = ' was'
        else:
            plur = 's were'
        final_error_count = (
            f'A total of {self.error_count} error{plur} recorded.'
        )
        end_text = 'Finishing ShapePipe Run'

        self.log.info(final_error_count)
        self.log.info(end_text)
        self.log.info(line())
        close_log(self.log, verbose=False)

        if self.verbose:
            print(final_error_count)
            print(end_text)
            print(line())

        if self.error_count > 0:
            raise RuntimeError(final_error_count)

    def _get_module_depends(self, property):
        """Get Module Dependencies.

        List the Python packages and executables needed to run the modules.

        Parameters
        ----------
        property : str
            Module property to be checked

        Returns
        -------
        tuple
            List of python dependencies, list of system executables

        """
        prop_list = []

        module_runners = self.filehd.module_runners

        for module in module_runners.keys():

            if self.config.has_option(module.upper(), property.upper()):
                prop_list += self.config.getlist(
                    module.upper(),
                    property.upper(),
                )
            else:
                prop_list += getattr(module_runners[module], property)

            if self.filehd.get_add_module_property(module, property):
                prop_list += self.filehd.get_add_module_property(
                    module,
                    property,
                )

        return prop_list

    def _check_dependencies(self):
        """Check Dependencies.

        Check that all pipeline dependencies have been installed.

        """
        module_dep = self._get_module_depends('depends') + __installs__
        module_exe = self._get_module_depends('executes')

        module_dep += ['mpi4py'] if import_mpi else module_dep

        dh = DependencyHandler(module_dep, module_exe)

        dep_text = 'Checking Python Dependencies:'
        exe_text = 'Checking System Executables:'

        self.log.info(dep_text)
        if self.verbose:
            print(dep_text)

        for dep in dh.check_dependencies():

            self.log.info(dep)

            if self.verbose:
                print(dep)

        self.log.info('')
        if self.verbose:
            print('')

        self.log.info(exe_text)
        if self.verbose:
            print(exe_text)

        for exe in dh.check_executables():

            self.log.info(exe)

            if self.verbose:
                print(exe)

        self.log.info('')
        if self.verbose:
            print('')

    def _check_module_versions(self):
        """Check Module Version.

        Check versions of the modules.

        """
        ver_text = 'Checking Module Versions:'

        self.log.info(ver_text)
        if self.verbose:
            print(ver_text)

        for module in set(self.modules):

            module_txt = (
                f' - {module} {self.filehd.module_runners[module].version}'
            )

            self.log.info(module_txt)
            if self.verbose:
                print(module_txt)

        self.log.info('')
        if self.verbose:
            print('')

    def _check_system_setup(self):
        """Check System Set Up.

        Check the set up of the machine on which the pipeline is running.

        """
        setup_text = 'Checking System Set Up:'
        cpu_info = f' - Number of available CPUs: {cpu_count()}'

        self.log.info(setup_text)
        self.log.info(cpu_info)
        self.log.info('')

        if self.verbose:
            print(setup_text)
            print(cpu_info)
            print('')

    def _get_module_run_methods(self):
        """Get Module Run Method.

        Create a dictionary of modules with corresponding run methods.

        """
        self.run_method = {}

        for module in self.modules:

            self.run_method[module] = (
                self.filehd.module_runners[module].run_method
            )

    def _prep_run(self):
        """Prepare Run.

        Prepare to run the pipeline.

        """
        # Make output directories for the pipeline run
        self.filehd.create_global_run_dirs()

        # Make a log for the pipeline run
        self._create_pipeline_log()

        # Check the pipeline dependencies
        self._check_dependencies()

        # Check the versions of the modules
        self._check_module_versions()

        # Check the system set up
        self._check_system_setup()

        # Get run method for each module
        self._get_module_run_methods()

    def record_mode(self):
        """Record Mode.

        Log mode in which ShapePipe is running.

        """
        mode_text = f'Running ShapePipe using {self.mode}'

        self.log.info(mode_text)
        self.log.info('')
        if self.verbose:
            print(mode_text)
            print('')


def run_smp(pipe):
    """Run SMP.

    Run ShapePipe using SMP.

    Parameters
    ----------
    pipe : ShapePipe
        ShapePipe instance

    """
    # Loop through modules to be run
    for module in pipe.modules:

        # Create a job handler for the current module
        jh = JobHandler(
            module,
            filehd=pipe.filehd,
            config=pipe.config,
            log=pipe.log,
            job_type=pipe.run_method[module],
            exclusive=pipe.exclusive,
            verbose=pipe.verbose,
        )

        # Submit jobs
        jh.submit_jobs()

        # Update error count
        pipe.error_count += jh.error_count

        # Delete job handler
        del jh

    # Finish and close the pipeline log
    pipe.close_pipeline_log()


def run_mpi(pipe, comm):
    """Run MPI.

    Run ShapePipe using MPI.

    Parameters
    ----------
    pipe : ShapePipe
        ShapePipe instance
    comm : MPI.COMM_WORLD
        MPI common world instance

    """
    # Assign master node
    master = comm.rank == 0

    # Get the module to be run
    modules = pipe.modules if master else None
    modules = comm.bcast(modules, root=0)

    # Get ShapePipe objects
    if master:
        config = pipe.config
        verbose = pipe.verbose
    else:
        config = verbose = None
    config = comm.bcast(config, root=0)
    verbose = comm.bcast(verbose, root=0)

    # Loop through modules to be run
    for module in modules:

        # Run set up on master
        if master:
            # Create a job handler for the current module
            jh = JobHandler(
                module,
                filehd=pipe.filehd,
                config=config,
                log=pipe.log,
                job_type=pipe.run_method[module],
                parallel_mode='mpi',
                exclusive=pipe.exclusive,
                verbose=verbose,
            )

            # Get job type
            job_type = jh.job_type

            # Handle serial jobs
            if job_type == 'serial':
                jh.submit_jobs()

            # Handle parallel jobs
            else:
                # Get JobHandler objects
                timeout = jh.timeout
                # Get file handler objects
                run_dirs = jh.filehd.module_run_dirs
                module_runner = jh.filehd.module_runners[module]
                worker_log = jh.filehd.get_worker_log_name
                # Define process list
                process_list = jh.filehd.process_list
                # Define job list
                jobs = split_mpi_jobs(process_list, comm.size)
                del process_list
        else:
            job_type = module_runner = worker_log = timeout = \
                jobs = run_dirs = None

        # Broadcast job type to all nodes
        job_type = comm.bcast(job_type, root=0)

        if job_type == 'parallel':

            # Broadcast objects to all nodes
            run_dirs = comm.bcast(run_dirs, root=0)

            module_runner = comm.bcast(module_runner, root=0)
            worker_log = comm.bcast(worker_log, root=0)
            timeout = comm.bcast(timeout, root=0)
            jobs = comm.scatter(jobs, root=0)

            # Submit the MPI jobs and gather results
            results = comm.gather(
                submit_mpi_jobs(
                    jobs,
                    config,
                    timeout,
                    run_dirs,
                    module_runner,
                    worker_log,
                    verbose
                ),
                root=0,
            )

            # Delete broadcast objects
            del module_runner, worker_log, timeout, jobs

            # Finish up parallel jobs
            if master:
                # Assign worker dictionaries
                jh.worker_dicts = jh.filehd.flatten_list(results)
                # Finish up job handler session
                jh.finish_up()
                # Delete results
                del results

        if master:
            # Update error count
            pipe.error_count += jh.error_count
            # Delete job handler
            del jh

    # Finish and close the pipeline log
    pipe.close_pipeline_log() if master else None


def run(*args):
    """Run ShapePipe.

    This function runs ShapePipe.

    """
    try:

        if import_mpi:
            comm = MPI.COMM_WORLD
            master = comm.rank == 0
        else:
            master = True

        if master:
            pipe = ShapePipe()
            pipe.set_up()
            mode = pipe.mode
        else:
            pipe = None
            mode = None

        mode = comm.bcast(mode, root=0) if import_mpi else 'smp'

        if master:
            pipe.mode = mode
            pipe.record_mode()

        if mode == 'mpi':
            run_mpi(pipe, comm)
        else:
            run_smp(pipe)

    except Exception as err:
        if master:
            catch_error(err, log=pipe.log)
            return 1
