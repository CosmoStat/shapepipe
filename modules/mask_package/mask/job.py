# -*- coding: utf-8 -*-

"""JOB MODULE

This module contain classes for managing the job processes.

:Authors: Samuel Farrens and Marc Gentile

:Date: 31/10/2017 (Happy Halloween!)

"""

# -- Python imports
import os
import sys

# -- External imports
from mpfx import mpfx_job

# --- Module-specific imports
from execute import PackageRunner
from helper import PackageHelper


class PackageJobProcessor(mpfx_job.MpfxJobProcessor):

    """Package Job Processor

    This class contains methods to submit jobs and process associated job
    results.

    Parameters
    ----------
    master : class
        Master class instance (PackageMasterMPI or PackageMasterSMP)

    Notes
    -----
    This class is inherits the methods from MpfxJobProcessor() in mpfx_job.py.

    """

    def __init__(self, master):

        mpfx_job.MpfxJobProcessor.__init__(self, master)

        self._helper = PackageHelper()
        self._master = master

    def create_dataset(self, master, dataset_name, dataset_type,
                       dataset_base_dir, dataset_dir_list, dataset_recurse):

        """Create Dataset

        This method creates the primary dataset.

        Parameters
        ----------
        master : class
            Master process class instance
        dataset_name : str
            Dataset name
        dataset_type : str
            Prefix of a Dataset class
        dataset_base_dir : str
            Dataset base directory
        dataset_dir_list : list, optional
            List of the specific directories to search
        dataset_recurse : bool, optional
            Option to walk down directories (default is 'True')

        Returns
        -------
        class PackageDataSet instance

        Notes
        -----
        This method is called by MpfxJobProcessor.create_jobs() in mpfx_job.py.

        """

        return PackageDataSet(self._master, dataset_name, dataset_base_dir,
                              dataset_dir_list, dataset_recurse)

    def process_job(self, job, worker):

        """Process Job

        This mthod processes a job of class PackageJob and returns the
        corresponding results in the form of a object to the Master.

        Parameters
        ----------
        job : class
            PackageJob class instance
        worker : class
            Worker class instance

        Returns
        -------
        class MpfxJobResult class instance

        Notes
        -----
        This method is called by MasterSMP.process_jobs() in mp_calc_SMP.py.

        """

        results_dict = {}

        try:

            runner = PackageRunner(self._helper, job, worker)

            file_types = runner.file_types

            filepath = [job.get_file_path(file_type) for file_type in
                        file_types]

            if file_types[0] not in results_dict:
                    results_dict[file_types[0]] = {}

            if worker.logging_enabled():
                temp_string = ('{0} - mask: Processing Input Files {1}')
                worker.logger.log_info_p(temp_string.format(worker.name,
                                         filepath))
                worker.logger.flush()

            # -- Run the package executable
            run_dict = runner.run_executable()

            # -- Update the results ditionary
            results_dict[file_types[0]]['run_dict'] = run_dict

            # --- Job result to return
            return mpfx_job.MpfxJobResult(worker, job, results_dict)

        except Exception:

            if worker.logging_enabled():
                temp_string = ('{0} - Some error occurred while processing '
                               'job: {1} ({2})')
                worker.logger.log_error_p(temp_string.format(worker.name, job,
                                          sys.exc_info()[1]))
                worker.logger.flush()

    def process_job_result(self, job_result, master):

        """Process Job Results

        This method updates the output log with the results associated with a
        given job.

        Parameters
        ----------
        job_result : class
            MpfxJobResult class instance
        master : class
            Master class instance

        Notes
        -----
        This method is called by mp_calc_SMP.py MasterSMP.process_jobs()

        """

        # --- Log processing statistics
        if self._master.logging_enabled():
            msg = ('{0} - /{1}/run-{2:03}-{3:1d} - Catalog '
                   'generated')
            msg = msg.format(self._master.name,
                             job_result.job.get_branch_tree(),
                             job_result.job.img_no,
                             job_result.job.epoch)
            self._master.logger.log_info_p(msg)
            self._master.logger.flush()

    def all_jobs_processed(self, *args):

        """All Jobs Processed

        This method is called by the Master once all the jobs have been
        processed.

        Notes
        -----
        This method simply overrides the mpfx_job.py
        MpfxJobProcessor.all_jobs_processed() method.

        """

        pass


class PackageDataSet(mpfx_job.MpfxDataset):

    """Package Data Set

    Notes
    -----
    This class was added to overwrite the is_catalog() method in
    the MpfxDataset class in mpfx_data.py.

    """

    def is_catalog(*args):
        """Is Catalogue

        This method returns true for all is_catalog() calls.

        Returns
        -------
        bool True

        """

        return True
