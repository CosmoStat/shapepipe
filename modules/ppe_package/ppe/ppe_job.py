"""!
   @package ppe.ppe_job Job Management
   @author Marc Gentile
   @file ppe_job.py
   Job Management
"""

# -- Python imports
import os
import sys

# -- External imports
from mpfx.mpfx_job import *

# --- Module-specific imports
from ppe_psfex import *
from ppe_help import *


# ----------------------------------------------------------------------------
class PpeJobProcessor(MpfxJobProcessor):
    """!
    Job processor: submit jobs and process associated job results. Based on
    mpf.JobProcessor.

    """

    def __init__(self, master):
        """!
        Job Processor constructor

        """
        MpfxJobProcessor.__init__(self, master)

        self._helper = PpeHelper()

        self._pe_runner = PSFExRunner(self)
        #   self._pe_processor = PSFExProcessor(self)

    # ~~~~~~~~~~
    # Properties
    # ~~~~~~~~~~

    @property
    def pe_runner(self):
        """!
        @return the PSFEx runner instance.

        """

        return self._pe_runner

    #
    # @property
    # def pe_processor(self):
    #    """! @return the PSFEx processor instance. """
    #    return self._pe_processor

    @property
    def helper(self):
        """!
        @return the PpeHelper instance.

        """

        return self._helper

    # ~~~~~~~~~~~~~~~
    # Public methods
    # ~~~~~~~~~~~~~~~

    # ------------------------------------------------------------------------
    def create_dataset(self, master, dataset_name, dataset_type,
                       dataset_base_dir, dataset_dir_list, dataset_recurse):
        """!
        Create a primary Dataset object that represent the data source for
        images and catalogs.
        @param master master process instance
        @param dataset_name dataset name
        @param dataset_type prefix of a Dataset class, assumed to be of the
        form
        @code <dataset_type>Dataset @endcode, like @c MpfxDataset
        @param dataset_base_dir dataset base directory
        @param dataset_dir_list [optional] a list of specific dirs to search
        under base directory
        @param dataset_recurse [optional] tell whether to walk down directories
        (default @c True)
        @return the Dataset instance

        """

        # SF NOTE: this method is called by mpfx_job.py
        # MpfxJobProcessor.create_jobs(). No apparent connection to specific
        # package, could be part of template.

        dataset_class = self.get_dataset_module(master, dataset_name,
                                                dataset_type, dataset_base_dir)

        if not isinstance(dataset_class, type(None)):
            return dataset_class(master, dataset_name, dataset_base_dir,
                                 dataset_dir_list, dataset_recurse)
        else:
            return None

    # ------------------------------------------------------------------------
    def get_dataset_module(self, master, dataset_name, dataset_type,
                           dataset_base_dir):

        """
        Return an instance of PpeDataSet

        """

        # SF NOTE: this method was added to overwrite the defualt mpfx method.

        return PpeDataSet

    # ------------------------------------------------------------------------
    def create_jobs(self, master):
        """!
        Locate all objects to process and create the corresponding jobs.
        @param master Master object instance
        @return the list of created jobs

        """

        # SF NOTE: This method is called by mp_calc.py Master.run() and the
        # JobProcessor is defined in the ..._SMP.py module.
        # No apparent connection to specific package, could be part of
        # template.

        job_list = MpfxJobProcessor.create_jobs(self, master)

        # --- Create output directories for job outputs:
        output_dir_dico = master.config.get_section_data("DIR.OUTPUT")
        for (dir_key, dir_value) in output_dir_dico.items():
            if dir_key.startswith("OUTPUT_"):
                if dir_key.find("LOG") == -1:  # no tree for log dirs
                    base_dir = os.path.join(master.run_output_dir, dir_value)
                    self._create_dir_tree(base_dir, self._job_branches)

        return job_list

    # ------------------------------------------------------------------------
    def create_job(self, master, dataset, *args):
        """!
        Factory method for creating a Job object. Each subclass can define a
        job class derived from mpf.Job and instanciate it here.
        @param master the master process
        @param dataset data source used to create the job
        @param args a list of arguments chosen by the caller and whose nature
        depends on the data necessary to process the job
        @return a new PpeJob job object
        @see Job, PpeJob

        """

        # SF NOTE: This method is called by mpfx_job.py
        # MpfxJobProcessor.create_jobs() and Simply returns an instance of
        # PpeJob(). Maybe package specific.

        return PpeJob(master, dataset, self.helper, *args)

    # ------------------------------------------------------------------------
    def process_job(self, job, worker):
        """!
        Process a job of class PpeJob and return the corresponding results in
        the form of a PpeJobResult object to the Master.

        @param job object of class PpeJob with processed data
        @param worker Worker object instance

        @return an object of class PpeJobResult containing the data of the
        processed job

        @note overrides MpfxJobProcessor.process_job()
        @see Job, MpfxJob, JobResult, MpfxJobResult

        """

        # SF NOTE: This method is called by mp_calc_SMP.py
        # MasterSMP.process_jobs() and is package specific.

        object_per_type_dico = {}

        try:

            for file_type in job.get_file_types():
                filepath = job.get_file_path(file_type)

                if file_type not in object_per_type_dico:
                    object_per_type_dico[file_type] = {}

                if worker.logging_enabled():
                    temp_string = ('{0} - Processing SExtractor Output File '
                                   '{1}...')
                    worker.logger.log_info_p(temp_string.format(worker.name,
                                             filepath))
                    worker.logger.flush()

                # SF NOTE: The following line calls PSFEx
                pe_run_dico = self.pe_runner.run_PSFEx(file_type, job, worker)
                pe_output_cat_filepath = pe_run_dico['pe_output_cat_filepath']
                object_per_type_dico[file_type]['pe_run_dico'] = pe_run_dico

                initial_results = PpeJobResult(object_per_type_dico, job,
                                               worker)

                # SF NOTE: NEED TO FIX THIS FOR OUTPUT FILE TESTING
                fake_dict = {'initial_object_count': 1,
                             'final_object_count': 1, 'elapsed_time': 1}
                object_per_type_dico[file_type]['pe_xform_dico'] = fake_dict

                #   --- Process the generated PSFEx catalog
                #   object_per_type_dico[file_type]["pe_xform_dico"] =\
                #                        self.pe_processor.process_catalog(
                #                                       pe_output_cat_filepath,
                #                                    pe_output_check_filepath,
                #                                       file_type, job, worker)

                # --- Job result to return
                job_result = PpeJobResult(object_per_type_dico, job, worker)

        except Exception:

            if worker.logging_enabled():
                temp_string = ('{0} - Some error occurred while processing '
                               'job: {1} ({2})')
                worker.logger.log_error_p(temp_string.format(worker.name, job,
                                          sys.exc_info()[1]))
                worker.logger.flush()

        # --- Create a PpeJobResult object with the results from all catalogues
        print "PSFEx JOB PROCESSED"
        return PpeJobResult(object_per_type_dico, job, worker)

    # ------------------------------------------------------------------------
    def process_job_result(self, job_result, master):
        """!
        Process the result associated with a processed job.

        @param job_result object of class PpeJobResult with processed data
        @param master Master object instance
        @note overrides MpfxJobProcessor.process_job_result()
        @see Job, JobResult, MpfxJobResult
        """

        # SF NOTE: This method is called by mp_calc_SMP.py
        # MasterSMP.process_jobs() and is package specific

        # --- Iterate through the results of all file types in the job...
        object_per_type_dico = job_result.result
        for file_type in object_per_type_dico.keys():

            pe_xform_dico = object_per_type_dico[file_type]['pe_xform_dico']

            if len(pe_xform_dico) > 0:

                # --- Log processing statistics
                if master.logging_enabled():
                    msg = ('{0} - /{1}/cat-{2:03}-{3:1d} - {4} - Catalog '
                           'generated - Object count: {5} -> {6} - Processor '
                           'time: {7:.2f} sec')
                    msg = msg.format(master.name,
                                     job_result.job.get_branch_tree(),
                                     job_result.job.img_no,
                                     job_result.job.epoch, file_type,
                                     pe_xform_dico["initial_object_count"],
                                     pe_xform_dico["final_object_count"],
                                     pe_xform_dico["elapsed_time"])
                    master.logger.log_info_p(msg)
                    master.logger.flush()

    # ------------------------------------------------------------------------
    def all_jobs_processed(self, master):
        """!
        This method is called by the Master once all the jobs have been
        processed.
        @param master instance of the Master
        @note Job results can be selectively obtained by calling the
        get_job_results() methods, specifying the relevant query criteria
        @note the entire list of Job results can aldo be obtained with by
        calling JobProcessor.job_result_list()
        @see Job, JobResult, get_job_result()
        @note PpeJobResult objects can also be directly queried from the
        ppe.ppe_job.PpeJobProcessor.job_result_dico dictionary

        """

        # SF NOTE: This method simply overrides the mpfx_job.py
        # MpfxJobProcessor.all_jobs_processed() method. Could certainly be
        # part of a template

        pass

    # ~~~~~~~~~~~~~~~
    # Private methods
    # ~~~~~~~~~~~~~~~

    # ------------------------------------------------------------------------
    def _create_dir_tree(self, base_dir, branches):
        """!
        Create the directory tree for storing the results.
        @param base_dir base directory
        @param branches list of ordered directory names, making a path tree
        @note example of directory layout: control/ground/constant
        @see configuration file
        """

        # SF NOTE: this method does not appear to package specific, could be
        # part of a template

        for branch in branches:
            dir_names = branch.split('/')
            branch_path = base_dir
            for dir_name in dir_names:
                branch_path = os.path.join(branch_path, dir_name)
                self.helper.make_dir(branch_path)


# ---------------------------------------------------------------------------
class PpeJob(MpfxJob):
    """!
    Represents a Ppe job. Based on mpf.Job and mpfcs82.MpfxJob.

    """

    def __init__(self, master, dataset, helper, img_path_dico, *args):
        """!
        Construct a new job for identiofiying a set of files in the dataset to
        process together.

        To obtain the path associated with a given type, just use:

        @code self.img_path_dico[type]@endcode with one of the aforementioned
        types
        @param master instance of the Master
        @param dataset data source used to create the job
        @param helper helper object for this class
        @param img_path_dico file absolute path dictionary  for every file
        types referenced in the job
        @param args a list of extra parameters

        @return an initialized PpeJob object

        """

        # SF NOTE: Adds methods to the MpfxJob() class.
        # Package specific

        # --- Construct a new Job
        MpfxJob.__init__(self, master, dataset, img_path_dico, *args)
        self._helper = helper

    # ~~~~~~~~~~
    # Properties
    # ~~~~~~~~~~

    @property
    def helper(self):
        """! @return the class' helper object

        """

        return self._helper


# ----------------------------------------------------------------------------
class PpeJobResult(MpfxJobResult):
    """!
    Represents the result of a GREAT3 job. See mpfcs82.MpfxJobResult and
    mpf.JobResult.

    """

    # SF NOTE: Not sure if this class is package specific

    def __init__(self, result, job, worker):
        """!
        Construct a PpeJobResult job result object

        """

        MpfxJobResult.__init__(self, worker, job, result)


# ----------------------------------------------------------------------------
class PpeDataSet(MpfxDataset):

    # SF NOTE: This class was added to overwrite the is_catalog() method in
    # the mpfx_data.py MpfxDataset class.

    def __init__(self, master, dataset_name, dataset_base_dir,
                 dataset_dir_list, dataset_recurse_dirs):
        """!
        Construct a PpeDataSet object

        """

        MpfxDataset.__init__(self, master, dataset_name, dataset_base_dir,
                             dataset_dir_list, dataset_recurse_dirs)

    def is_catalog(*args):
        """
        Return true for is_catalog() calls.

        """

        return True

# -- EOF ppe_job.py
