"""!
   ppe_psfex.py - PSFEx catalog generation and processing
"""

# -- Python imports
import os
import time
import shutil


# -------------------------------------------------------------------------------------------------
class PSFExRunner(object):

    """!
      Run PSFEx and create a PSFEx catalog from an input SExtractor catalogue
    """

    # ------------------------------------------------------------------------
    def __init__(self, job_processor):

        self._job_processor = job_processor    # the JobProcessor object
        self._helper = job_processor.helper   # utility class for convenience

    # ~~~~~~~~~~~
    # Properties
    # ~~~~~~~~~~~

    # --- Getters

    @property
    def helper(self):
        """! @return the Helper class
        """

        return self._helper

    @property
    def job_processor(self):
        """! @return the JobProcessor object.
        """

        return self._job_processor

    # ~~~~~~~~~~~~~~
    # Public methods
    # ~~~~~~~~~~~~~~

    # -----------------------------------------------------------------------
    def run_PSFEx(self, file_type, job, worker):

        # SF NOTE: This is the method that actually calls PSFEx from the
        # system. It is called by PpeJobProcessor.process_job().

        pe_object_dico = {}     # Result dictionary to return populated

        try:

            # --- Keep track of execution time
            start_time = time.clock()

            # --- SExtractor Catalogue path to analyse
            sexcat_filepath = job.get_file_path(file_type)

            # --- Get the relevant configuration information
            config = worker.config

            # --- Input directory
            input_dir = os.path.abspath(worker.base_input_dir)

            # PSFEx configuration
            pe_config_filepath = (self._get_psfex_config_filepath(
                                  sexcat_filepath, file_type, job, worker))

            if not isinstance(pe_config_filepath, type(None)):
                if worker.logging_enabled():
                    temp_string = ('{0} - /{1}/cat-{2:03}-{3:1d} - {4} - '
                                   'Using PSFEx .psfex configuration file: '
                                   '{5} ...')
                    worker.logger.log_info_p(temp_string.format(worker.name,
                                             job.get_branch_tree(), job.img_no,
                                             job.epoch, file_type,
                                             pe_config_filepath))
                    worker.logger.flush()

            # --- Make a copy of the used .psfex file to the log directory for
            # record
            shutil.copy(pe_config_filepath, worker.log_output_dir)

            # --- Target directory where to store files
            pe_output_path = os.path.join(worker.result_output_dir,
                                          job.get_branch_tree())

            # --- PSFEx catalog output file name
            pe_output_cat_filename = (self._get_output_pe_catalog_filename(
                                      sexcat_filepath, job, config))

            # --- PSFEx catalog output file path
            pe_output_cat_filepath = os.path.abspath(os.path.join(
                                                     pe_output_path,
                                                     pe_output_cat_filename))

            if worker.logging_enabled():
                temp_string = ('{0} - /{1}/cat-{2:03}-{3:1d} - {4} - '
                               'Generating PSFEx catalog {5} ...')
                worker.logger.log_info_p(temp_string.format(worker.name,
                                         job.get_branch_tree(), job.img_no,
                                         job.epoch, file_type,
                                         pe_output_cat_filepath))
                worker.logger.flush()

            # --- Verbose type
            pe_verbose_type = config.get_as_string('PE_VERBOSE_TYPE', 'PSFEX')

            # --- PSFEx Execution line
            pe_exec_path = config.get_as_string('PE_EXEC_PATH', 'PSFEX')
            pe_exec_line = ('{0} {1} -PSF_DIR {2} -VERBOSE_TYPE {3}').format(
                            pe_exec_path, sexcat_filepath, pe_output_path,
                            pe_verbose_type)

            # --- Execute PSFEx
            cur_dir = os.getcwd()
            os.chdir(input_dir)
            os.system(pe_exec_line)
            os.chdir(cur_dir)

            # --- Dictionary with the relevant information for later processing
            pe_object_dico = {}
            pe_object_dico['pe_output_cat_filepath'] = pe_output_cat_filepath
            pe_object_dico['elapsed_time'] = time.clock() - start_time

            if worker.logging_enabled():
                if os.path.exists(pe_output_cat_filepath):
                    temp_string = ('{0} - /{1}/cat-{2:03}-{3:1d} - {4} - '
                                   'Catalog {5} generated successfully')
                    worker.logger.log_info_p(temp_string.format(worker.name,
                                             job.get_branch_tree(), job.img_no,
                                             job.epoch, file_type,
                                             pe_output_cat_filepath))
                else:
                    temp_string = ('{0} - An error occurred while generating '
                                   'PSFEx catalog: {1}')
                    worker.logger.log_error_p(temp_string.format(worker.name,
                                              job))
                worker.logger.flush()

        except Exception as detail:

            if worker.logging_enabled():
                temp_string = ('{0} - An error occurred while generating '
                               'PSFEx catalog: {1} ({2})')
                worker.logger.log_error_p(temp_string.format(worker.name, job,
                                          detail))
                worker.logger.flush()

        return pe_object_dico

    # ~~~~~~~~~~~~~~~
    # Private methods
    # ~~~~~~~~~~~~~~~

    # ------------------------------------------------------------------------
    def _get_psfex_config_filepath(self, sexcat_file_path, file_type, job,
                                   worker):
        """!
        Find and return the .psfex PSFEx cofiguration filepath to use.
        - If the search fails, return @c None.
        """

        psfex_filepath = None

        default_psfex_filename = (worker.config.get_as_string(
                                  'PE_DEFAULT_PSFEX_FILENAME', 'PSFEX'))
        found_files = (self.helper.locate_files(worker,
                       [default_psfex_filename], worker.base_input_dir))
        if len(found_files) > 0:
            psfex_filepath = found_files[0]
        else:
            psfex_filepath = (os.path.abspath(os.path.join(
                              worker.base_input_dir, 'default.psfex')))

            if worker.logging_enabled():
                temp_string = ('{0} - /{1}/img-{2:03}-{3:1d} - {4} - '
                               'Could not find PSFEx .psfex file {5}- '
                               'Trying {6}')
                worker.logger.log_warning_p(temp_string.format(worker.name,
                                            job.get_branch_tree(),
                                            job.img_no, job.epoch,
                                            file_type,
                                            default_psfex_filename,
                                            psfex_filepath))
                worker.logger.flush()

        return psfex_filepath

    # ------------------------------------------------------------------------
    def _get_output_pe_catalog_filename(self, sexcat_filepath, job, config):
        """!
        Build the PE output catalog filename based on the SExtractor catalogue
        filepath
        """

        _, sexcat_filename = os.path.split(sexcat_filepath)
        pe_catalog_prefix, _ = os.path.splitext(sexcat_filename)
        pe_catalog_filename_pattern = (config.get_as_string(
                                       'PE_OUTPUT_CATALOG_FILENAME', 'PSFEX'))
        return pe_catalog_filename_pattern.format(pe_catalog_prefix)
