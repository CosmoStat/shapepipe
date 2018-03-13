# -*- coding: utf-8 -*-

"""EXECUTE MODULE

This module contain a class for executing the package specific code.

:Authors: Samuel Farrens and Marc Gentile

:Date: 31/10/2017 (Happy Halloween!)

"""

# -- Python imports
import os
from time import clock
from shutil import copy

import SETools_script as setools


class PackageRunner(object):

    """Package Runner Class

    This class contains the method run_executable() that runs the relevant
    executable code for this package.

    Parameters
    ----------
    helper : class
        Helper class instance (see helper.py)

    Notes
    -----
    This class is initialised by PackageJobProcessor() in job.py.

    """

    def __init__(self, helper, job, worker):

        self._helper = helper
        self._job = job
        self._worker = worker
        self._fnames = {}
        self.file_types = []
        self._get_file_types()

    def _get_file_types(self):

        """Get File Types

        This method extracts the files types needed to run the executable.

        """

        try:
            self.file_types = (self._worker.config.get_as_list(
                               'INPUT_FILENAME_FORMATS', 'CODE'))

        except Exception:
            pass

        if not self.file_types:
            self.file_types = self._job.get_file_types()

    def run_executable(self):

        """Run Executable

        This method runs the corresponding executable for this package.

        Returns
        -------
        results_dict ditionary of output results

        Notes
        -----
        This method is called by PackageJobProcessor.process_job() in job.py.

        """

        # Dictionary of results
        results_dict = {}

        try:

            # --- Keep track of execution time
            start_time = clock()

            # --- Set input and output files names required to run code
            self._set_filenames()

            # --- Execute the code
            self._exec_code()

            # --- Dictionary with the relevant information for later processing
            results_dict['output_cat_filepath'] = \
                self._fnames['output_filepath']
            results_dict['elapsed_time'] = clock() - start_time

            self._log_output_success(self._fnames['output_filepath_exp'])

        except Exception as detail:

            if self._worker.logging_enabled():
                temp_string = ('{0} - An error occurred while generating '
                               'catalog: {1} ({2})')
                self._worker.logger.log_error_p(temp_string.format(
                                                self._worker.name,
                                                self._job, detail))
                self._worker.logger.flush()

        return results_dict

    def _set_filenames(self):

        """ Set File Names

        This method sets all of the input and output files names needed for
        running the code.

        """

        # --- Input files to be read
        self._fnames['input_filepath'] = [self._job.get_file_path(file_type)
                                          for file_type in self.file_types]

        input_filename = (os.path.splitext(os.path.split(
                          self._fnames['input_filepath'][0])[1])[0])

        # --- Executable configuration file
        self._fnames['config_filepath'] = self._get_exec_config_filepath()

        # --- Target directory where to store files
        output_path = os.path.join(self._worker.result_output_dir,
                                   self._job.get_branch_tree())

        # --- Outpu file name
        self._fnames['output_filepath'] = os.path.join(output_path,
                                                       input_filename)

        # --- Expected output catalog file name
        output_cat_filename_exp = (self._get_output_catalog_filename(
                                   input_filename))

        # --- Output catalog file path
        self._fnames['output_filepath_exp'] = (os.path.abspath(os.path.join(
                                               output_path,
                                               output_cat_filename_exp)))

        self._log_exp_output(self._fnames['output_filepath_exp'])

    def _exec_code(self):

        """Execute the Code

        This method executes the script SETools_script.

        """

        r=setools.SETools(cat_filepath=self._fnames['input_filepath'][0],
                          config_filepath=self._fnames['config_filepath'],
                          output_dir=self._worker.result_output_dir,
                          stat_output_dir=self._worker.stat_output_dir,
                          plot_output_dir=self._worker.plot_output_dir)
        r.process()


    def _get_exec_config_filepath(self):

        """ Get Executable Configuration File

        This method finds and returns the cofiguration file (with full path) to
        use.

        Returns
        -------
        str configuration file name with full path

        """

        default_filename = (self._worker.config.get_as_string(
                            'DEFAULT_FILENAME', 'CODE'))
        found_files = (self._helper.locate_files([default_filename],
                       self._worker.base_input_dir))

        if len(found_files) > 0:

            config_filepath = found_files[0]

            if self._worker.logging_enabled():
                temp_string = ('{0} - /{1}/run-{2:03}-{3:1d} - '
                               'Using configuration file: '
                               '{4}')
                self._worker.logger.log_info_p(temp_string.format(
                                               self._worker.name,
                                               self._job.get_branch_tree(),
                                               self._job.img_no,
                                               self._job.epoch,
                                               config_filepath))
                self._worker.logger.flush()

            # --- Make a copy of the used config file to the log directory
            # for record
            copy(config_filepath, self._worker.log_output_dir)

            return config_filepath

        else:
            if self._worker.logging_enabled():
                temp_string = ('{0} - /{1}/img-{2:03}-{3:1d} - '
                               'Could not find config file {4}')
                self._worker.logger.log_warning_p(temp_string.format(
                                                  self._worker.name,
                                                  self._job.get_branch_tree(),
                                                  self._job.img_no,
                                                  self._job.epoch,
                                                  default_filename))
                self._worker.logger.flush()
            return None

    def _get_output_catalog_filename(self, filepath):

        """ Get Output Catalogue Filename

        This method builds the output catalog file name based on the input
        catalogue file name.

        Parameters
        ----------
        filepath : str
            Input file name with full path

        Returns
        -------
        str output filename with full path

        """

        catalog_prefix = os.path.splitext(os.path.split(filepath)[1])[0]
        output_file_ending = (self._worker.config.get_as_string(
                              'OUTPUT_CATALOG_FILE_ENDING', 'CODE'))

        return catalog_prefix + output_file_ending

    def _log_exp_output(self, file_path):

        """ Update log with expected output catalogue name

        Parameters
        ----------
        file_path : str
            Output file name with full path

        """

        if self._worker.logging_enabled():
            temp_string = ('{0} - /{1}/run-{2:03}-{3:1d} - '
                           'Generating catalog {4}')
            self._worker.logger.log_info_p(temp_string.format(
                                           self._worker.name,
                                           self._job.get_branch_tree(),
                                           self._job.img_no,
                                           self._job.epoch,
                                           file_path))
            self._worker.logger.flush()

    def _log_output_success(self, file_path):

        """ Update log with status of output success

        Parameters
        ----------
        file_path : str
            Output file name with full path

        """

        if self._worker.logging_enabled():
            if os.path.exists(file_path):
                temp_string = ('{0} - /{1}/run-{2:03}-{3:1d} - '
                               'Catalog {4} generated successfully')
                self._worker.logger.log_info_p(temp_string.format(
                                               self._worker.name,
                                               self._job.get_branch_tree(),
                                               self._job.img_no,
                                               self._job.epoch,
                                               file_path))
            else:
                temp_string = ('{0} - An error occurred while generating '
                               'catalog: {1}')
                self._worker.logger.log_error_p(temp_string.format(
                                                self._worker.name,
                                                self._job))
            self._worker.logger.flush()
