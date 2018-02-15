# -*- coding: utf-8 -*-

"""EXECUTE MODULE

This module contain a class for executing the package specific code.

:Authors: Samuel Farrens, Marc Gentile and Axel Guinot

:Date: 18/01/2018

"""

# -- Python imports
import os
from time import clock
from shutil import copy

import scatalog as sc
import re


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


        self._fnames['SEXTRACTOR'] = {'input_dir' : self._worker.log_output_dir + '/sextractor_config_inputs/'}

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

            # --- Set the execution command for the code
            self._set_exec_line()

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

        """Set File Names

        This method sets all of the input and output files names needed for
        running the code.

        """

        # --- Input files to be read
        self._fnames['input_filepath'] = [self._job.get_file_path(file_type)
                                          for file_type in self.file_types]

        input_filename = (os.path.splitext(os.path.split(
                          self._fnames['input_filepath'][0])[1])[0])

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

        # -- Load an extra code option from config file
        # self._fnames['extra_option'] = (self._worker.config.get_as_string(
        #                                 'EXTRA_CODE_OPTION', 'CODE'))

    def _set_exec_line(self):

        """Set the Command Line to be Executed

        This method defines the command line for the code corresponding to this
        package.

        Notes
        -----
        In the first time it will create config files needed for SExtractor
        (*.sex and *.param). Then setup the command line.

        """

        # --- Execution line
        exec_path = self._worker.config.get_as_string('EXEC_PATH', 'CODE')

        # --- SExtractor configuration
        self._config_sex_input()
        self._config_sex_output()

        self._exec_line = ('{0} {1} -c {2}').format(
                           exec_path,
                           self._fnames['input_filepath'][0],
                           self._fnames['SEXTRACTOR']['input_config_path'])

        self._log_exec_line()


    def _exec_code(self):

        """Execute the Code

        This method executes the command line defined by _set_exec_line().

        """

        os.system(self._exec_line)


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
        str
            Output filename with full path

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

    def _log_exec_line(self):

        """ Update log with expected output catalogue name

        Parameters
        ----------
        exec_line : str
            Execution line to be run

        """

        if self._worker.logging_enabled():
            temp_string = ('{0} - /{1}/run-{2:03}-{3:1d} - '
                           'Executing command: {4}')
            self._worker.logger.log_info_p(temp_string.format(
                                           self._worker.name,
                                           self._job.get_branch_tree(),
                                           self._job.img_no,
                                           self._job.epoch,
                                           self._exec_line))
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


    def _config_sex_input(self):

        """ Create the *.sex input file for SEXTRACTOR

        This function mix parameters from the package config file with thus
        provide in the default config file.

        Note
        ----
        Parameters can be set using header's parameters value (possibility to
        to make operation).

        """

        params = self._worker.config.get_section_data('SEXTRACTOR_INPUT')
        if params is None:
            return None

        self._sex_input_params = {}


        s=re.split("\-([0-9]{3})\-([0-9]+)\.", self._fnames['input_filepath'][0])
        img_num = '-{0}-{1}'.format(s[1],s[2])
        self._sex_input_params['CATALOG_NAME'] = self._worker.result_output_dir + '/' + re.split('\.', self.file_types[0])[0] + img_num + '.cat'

        self._extra_file()

        self._sex_input_params['PARAMETERS_NAME'] = self._fnames['SEXTRACTOR']['input_dir'] + 'default.param'

        f = sc.FITSCatalog(self._fnames['input_filepath'][0], hdu_no = 0)
        f.open()

        for i in params.keys():
            if i in ['WEIGHT_IMAGE', 'FLAG_IMAGE', 'PSF_NAME', 'ASSOC_NAME', 'PARAMETERS_NAME']:
                continue
            ds_param = re.split('@',params[i])
            if len(ds_param) == 1:
                self._sex_input_params[i] = params[i]
            elif len(ds_param) == 2:
                self._sex_input_params[i] = f.get_header_value(ds_param[1])
            else:
                raise Exception('Value in a wrong format for key : {0}'.format(i))

        self._set_default_input()
        f.close()
        del(f)
        self._save_sex_input_config()


    def _extra_file(self):
        """Configure extra input

        This functions set coresponding parameters to use optional input file.

        Note
        ----
        Allow one to provide weight and/or flag images.
        Handle *.psf file from PSFEx.
        Allow one to provide association catalog.

        """

        if len(self._fnames['input_filepath']) == 1:
            return None

        i = 1
        file_type1 = ['WEIGHT', 'FLAG']
        file_type2 = ['PSF', 'ASSOC']
        while i+1 <= len(self._fnames['input_filepath']):
            for j in file_type1 + file_type2:
                if self._worker.config.get_as_boolean(j, 'CODE'):
                    if j in file_type1:
                        self._sex_input_params['{0}_IMAGE'.format(j)] = self._fnames['input_filepath'][i]
                        file_type1.remove(j)
                        break
                    elif j in file_type2:
                        self._sex_input_params['{0}_NAME'.format(j)] = self._fnames['input_filepath'][i]
                        file_type2.remove(j)
                        break
            i+=1

        # Consistency checks on flags and file types
        want_WEIGHT  = self._worker.config.get_as_boolean('WEIGHT', 'CODE')
        want_FLAG    = self._worker.config.get_as_boolean('FLAG', 'CODE')
        n_file_types = len(self._worker.config.get_as_list('INPUT_FILENAME_FORMATS', 'CODE'))
        if int(want_WEIGHT) + int(want_FLAG) != n_file_types - 1:
            raise Exception('Boolean flags in package config file WEIGHT={} and FLAG={} not consistent '
                            'with number of input file types {}'.format(want_WEIGHT, want_FLAG, n_file_types))


    def _set_default_input(self):
        """Set default input parameters

        This function read the default.sex file and fill not specified
        parameters to run SExtractor

        """

        dir_path = self._worker.base_input_dir + 'SExtractor_default'

        self._read_sex_default(dir_path + '/default.sex')

        for i in self._default_input_params.keys():
            if not i in self._sex_input_params.keys():
                self._sex_input_params[i] = self._default_input_params[i]


    def _read_sex_default(self, path):
        """Read default.sex

        This function read the default.sex file and put them in a dictionary.

        Parameters
        ----------
        path : str
            Path to default.sex file.

        """

        self._default_input_params = {}
        f=open(path, 'r')

        while True:
            line = f.readline()

            if line == '':
                break

            line = self._clean_line(line)

            if line != None:
                s = re.split('\s+', line)
                if len(s) >= 2:
                    self._default_input_params[s[0]] = ''.join(s[1:])

        f.close()


    def _clean_line(self, line_tmp):
        """Clean Lines

        This function is called during the reading process to clean empty lines
        and ignore comments.

        Parameters
        ----------
        line_tmp : str
            The string to clean up.

        Returns
        -------
        str
            If the line is not empty or a comment return the contents and
            None otherwise.

        """

        if re.split('#',line_tmp)[0] == '':
            return None

        line_tmp = line_tmp.replace('\n','')
        line_tmp = re.split('#', line_tmp)[0]

        if line_tmp != '':
            return line_tmp
        else:
            return None


    def _save_sex_input_config(self):
        """Save the input *.sex file

        This function write the *.sex file for this run and save it in the log
        output directory.

        Note
        ----
        Name of the file : config-***-*.sex

        """

        dir_path = self._fnames['SEXTRACTOR']['input_dir']
        s=re.split("\-([0-9]{3})\-([0-9]+)\.", self._fnames['input_filepath'][0])
        img_num = '-{0}-{1}'.format(s[1],s[2])
        self._fnames['SEXTRACTOR']['input_config_path'] = dir_path + '/config{0}.sex'.format(img_num)

        if not os.path.isdir(dir_path):
            try:
                os.system('mkdir {0}'.format(dir_path))
            except:
                raise Exception('Impossible to create the directory for sextractor input config files in {0}'.format(self._worker.log_output_dir))

        self._file_dependencies()

        f = open(self._fnames['SEXTRACTOR']['input_config_path'],'w')

        if len(self._sex_input_params) == 0:
            raise ValueError('No parameters specified for the SExtractor input config file')

        for i in self._sex_input_params.keys():
            f.write('{0}\t{1}\n'.format(i, self._sex_input_params[i]))

        f.close()

    def _config_sex_output(self):
        """Configure *.param output config for SExtractor

        This function will create the *.param file from parameters provide in
        the package input config file.

        Note
        ----
        At this point it's not possible to provide an externat *.param file.

        """

        if os.path.isfile(self._fnames['SEXTRACTOR']['input_dir'] + 'default.param'):
            return None

        params = self._worker.config.get_section_data('SEXTRACTOR_OUTPUT')
        self._sex_output_params = []

        if params is None:
            return None

        for i in params.keys():
            self._sex_output_params.append(i)

        self._save_sex_output_config()

    def _save_sex_output_config(self):
        """Save the output *.param file

        This function write the *.param file for this run and save it in the log
        output directory.

        """

        dir_path = self._fnames['SEXTRACTOR']['input_dir']
        filename = 'default.param'
        self._fnames['SEXTRACTOR']['output_config_path'] = dir_path + filename

        f = open(dir_path + filename, 'w')

        for i in self._sex_output_params:
            f.write('{0}\n'.format(i))

        f.close()

    def _file_dependencies(self):
        """Handle file dependencies

        This function copy and set parameters corresponding to input file
        necessary for filtering and neural-network.

        Note
        ----
        The file are copied to the output log directory.

        """

        for i in ['FILTER_NAME', 'STARNNW_NAME']:
            try:
                copy(self._worker.base_input_dir + '/SExtractor_default/' + self._sex_input_params[i], self._worker.log_output_dir + '/sextractor_config_inputs/')
                self._sex_input_params[i] = self._worker.log_output_dir + '/sextractor_config_inputs/' + self._sex_input_params[i]
            except:
                pass
