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

import numpy as np
import re
import scatalog as sc
from math import ceil, floor


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
        self._fnames['col_name'] = (self._worker.config.get_as_list(
                                        'COL_NAME', 'CODE'))

        # -- Load an extra code option from config file
        self._fnames['vign_size'] = (self._worker.config.get_as_list(
                                        'VIGN_SIZE', 'CODE'))


    def _exec_code(self):

        """Execute the Code

        This method executes the command line defined by _set_exec_line().

        Notes
        -----
        This method need only be modified if you wish to execute the command
        line differently. e.g. using subprocess, etc.

        """

        image = self._load_data(self._fnames['input_filepath'][0], image= True)
        cat = self._load_data(self._fnames['input_filepath'][1], SEx_cat= True)

        pos = [[x,y] for x,y in zip(cat[self._fnames['col_name'][0]], cat[self._fnames['col_name'][1]])]

        vign = self._extract_vignet(image= image, pos= pos, size= self._fnames['vign_size'])
        self._save(vign)

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

    def _load_data(self, file_path, SEx_cat= False, image= False):
        """Load data

        Load the image and the catalog (assume a SExtractor catalog for position)

        Parameters
        ----------
        file_path : str
            Path to the image/catalog
        SEx_cat : bool
            If a SExtractor catalog input set 'True'
        image : bool
            If an image input set 'True'

        """
        if SEx_cat:
            f = sc.FITSCatalog(file_path, SEx_catalog= True)
        elif image:
            f = sc.FITSCatalog(file_path, hdu_no= 0)
        else:
            raise ValueError("SEx_cat or image not set. One of them has to be set on 'True' in accordingly to the file to open.")

        f.open()
        data = f.get_data()

        f.close()

        return data

    def _extract_vignet(self, image, pos, size):
        """Extract vignet

        Create vignets at the positions given with the given size

        Parameters
        ----------
        image : nump.ndarray
            Array containing the image
        pos : list
            List of position with the following format : [[x0, y0], [x1, y1]]
        size : list
            List of the vignet size : [size_x, size_y]

        Returns
        -------
        numpy.ndarray
            Array containing the vignets

        """

        # p : parity
        # d : delta
        # ru : delta right/up
        # ld : delta left/down
        # rs : resize
        p = [0, 0]
        d = [0, 0]
        rd = [0, 0]
        lu = [0, 0]
        self._check_pos_rs = [0, 0]
        im_shape = image.shape

        if (size[0] > im_shape[0]) | (size[1] > im_shape[1]):
            raise ValueError('Vignet size : {0} must be smaller than the image size : {1}'.format(tuple(size), im_shape))

        for i in range(2):
            if size[i]%2 == 0:
                p[i] = 1
                d[i] = int(size[i]/2.)
            else:
                p[i] = 0
                d[i] = int((size[i]-1)/2.)

        vign = []
        for pos_temp in pos:
            for j in range(2):
                if p[j] == 1:
                    if pos_temp[j] % 1 < 0.5:
                        rd[j] = self._check_pos(int(pos_temp[j] // 1 + d[j]), im_shape, j)
                        lu[j] = self._check_pos(int(pos_temp[j] // 1 - d[j]), im_shape, j)
                    else:
                        rd[j] = self._check_pos(int(pos_temp[j] // 1 + 1 + d[j]), im_shape, j)
                        lu[j] = self._check_pos(int(pos_temp[j] // 1 + 1 - d[j]), im_shape, j)
                else:
                    rd[j] = self._check_pos(int(ceil(pos_temp[j]) + d[j]), im_shape, j)
                    lu[j] = self._check_pos(int(floor(pos_temp[j]) - d[j]), im_shape, j)

            if self._check_pos_rs == [0, 0]:
                vign.append(image[lu[1]:rd[1],lu[0]:rd[0]])
            else:
                tmp_vign = np.zeros(tuple(size))
                if self._check_pos_rs[0] == -1:
                    if self._check_pos_rs[1] == 0:
                        tmp_vign[:,size[0]-rd[0]-lu[0]:] = image[lu[1]:rd[1],lu[0]:rd[0]]
                        vign.append(tmp_vign)
                    elif self._check_pos_rs[1] == -1:
                        tmp_vign[size[1]-lu[1]-rd[1]:,size[0]-rd[0]-lu[0]:] = image[lu[1]:rd[1],lu[0]:rd[0]]
                        vign.append(tmp_vign)
                    elif self._check_pos_rs[1] == 1:
                        tmp_vign[:rd[1]-lu[1],size[0]-rd[0]-lu[0]:] = image[lu[1]:rd[1],lu[0]:rd[0]]
                        vign.append(tmp_vign)
                elif self._check_pos_rs[0] == 1:
                    if self._check_pos_rs[1] == 0:
                        tmp_vign[:,:rd[0]-lu[0]] = image[lu[1]:rd[1],lu[0]:rd[0]]
                        vign.append(tmp_vign)
                    elif self._check_pos_rs[1] == -1:
                        tmp_vign[size[1]-lu[1]-rd[1]:,:rd[0]-lu[0]] = image[lu[1]:rd[1],lu[0]:rd[0]]
                        vign.append(tmp_vign)
                    elif self._check_pos_rs[1] == 1:
                        tmp_vign[:rd[1]-lu[1],:rd[0]-lu[0]] = image[lu[1]:rd[1],lu[0]:rd[0]]
                        vign.append(tmp_vign)

        return np.array(vign)

    def _check_pos(self, pos, im_shape, axis):
        """Check position

        Check if the vignet is entirely in the image

        Parameters
        ----------
        pos : float
            Position to check
        im_shape : list
            Shape of the image
        axis : int
            Axis to check in [0, 1]

        Returns
        -------
        float
            Return the nearest border position or the inputed position if it is inside the image

        """

        if pos < 0:
            self._check_pos_rs[axis] = -1
            return 0
        elif pos > im_shape[axis]:
            self._check_pos_rs[axis] = 1
            return im_shape[axis]
        else:
            return pos

    def _save(self, vign):
        """Save

        Save the vignets array in a SExtractor fits catalog format

        Parameters
        ----------
        vign : numpy.ndarray
            Array containing the vignets

        """

        s=re.split("\-([0-9]{3})\-([0-9]+)\.",self._fnames['input_filepath'][1])
        output_path = self._worker.result_output_dir + '/' + 'vignets-{0}-{1}.fits'.format(s[1],s[2])

        f = sc.FITSCatalog(output_path, SEx_catalog= True, open_mode= sc.BaseCatalog.OpenMode.ReadWrite)
        f.save_as_fits(vign, ['VIGNET'], sex_cat_path= self._fnames['input_filepath'][1])
