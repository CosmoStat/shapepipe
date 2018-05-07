# -*- coding: utf-8 -*-

"""EXECUTE MODULE

This module contain a class for executing the package specific code.

:Authors: Samuel Farrens, Marc Gentile, Axel Guinot and Morgan Schmitz

:Date: 07/05/2018

"""

# -- Python imports
import os
import re
from time import clock
from shutil import copy

import numpy as np
import scatalog as sc


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

        self._fnames['world_coord'] = self._load_world_coord()

        self._fnames['rect_file'] = self._get_rect_filepath()

        # self._fnames['make_plot'] = self._worker.config.get_as_boolean('MAKE_PLOT', 'CODE')

    def _exec_code(self):

        """Execute the Code

        This method executes the command line defined by _set_exec_line().

        Notes
        -----
        This method need only be modified if you wish to execute the command
        line differently. e.g. using subprocess, etc.

        """

        self.process()

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

    def process(self):
        """Process

        Process the tile to create region's flags file.

        """

        cat = sc.FITSCatalog(self._fnames['input_filepath'][0], SEx_catalog= True)
        cat.open()
        cat_header = cat.get_data(1)[0][0]
        # print(np.array([cat.get_data()[self._fnames['world_coord'][0]], cat.get_data()[self._fnames['world_coord'][1]]]).T)

        try:
            field_coords = np.array([cat.get_data()[self._fnames['world_coord'][0]], cat.get_data()[self._fnames['world_coord'][1]]]).T
        except:
            raise ValueError('World coordinates not found : {}'.format(self._fnames['world_coord']))

        cat.close()

        exprs_names = self._get_exp_name(cat_header)
        rect = self._load_rectangles()
        flag_exp = self._make_flags(field_coords, exprs_names, rect)
        self._save_flags(flag_exp)

    def _load_world_coord(self):
        """Load world coordinates paramaters

        Load the world coordinates to use from catalog.

        Returns
        -------
        list of the world coordinates

        """

        world_params = self._worker.config.get_as_string('WORLD_COORD', 'CODE').split()

        if len(world_params) != 2:
            raise ValueError('You has to provide 2 World coordinates')

        return world_params

    def _get_rect_filepath(self):
        """Load rectangle file

        Load the file containing the coordinates of rectangle defining each single exposures.

        Returns
        -------
        str fullpath to the rectangle file

        """

        file_name = self._worker.config.get_as_string('RECT_FILE', 'CODE')
        config_dir = self._worker.base_input_dir
        full_path = config_dir + '/' + file_name

        if not os.path.isfile(full_path):
            raise ValueError('Rect file not found : {}'.format(full_path))

        return full_path

    def _get_exp_name(self, header):
        """Get exposures name

        Get the name of the single exposures used for the tiles.

        Parameters
        ----------
        header : chararray
            Tile's header from a SExtractor catalog (hdu = 1)

        Returns
        -------
        list of the single exposures

        """

        exp_name = []
        for i in header:
            if i.split()[0] == 'HISTORY':
                exp_name.append(i.split()[3][:-6])

        return exp_name

    def _check_obj_in_exp(self, pos, hori_rect, vert_rect):
        """Check objects in exposures

        Check if an objects belongs to a single exposure.

        Parameters
        ----------
        pos : numpy.ndarray
            Position of the object
        hori_rect : numpy.ndarray
            Limit of the horizontal rectangle
        vert_rect : numpy.ndarray
            Limit of the vertical rectangle

        Returns
        -------
        Boolean
            True if the object is in the single exposure, False otherwise

        """

        belong = False

        if (vert_rect[0,0] < pos[0] < vert_rect[1,0]) and (vert_rect[0,1] < pos[1] < vert_rect[1,1]):
            belong = True
        elif (hori_rect[0,0] < pos[0] < hori_rect[1,0]) and (hori_rect[0,1] < pos[1] < hori_rect[1,1]):
            belong = True

        return belong

    def _load_rectangles(self):
        """Load rectangles

        Load the rectangles of defining each single exposures

        Returns
        -------
        list containing 2 dictionaries for each single exposures defining limits of the rectangles

        """

        try:
            rect = np.load(self._fnames['rect_file'])
        except:
            raise ValueError('Impossible to load rectangles : {}'.format(self._fnames['rect_file']))

        return rect

    def _make_flags(self, field_coords, exp_name, rect):
        """Make flags

        Flag each objects depending of the combination of each single exposures.

        Parameters
        ----------
        field_coords : numpy.ndarray
            Position of each objects
        exp_name : list
            List with the name of each single exposures contributing to the tile
        rect : list
            List containing 2 dictionaries with the position of the rectangles for each single exposures

        Returns
        -------
        numpy.ndarray
            Array with a flag for each object
        """

        nb_obj = field_coords.shape[0]
        flag_exp = np.zeros(nb_obj, dtype= int)
        for i, exp in enumerate(exp_name):
            hori_rect = rect[0][exp]
            vert_rect = rect[1][exp]
            flag_exp += np.array([2**i if self._check_obj_in_exp(pos, hori_rect, vert_rect) else 0 for pos in field_coords])

        for i, flag_val in enumerate(sorted(list(set(flag_exp)))):
            flag_exp[flag_exp == flag_val] = i

        return flag_exp

    def _save_flags(self, exp_flag):
        """Save flags

        Save flags into a numpy file.

        Parameters
        ----------
        exp_flag : numpy.ndarray
            Array with the flag for each object

        """

        s = re.split('-', os.path.split(self._fnames['input_filepath'][0])[1])
        output_name = self._worker.result_output_dir + '/' + s[0] + '_exp_flag' + '-' + s[1] + '-' + s[2]

        np.save(output_name, exp_flag)

    def _make_plot(self):
        """Make plot
        """

        pass
