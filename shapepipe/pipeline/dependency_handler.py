# -*- coding: utf-8 -*-

"""DEPENDENCY HANDLER

This module defines a class for handling pipeline dependencies.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
import re
import subprocess
import importlib


class DependencyHandler(object):
    """ Dependency Handler

    This class manages the required Python packages and system executables
    required to run the pipeline.

    Parameters
    ----------
    dependencies : list
        List of Python packages names, optionally with required versions
    executables : list
        List of system executables

    """

    def __init__(self, dependencies=[], executables=[]):

        self.depend = dependencies
        self.execute = executables
        self._greq = '>='
        self._equal = '=='
        self._great = '>'
        self._less = '<'
        self.dependency_list = []
        self.executable_list = list(set(self.execute))

        if dependencies:
            self._split_strings()
            self._unique_dependencies()

    @property
    def depend(self):
        """ Input Dependency List
        """

        return self._depend

    @depend.setter
    def depend(self, value):

        if not isinstance(value, list):
            raise TypeError('Input must be list type.')

        if not all(isinstance(x, str) for x in value):
            raise ValueError('List elements must be strings.')

        self._depend = value

    @property
    def execute(self):
        """ Input Executable List
        """

        return self._execute

    @execute.setter
    def execute(self, value):

        if not isinstance(value, list):
            raise TypeError('Input must be list type.')

        if not all(isinstance(x, str) for x in value):
            raise ValueError('List elements must be strings.')

        self._execute = value

    @staticmethod
    def _convert_to_float(string):
        """ Convert String to Float

        This method converts numerical strings to floats.

        Parameters
        ----------
        string : str
            Input string

        Returns
        -------
        float
            Converted value

        """

        try:
            val = float(string)
        except Exception:
            val = 0.0

        return val

    @staticmethod
    def _slice_1d(array, indices):
        """ Slice 1D

        Slice 1D list by indices.

        Parameters
        ----------
        array : list
            List of values
        indices : list
            List of inidices

        Returns
        -------
        list
            Sliced list

        """

        return [array[index] for index in indices]

    @classmethod
    def _slice_2d(cls, array, indices):
        """ Slice 2D

        Slice a list of lists by indices.

        Parameters
        ----------
        array : list
            List of lists
        indices : list
            List of inidices

        Returns
        -------
        list
            Sliced list

        """

        return [cls._slice_1d(sublist, indices) for sublist in array]

    @staticmethod
    def _get_indices(array, value):
        """ Get Indices

        Get indices of array elements equal to input value.

        Parameters
        ----------
        array : list
            List of values
        value : str
            Value string

        Returns
        -------
        list
            List of indices

        """

        return [index for index, element in enumerate(array) if
                element == value]

    @classmethod
    def _slice_col_val(cls, array, col, value):
        """ Slice by Column and Value

        Slice a list of lists by elements in a given column equal to a given
        value.

        Parameters
        ----------
        array : list
            List of lists
        col : int
            Column number
        value : str
            Value string

        Returns
        -------
        list
            Slices list

        """

        return cls._slice_2d(array, cls._get_indices(array[col], value))

    @staticmethod
    def _check_executable(exe_name):
        """Check if Input is Executable

        This method checks if the input executable exists.

        Parameters
        ----------
        exe_name : str
            Executable name

        Returns
        -------
        Bool result of test

        Raises
        ------
        TypeError
            For invalid input type

        """

        if not isinstance(exe_name, str):

            raise TypeError('Executable name must be a string.')

        def is_exe(fpath):

            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(exe_name)

        if not fpath:

            res = any([is_exe(os.path.join(path, exe_name)) for path in
                       os.environ["PATH"].split(os.pathsep)])

        else:

            res = is_exe(exe_name)

        if not res:
            raise IOError('{} does not appear to be a valid executable on '
                          'this system.'.format(exe_name))

    def _split_string(self, string):
        """ Split String

        This method splits the version number from the input module string.

        Parameters
        ----------
        string : str
            Input module string

        Returns
        -------
        np.ndarray
            Array of string components

        """

        if self._greq in string:
            val = re.split('({})'.format(self._greq), string)

        elif self._equal in string:
            val = re.split('({})'.format(self._equal), string)

        elif self._great in string:
            val = re.split('({})'.format(self._great), string)

        elif self._less in string:
            raise ValueError('"<" not permitted in package version string.')

        else:
            val = [string, '', '']

        return val

    def _split_strings(self):
        """ Split Strings

        This method splits the input dependency modules strings.

        """

        self._depend_arr = list(map(list, zip(*[self._split_string(string)
                                for string in self.depend])))
        self._dependency_set = set(self._depend_arr[0])

    def _unique_dependencies(self):
        """ Unique Dependencies

        This method creates a unique list of depencies.

        """

        for package_name in self._dependency_set:

            subset = self._slice_col_val(self._depend_arr, 0, package_name)

            if any(self._equal in element for element in subset):

                subset = self._slice_col_val(subset, 1, self._equal)

            if any([ver != '' for ver in subset[2]]):

                subset = (self._slice_col_val(subset, 2,
                          str(max([self._convert_to_float(ver)
                                   for ver in subset[2]]))))

            subset = [element[0] for element in self._slice_2d(subset, [0])]

            self.dependency_list.append(''.join(subset))

    def check_dependencies(self):
        """ Check Dependencies

        This method checks that the required dependencies are installed.

        Returns
        -------
        list
            List of depenecies with versions and paths

        """

        dependency_status_list = []

        for dependency in self._dependency_set:

            try:
                package = importlib.import_module(dependency)
            except Exception:
                raise ImportError('Could not import pipeline dependency '
                                  '{}'.format(dependency))

            if hasattr(package, '__version__'):
                version = package.__version__
            else:
                version = 'N/A'

            if hasattr(package, '__path__'):
                path = package.__path__[0]
            elif hasattr(package, '__file__'):
                path = package.__file__
            else:
                path = 'N/A'

            dependency_status_list.append(' - {} {} {}'.format(
                                          package.__name__,
                                          version, path))

        return dependency_status_list

    def check_executables(self):
        """ Check Executables

        This method checks that the required executables are installed.

        Returns
        -------
        list
            List of executables with paths

        """

        executable_status_list = []

        for executable in self.executable_list:

            self._check_executable(executable)

            exe_path, err = (subprocess.Popen('which {0}'.format(executable),
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE).communicate())

            string = ' - {} {}'.format(executable,
                                       exe_path.rstrip().decode('utf-8'))

            executable_status_list.append(string)

        return executable_status_list
