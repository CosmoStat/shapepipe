"""PYTHON EXAMPLE.

This module contains an example Python class.

:Author: Samuel Farrens

.. warning::
  This is a example module. It should only be used for running tests or as a
  reference for developing new modules.

"""

import time

from numpy.random import randint


class PythonExample:
    """Python Example.

    An example Python class.

    Parameters
    ----------
    sleep_time : int
        Sleep time in seconds

    """

    def __init__(self, sleep_time=None):

        if not isinstance(sleep_time, type(None)):
            self.sleep_time = sleep_time

        else:
            self.sleep_time = randint(1, 10)

    def _wait(self):
        """Wait.

        Wait for :math:`n` seconds.

        """
        time.sleep(self.sleep_time)

    def _read_file(self, file_name):
        """Read File.

        Read input file content.

        Parameters
        ----------
        file_name : str
            Name of file to read

        Returns
        -------
        str
            Content of file

        """
        with open(file_name) as data_file:
            content = data_file.read().replace("\n", "")

        return content

    def read_files(self, file_name1, file_name2):
        """Read Files.

        Read two input files.

        Parameters
        ----------
        file_name1 : str
            Name of first file
        file_name2 : str
            Name of second file

        """
        self._wait()
        content1 = self._read_file(file_name1)
        content2 = self._read_file(file_name2)

        self.content = f"{content1} and {content2}"

    def write_file(self, file_name, message):
        """Write File.

        Write content to file.

        Parameters
        ----------
        file_name : str
            Name of output file
        Message : str
            Content to write to file

        """
        new_content = message + str(self.content)

        text_file = open(file_name, "w")
        text_file.write(new_content)
        text_file.close()
