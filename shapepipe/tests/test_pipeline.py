"""UNIT TESTS FOR PIPELINE.

This module contains unit tests for the shapepipe.pipeline module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from unittest import TestCase

import numpy as np
import numpy.testing as npt

from shapepipe.pipeline import *


class ExecuteTestCase(TestCase):

    def setUp(self):

        self.command_line = "echo 1"
        self.output_tuple = ("1\n", "")

    def tearDown(self):

        self.command_line = None
        self.output_tuple = None

    def test_execute(self):

        npt.assert_raises(TypeError, execute.execute, 1)
        npt.assert_equal(execute.execute(self.command_line), self.output_tuple)

    def test_check_executable(self):

        npt.assert_raises(TypeError, execute.check_executable, 1)
        npt.assert_raises(OSError, execute.check_executable, "")
        self.assertIsNone(execute.check_executable("/bin/ls"))
