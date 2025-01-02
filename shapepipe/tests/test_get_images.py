"""UNIT TESTS FOR MODULE PACKAGE: GET_IMAGES.

This module contains unit tests for the module package
shapepipe.modules.get_images_package.get_images

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from unittest import TestCase

import numpy as np
import numpy.testing as npt

from shapepipe.modules.get_images_package import get_images


class GetImagesTestCase(TestCase):

    def setUp(self):

        self.number_tile = "123.456"
        self.number_exp = "2490092p"
        self.number_int = 123456

    def tearDown(self):

        self.number_tile = None
        self.number_exp = None
        self.number_int = None

    def test_in2out_pattern(self):

        npt.assert_string_equal(
            get_images.in2out_pattern(self.number_tile), "123-456"
        )

        npt.assert_string_equal(
            get_images.in2out_pattern(self.number_exp), "2490092"
        )

        npt.assert_raises(
            TypeError,
            get_images.in2out_pattern,
            self.number_int,
        )
