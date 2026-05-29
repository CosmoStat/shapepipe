"""UNIT TESTS FOR MODULE PACKAGE: GET_IMAGES.

This module contains unit tests for the module package
shapepipe.modules.get_images_package.get_images

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from unittest import TestCase

from hypothesis import given
from hypothesis import strategies as st
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

        npt.assert_string_equal(get_images.in2out_pattern(self.number_tile), "123-456")

        npt.assert_string_equal(get_images.in2out_pattern(self.number_exp), "2490092")

        npt.assert_raises(
            TypeError,
            get_images.in2out_pattern,
            self.number_int,
        )


@given(st.text(alphabet=st.characters(blacklist_categories=["Cs"])))
def test_in2out_pattern_is_idempotent(number):

    npt.assert_string_equal(
        get_images.in2out_pattern(get_images.in2out_pattern(number)),
        get_images.in2out_pattern(number),
    )


@given(st.text(alphabet=st.characters(blacklist_categories=["Cs"])))
def test_in2out_pattern_normalizes_filename_separators(number):

    transformed = get_images.in2out_pattern(number)

    assert "." not in transformed
    assert "_" not in transformed
