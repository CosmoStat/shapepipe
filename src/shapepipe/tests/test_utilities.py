"""UNIT TESTS FOR UTILITIES.

This module contains unit tests for the shapepipe.pipeline module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from unittest import TestCase

import numpy as np
import numpy.testing as npt

from shapepipe.utilities import galaxy


class GalaxyTestCase(TestCase):

    def setUp(self):

        self.sigma_float = 5.5
        self.sigma_array = np.arange(1, 4) * 0.1
        self.sigma_int = 1
        self.sigma_array_int = np.arange(3)
        self.pixel_scale = 2.0
        self.pixel_int = 1
        self.sigma_float_exp = 12.9515102
        self.sigma_float_ps_exp = 25.9030204
        self.sigma_array_exp = np.array([0.235482, 0.47096401, 0.70644601])

    def tearDown(self):

        self.sigma_float = None
        self.sigma_array = None
        self.sigma_int = None
        self.sigma_array_int = None
        self.pixel_scale = None
        self.pixel_int = None
        self.sigma_float_exp = None
        self.sigma_float_ps_exp = None
        self.sigma_array_exp = None

    def test_sigma_to_fwhm(self):

        npt.assert_almost_equal(
            galaxy.sigma_to_fwhm(self.sigma_float),
            self.sigma_float_exp,
            err_msg="sigma_to_fwhm gave invalid result for float input",
        )

        npt.assert_almost_equal(
            galaxy.sigma_to_fwhm(self.sigma_float, self.pixel_scale),
            self.sigma_float_ps_exp,
            err_msg=(
                "sigma_to_fwhm gave invalid result for float input with "
                + "non-default pixel scale"
            ),
        )

        npt.assert_allclose(
            galaxy.sigma_to_fwhm(self.sigma_array),
            self.sigma_array_exp,
            err_msg="sigma_to_fwhm gave invalid result for array input",
        )

        npt.assert_raises(TypeError, galaxy.sigma_to_fwhm, self.sigma_int)

        npt.assert_raises(
            TypeError,
            galaxy.sigma_to_fwhm,
            self.sigma_array_int,
        )

        npt.assert_raises(
            TypeError,
            galaxy.sigma_to_fwhm,
            self.sigma_float,
            self.pixel_int,
        )

        npt.assert_raises(ValueError, galaxy.sigma_to_fwhm, -self.sigma_float)

        npt.assert_raises(
            ValueError,
            galaxy.sigma_to_fwhm,
            np.negative(self.sigma_array),
        )

        npt.assert_raises(
            ValueError,
            galaxy.sigma_to_fwhm,
            self.sigma_float,
            -self.pixel_scale,
        )
