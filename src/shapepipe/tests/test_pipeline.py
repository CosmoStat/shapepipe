"""UNIT TESTS FOR PIPELINE.

This module contains unit tests for the shapepipe.pipeline module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from unittest import TestCase

from astropy.io import fits
from astropy import units
from hypothesis import given
from hypothesis import strategies as st
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


@given(
    st.lists(
        st.text(
            alphabet=st.characters(min_codepoint=65, max_codepoint=90) | st.just("_"),
            min_size=1,
        ),
        min_size=1,
        max_size=10,
        unique=True,
    ),
    st.sampled_from(["deg", "arcsec", "pixel"]),
    st.data(),
)
def test_get_unit_from_fits_header_matches_ttype_column(column_names, unit, data):

    key = data.draw(st.sampled_from(column_names))
    header = fits.Header(
        {
            card: value
            for idx, column_name in enumerate(column_names, start=1)
            for card, value in (
                (f"TTYPE{idx}", column_name),
                (f"TCUNI{idx}", unit),
            )
        }
    )

    assert file_io.get_unit_from_fits_header(header, key) == units.Unit(unit)


def test_get_unit_from_fits_header_raises_for_missing_column():

    header = fits.Header({"TTYPE1": "RA", "TCUNI1": "deg"})

    npt.assert_raises(IndexError, file_io.get_unit_from_fits_header, header, "DEC")


def test_get_unit_from_fits_header_raises_for_missing_unit():

    header = fits.Header({"TTYPE1": "RA"})

    npt.assert_raises(IndexError, file_io.get_unit_from_fits_header, header, "RA")


def test_custom_parser_getlist_expands_env_and_strips_entries(monkeypatch):

    monkeypatch.setenv("SP_TEST_ROOT", "/tmp/shapepipe")
    parser = config.CustomParser()
    parser.add_section("PATHS")
    parser.set("PATHS", "FILES", "$SP_TEST_ROOT/a.fits, $SP_TEST_ROOT/b.fits")

    assert parser.getlist("PATHS", "FILES") == [
        "/tmp/shapepipe/a.fits",
        "/tmp/shapepipe/b.fits",
    ]


def test_custom_parser_getlist_honours_custom_delimiter():

    parser = config.CustomParser()
    parser.add_section("VALUES")
    parser.set("VALUES", "ITEMS", "alpha| beta |gamma")

    assert parser.getlist("VALUES", "ITEMS", delimiter="|") == [
        "alpha",
        "beta",
        "gamma",
    ]
