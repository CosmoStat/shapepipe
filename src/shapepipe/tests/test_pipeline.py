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
import pytest

from shapepipe.pipeline import *
from shapepipe.modules.module_decorator import module_runner


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
            alphabet=(
                st.characters(min_codepoint=65, max_codepoint=90)
                | st.just("_")
            ),
            min_size=1,
        ),
        min_size=1,
        max_size=10,
        unique=True,
    ),
    st.sampled_from(["deg", "arcsec", "pixel"]),
    st.data(),
)
def test_get_unit_from_fits_header_matches_ttype_column(
    column_names,
    unit,
    data,
):

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

    npt.assert_raises(
        IndexError,
        file_io.get_unit_from_fits_header,
        header,
        "DEC",
    )


def test_get_unit_from_fits_header_raises_for_missing_unit():

    header = fits.Header({"TTYPE1": "RA"})

    npt.assert_raises(
        IndexError,
        file_io.get_unit_from_fits_header,
        header,
        "RA",
    )


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


def test_fits_catalogue_table_roundtrips(tmp_path):

    path = tmp_path / "catalogue.fits"
    data = {
        "X": np.array([1.0, 2.0]),
        "FLAG": np.array([0, 1], dtype=np.int16),
    }
    cat = file_io.FITSCatalogue(
        str(path),
        open_mode=file_io.FITSCatalogue.OpenMode.ReadWrite,
    )
    cat.save_as_fits(data=data, ext_name="OBJECTS")

    cat = file_io.FITSCatalogue(str(path))
    cat.open()

    assert cat.get_nb_rows() == 2
    assert cat.get_col_names() == ["X", "FLAG"]
    npt.assert_allclose(cat.get_data()["X"], data["X"])
    npt.assert_array_equal(cat.get_data()["FLAG"], data["FLAG"])

    cat.close()


def test_fits_catalogue_image_roundtrips(tmp_path):

    path = tmp_path / "image.fits"
    image = np.arange(9, dtype=np.float32).reshape(3, 3)
    header = fits.Header({"CRPIX1": 1.0, "CRPIX2": 1.0})
    cat = file_io.FITSCatalogue(
        str(path),
        open_mode=file_io.FITSCatalogue.OpenMode.ReadWrite,
    )
    cat.save_as_fits(
        data=image,
        image=True,
        image_header=header,
        overwrite=True,
    )

    cat = file_io.FITSCatalogue(str(path), hdu_no=0)
    cat.open()

    npt.assert_array_equal(cat.get_data(0), image)
    assert cat.get_header(0)["CRPIX1"] == 1.0

    cat.close()


def test_module_runner_injects_runner_contract_attributes():

    @module_runner(
        version="1.2",
        input_module="sextractor",
        file_pattern=["cat", "weight"],
        file_ext=".fits",
        depends="source-extractor",
        executes=["python"],
        numbering_scheme="-0000000",
        run_method="serial",
    )
    def run(
        input_file_list,
        run_dirs,
        file_number_string,
        config,
        section,
        w_log,
    ):
        return "", ""

    assert run.version == "1.2"
    assert run.input_module == ["sextractor"]
    assert run.file_pattern == ["cat", "weight"]
    assert run.file_ext == [".fits", ".fits"]
    assert run.depends == ["source-extractor"]
    assert run.executes == ["python"]
    assert run.numbering_scheme == "-0000000"
    assert run.run_method == "serial"


def test_module_runner_rejects_mismatched_pattern_and_extension_lengths():

    with pytest.raises(ValueError, match="number of file_ext values"):
        module_runner(
            file_pattern=["cat", "weight"],
            file_ext=[".fits", ".sqlite", ".txt"],
        )


def test_module_runner_rejects_unknown_run_method():

    with pytest.raises(ValueError, match="parallel.*serial"):
        module_runner(run_method="distributed")
