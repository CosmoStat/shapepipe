"""UNIT / PROPERTY TESTS FOR THE COVERAGE-MASK FEATURE.

Covers the pure and lightly-fixtured logic behind the coverage nexp masks:
exposure-number parsing, the CCD-list -> unique-exposure reduction, the
MegaCam field-corner convention, the multi-HDU header -> WCS split, the
power-of-two ``nside`` validation, and the ``-h``/argv handling of the
console runners. The heavy paths (VOSpace download, HealSparse map building,
plotting, multiprocessing) are exercised end to end by the pipeline, not here.
"""

import sys

import numpy as np
import numpy.testing as npt
import pytest
from astropy import wcs
from hypothesis import given
from hypothesis import strategies as st

from shapepipe.utilities.coverage_map_builder import CoverageMapBuilder
from shapepipe.utilities.field_corners_extractor import (
    FieldCornersExtractor,
    _expnum_from_path,
    _megacam_field_corners,
    _parse_header_to_wcs,
)
from shapepipe.utilities.header_downloader import HeaderDownloader


def _tan_wcs(crval_ra, crval_dec=0.0):
    """Build a minimal 2-D TAN WCS centred on ``(crval_ra, crval_dec)``."""
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [crval_ra, crval_dec]
    w.wcs.crpix = [1040.0, 1.0]
    w.wcs.cd = [[-1e-5, 0.0], [0.0, 1e-5]]
    return w


# ---------------------------------------------------------------------------
# _expnum_from_path
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "path, expected",
    [
        ("1234567.txt", 1234567),
        ("/a/b/2143523.txt", 2143523),
        ("headers/0000042.txt", 42),
    ],
)
def test_expnum_from_path_extracts_trailing_number(path, expected):
    """The trailing ``<digits>.txt`` is parsed as the exposure number."""
    assert _expnum_from_path(path) == expected


@pytest.mark.parametrize("path", ["no_number.txt", "1234567.fits", "abc.txt"])
def test_expnum_from_path_raises_without_number(path):
    """A filename without a trailing numeric stem raises ``ValueError``."""
    with pytest.raises(ValueError):
        _expnum_from_path(path)


@given(st.integers(min_value=0, max_value=99999999))
def test_expnum_from_path_roundtrips(expnum):
    """Any exposure number round-trips through the filename convention."""
    assert _expnum_from_path(f"vos_headers/{expnum}.txt") == expnum


# ---------------------------------------------------------------------------
# HeaderDownloader.get_exposures
# ---------------------------------------------------------------------------

def test_get_exposures_reduces_to_unique_exposures(tmp_path):
    """A ``<exp>-<ccd>`` CCD list collapses to its unique exposure numbers."""
    ccd_list = tmp_path / "ccds.txt"
    ccd_list.write_text("2143523-0\n2143523-5\n2143524-3\n2143524-8\n")

    exps = HeaderDownloader().get_exposures(str(ccd_list))

    npt.assert_array_equal(exps, np.array([2143523, 2143524]))


def test_get_exposures_csv_matches_txt(tmp_path):
    """The CSV and text code paths yield the same unique exposures."""
    txt = tmp_path / "ccds.txt"
    txt.write_text("2143523-0\n2143523-5\n2143524-3\n")
    csv = tmp_path / "ccds.csv"
    csv.write_text("CCD\n2143523-0\n2143523-5\n2143524-3\n")

    dl = HeaderDownloader()

    npt.assert_array_equal(
        dl.get_exposures(str(txt)), dl.get_exposures(str(csv))
    )


# ---------------------------------------------------------------------------
# field-corner geometry
# ---------------------------------------------------------------------------

def test_megacam_field_corners_uses_expected_ccds():
    """The corners come from CCDs 0, 8, 35, 27, in that order."""
    # Centre each CCD at RA = 100 + index (away from the RA=0 seam) so the
    # returned corner RA reveals which CCD index each slot was drawn from.
    wcs_list = [_tan_wcs(100.0 + i) for i in range(40)]

    ra, dec = _megacam_field_corners(wcs_list)

    assert len(ra) == 4
    assert len(dec) == 4
    npt.assert_allclose(ra, [100, 108, 135, 127], atol=0.5)
    npt.assert_allclose(dec, [0, 0, 0, 0], atol=0.5)


def test_parse_header_to_wcs_returns_one_wcs_per_extension(tmp_path):
    """One WCS is returned per CCD HDU; the primary HDU is skipped."""
    ccd_wcs = [_tan_wcs(10.0), _tan_wcs(20.0), _tan_wcs(30.0)]

    blocks = ["SIMPLE  =                    T"]
    blocks += ["\n".join(str(card) for card in w.to_header().cards)
               for w in ccd_wcs]
    text = "".join(f"{block}\nEND     \n" for block in blocks)

    path = tmp_path / "1234567.txt"
    path.write_text(text)

    result = _parse_header_to_wcs(str(path))

    assert len(result) == len(ccd_wcs)


# ---------------------------------------------------------------------------
# CoverageMapBuilder.check_params (nside validation)
# ---------------------------------------------------------------------------

def test_check_params_accepts_power_of_two_nside(tmp_path):
    """Powers of two for both nside values pass validation."""
    infile = tmp_path / "corners.txt"
    infile.write_text("0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0\n")

    builder = CoverageMapBuilder()
    builder._params["input_file"] = str(infile)
    builder._params["nside_coverage"] = 32
    builder._params["nside"] = 2048

    builder.check_params()


@pytest.mark.parametrize(
    "nside_coverage, nside", [(33, 2048), (32, 100), (32, 3000)]
)
def test_check_params_rejects_non_power_of_two_nside(
    tmp_path, nside_coverage, nside
):
    """A non-power-of-two nside raises ``ValueError``."""
    infile = tmp_path / "corners.txt"
    infile.write_text("0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0\n")

    builder = CoverageMapBuilder()
    builder._params["input_file"] = str(infile)
    builder._params["nside_coverage"] = nside_coverage
    builder._params["nside"] = nside

    with pytest.raises(ValueError):
        builder.check_params()


def test_check_params_missing_input_file_raises(tmp_path):
    """A missing input file raises ``FileNotFoundError``."""
    builder = CoverageMapBuilder()
    builder._params["input_file"] = str(tmp_path / "does_not_exist.txt")

    with pytest.raises(FileNotFoundError):
        builder.check_params()


# ---------------------------------------------------------------------------
# console-runner argv handling (regression guard for the -h entry points)
# ---------------------------------------------------------------------------

def test_run_help_flag_exits_cleanly(monkeypatch):
    """``run()`` with no args reads ``sys.argv[1:]`` so ``-h`` exits 0.

    Guards the regression where the runners parsed the full ``sys.argv``
    (including ``argv[0]``), which made ``-h`` collide with the integer
    ``-c`` option and exit non-zero.
    """
    monkeypatch.setattr(sys, "argv", ["extract_field_corners", "-h"])

    with pytest.raises(SystemExit) as excinfo:
        FieldCornersExtractor().run()

    assert excinfo.value.code == 0
