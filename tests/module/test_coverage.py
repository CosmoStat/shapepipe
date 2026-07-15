"""UNIT / PROPERTY TESTS FOR THE COVERAGE-MASK FEATURE.

Covers the pure and lightly-fixtured logic behind the per-CCD coverage nexp
masks: exposure-number parsing, the CCD-list -> unique-exposure reduction, the
handler's missing-CCD subtraction, per-CCD corner extraction from a multi-HDU
header, ``--ccd_list`` filtering, the builder's per-CCD row parsing and RA-wrap
guard, the power-of-two ``nside`` validation, the atomic-download rename, and
the ``-h``/argv handling of the console runners. The heavy paths (real VOSpace
download, production-resolution map building, plotting, multiprocessing) are
exercised end to end by the pipeline, not here.
"""

import sys
from types import SimpleNamespace
from unittest import mock

import numpy as np
import numpy.testing as npt
import pytest
from astropy import wcs
from astropy.io.fits import Header
from hypothesis import given
from hypothesis import strategies as st

from shapepipe.utilities import ccd_psf_handler as ccd_psf_handler_mod
from shapepipe.utilities.ccd_psf_handler import CcdPsfHandler
from shapepipe.utilities.coverage_map_builder import (
    CoverageMapBuilder,
    unwrap_ra,
)
from shapepipe.utilities.field_corners_extractor import (
    FieldCornersExtractor,
    _ccd_corners,
    _expnum_from_path,
    _parse_header_to_wcs,
)
from shapepipe.utilities.header_downloader import HeaderDownloader


def _tan_wcs(crval_ra, crval_dec=0.0, nx=2080, ny=4612):
    """Build a minimal 2-D TAN WCS with a populated pixel shape.

    ``nx``/``ny`` set ``pixel_shape`` (populated by astropy from
    ``NAXIS1``/``NAXIS2``), so ``_ccd_corners`` can read the CCD pixel bounds.
    """
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [crval_ra, crval_dec]
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.cd = [[-1e-5, 0.0], [0.0, 1e-5]]
    w.pixel_shape = (nx, ny)
    return w


def _header_text(wcs_list):
    """Render a multi-HDU header text file: a primary HDU then one per WCS.

    Each CCD header carries ``NAXIS1``/``NAXIS2`` so the reconstructed WCS has a
    pixel shape.
    """
    blocks = ["SIMPLE  =                    T"]
    for w in wcs_list:
        nx, ny = w.pixel_shape
        header = w.to_header()
        header["NAXIS"] = 2
        header["NAXIS1"] = nx
        header["NAXIS2"] = ny
        blocks.append("\n".join(str(card) for card in header.cards))
    return "".join(f"{block}\nEND     \n" for block in blocks)


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
# HeaderDownloader.get_fits_header — atomic rename
# ---------------------------------------------------------------------------

def test_get_fits_header_writes_atomically(tmp_path):
    """A successful download copies to ``.part`` then renames to the dest."""
    dl = HeaderDownloader()
    dl._params["output_dir"] = str(tmp_path)
    dl._params["overwrite"] = False
    dl._params["dir_for_links"] = None

    dest = tmp_path / "42.txt"
    tmp_dest = tmp_path / "42.txt.part"

    def fake_copy(source, target, head=True):
        # The copy must land on the temp path, not the final destination.
        assert target == str(tmp_dest)
        with open(target, "w") as f:
            f.write("HEADER")

    client = SimpleNamespace(copy=fake_copy)

    assert dl.get_fits_header(42, client) is True
    assert dest.exists()
    assert not tmp_dest.exists()
    assert dest.read_text() == "HEADER"


def test_get_fits_header_failed_copy_leaves_no_dest(tmp_path):
    """A failed download leaves no destination file (only, if any, ``.part``)."""
    dl = HeaderDownloader()
    dl._params["output_dir"] = str(tmp_path)
    dl._params["overwrite"] = False
    dl._params["dir_for_links"] = None

    def failing_copy(source, target, head=True):
        raise RuntimeError("transfer interrupted")

    client = SimpleNamespace(copy=failing_copy)

    assert dl.get_fits_header(42, client) is False
    assert not (tmp_path / "42.txt").exists()
    assert not (tmp_path / "42.txt.part").exists()


# ---------------------------------------------------------------------------
# CcdPsfHandler.get_ccds_with_psf — missing-CCD subtraction
# ---------------------------------------------------------------------------

def test_get_ccds_with_psf_subtracts_missing(monkeypatch):
    """Valid CCDs are all exposure single-HDUs minus the missing set."""
    handler = CcdPsfHandler()

    # Two exposures, 3 CCDs each -> 6 candidate CCDs; two are missing.
    monkeypatch.setattr(
        handler, "get_exp", lambda patches: {"100", "200"}
    )
    monkeypatch.setattr(
        handler,
        "get_exp_shdu_missing",
        lambda patches: {"100-1", "200-2"},
    )
    monkeypatch.setattr(
        ccd_psf_handler_mod.summary,
        "get_all_shdus",
        lambda exposures, n_CCD: [
            f"{e}-{c}" for e in exposures for c in range(n_CCD)
        ],
    )

    result = handler.get_ccds_with_psf(["P1"], n_CCD=3)

    assert result == {"100-0", "100-2", "200-0", "200-1"}


# ---------------------------------------------------------------------------
# per-CCD corner geometry
# ---------------------------------------------------------------------------

def test_ccd_corners_returns_four_corners_around_centre():
    """The 4 corners bracket the CCD centre and cover its full pixel extent."""
    nx, ny = 2080, 4612
    w = _tan_wcs(100.0, 20.0, nx=nx, ny=ny)

    ra, dec = _ccd_corners(w)

    assert len(ra) == 4
    assert len(dec) == 4
    # The footprint straddles the CRVAL centre in both coordinates.
    assert min(ra) < 100.0 < max(ra)
    assert min(dec) < 20.0 < max(dec)
    # Extent matches the CD-scaled pixel span (with edge padding).
    npt.assert_allclose(max(dec) - min(dec), (ny) * 1e-5, rtol=1e-3)


def test_ccd_corners_raises_without_pixel_shape():
    """A WCS with no pixel shape raises a clear ``ValueError``."""
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [100.0, 20.0]
    assert w.pixel_shape is None

    with pytest.raises(ValueError, match="pixel_shape"):
        _ccd_corners(w)


def test_parse_header_to_wcs_returns_one_wcs_per_extension(tmp_path):
    """One WCS is returned per CCD HDU; the primary HDU is skipped."""
    ccd_wcs = [_tan_wcs(10.0), _tan_wcs(20.0), _tan_wcs(30.0)]

    path = tmp_path / "1234567.txt"
    path.write_text(_header_text(ccd_wcs))

    result = _parse_header_to_wcs(str(path))

    assert len(result) == len(ccd_wcs)
    # Reconstructed WCS carry the pixel shape from NAXIS1/NAXIS2.
    assert all(w.pixel_shape == (2080, 4612) for w in result)


# ---------------------------------------------------------------------------
# FieldCornersExtractor.process_single_header — per-CCD rows and filtering
# ---------------------------------------------------------------------------

def test_process_single_header_emits_one_row_per_ccd(tmp_path):
    """Without a CCD list, every HDU yields a ``<expnum>-<idx>`` row."""
    ccd_wcs = [_tan_wcs(10.0), _tan_wcs(20.0), _tan_wcs(30.0)]
    path = tmp_path / "1234567.txt"
    path.write_text(_header_text(ccd_wcs))

    rows = FieldCornersExtractor.process_single_header(
        (str(path), False, None)
    )

    assert [r[0] for r in rows] == [
        "1234567-0",
        "1234567-1",
        "1234567-2",
    ]
    assert all(len(r[1]) == 4 and len(r[2]) == 4 for r in rows)


def test_process_single_header_filters_to_ccd_list(tmp_path):
    """A ``valid_ccds`` set keeps only the listed CCDs of the exposure."""
    ccd_wcs = [_tan_wcs(10.0), _tan_wcs(20.0), _tan_wcs(30.0)]
    path = tmp_path / "1234567.txt"
    path.write_text(_header_text(ccd_wcs))

    rows = FieldCornersExtractor.process_single_header(
        (str(path), False, {"1234567-0", "1234567-2"})
    )

    assert [r[0] for r in rows] == ["1234567-0", "1234567-2"]


def test_load_ccd_list_reads_ids(tmp_path):
    """The CCD list loader returns the stripped, non-blank IDs as a set."""
    path = tmp_path / "ccds.txt"
    path.write_text("100-0\n100-1\n\n200-5\n")

    assert FieldCornersExtractor.load_ccd_list(str(path)) == {
        "100-0",
        "100-1",
        "200-5",
    }


def test_run_extract_writes_per_ccd_rows(tmp_path):
    """End to end: run() writes one per-CCD row filtered by the CCD list."""
    header_dir = tmp_path / "headers"
    header_dir.mkdir()
    ccd_wcs = [_tan_wcs(10.0), _tan_wcs(20.0), _tan_wcs(30.0)]
    (header_dir / "1234567.txt").write_text(_header_text(ccd_wcs))

    ccd_list = tmp_path / "ccds.txt"
    ccd_list.write_text("1234567-0\n1234567-2\n")

    out = tmp_path / "corners.txt"

    extractor = FieldCornersExtractor()
    extractor.run(
        args=[
            "-i", str(header_dir),
            "-l", str(ccd_list),
            "-o", str(out),
        ]
    )

    lines = out.read_text().splitlines()
    assert len(lines) == 2
    ids = [line.split()[0] for line in lines]
    assert ids == ["1234567-0", "1234567-2"]
    # Each row: 1 ID + 4 RA + 4 Dec = 9 columns.
    assert all(len(line.split()) == 9 for line in lines)


def test_get_done_exposures_keys_on_expnum(tmp_path):
    """Resume treats an exposure as done if any of its CCD rows are present."""
    out = tmp_path / "corners.txt"
    out.write_text(
        "1234567-0 10 10 10 10 20 20 20 20\n"
        "1234567-5 10 10 10 10 20 20 20 20\n"
        "7654321-3 30 30 30 30 40 40 40 40\n"
    )

    extractor = FieldCornersExtractor()
    extractor._params["output_file"] = str(out)

    assert extractor.get_done_exposures() == {1234567, 7654321}


# ---------------------------------------------------------------------------
# CoverageMapBuilder.check_params (nside validation)
# ---------------------------------------------------------------------------

def test_check_params_accepts_power_of_two_nside(tmp_path):
    """Powers of two for both nside values pass validation."""
    infile = tmp_path / "corners.txt"
    infile.write_text("100-0 0.0 1.0 1.0 0.0 0.0 0.0 1.0 1.0\n")

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
    infile.write_text("100-0 0.0 1.0 1.0 0.0 0.0 0.0 1.0 1.0\n")

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
# CoverageMapBuilder — RA-wrap guard and per-CCD row build
# ---------------------------------------------------------------------------

def test_unwrap_ra_shifts_straddling_seam():
    """A polygon straddling RA=0 is put on a common (negative) branch."""
    out = unwrap_ra([359.9, 0.1, 0.1, 359.9])
    npt.assert_allclose(sorted(out), [-0.1, -0.1, 0.1, 0.1], atol=1e-6)


def test_unwrap_ra_leaves_normal_polygon_unchanged():
    """A polygon well away from the seam is returned unchanged."""
    npt.assert_allclose(
        unwrap_ra([99.9, 100.1, 100.1, 99.9]),
        [99.9, 100.1, 100.1, 99.9],
    )


def test_build_map_wrapped_ccd_fills_small_patch(tmp_path):
    """A CCD straddling RA=0 fills a small footprint, not the ~360° complement.

    Verified empirically against a same-size CCD away from the seam: with the
    unwrap guard the two footprints have the same pixel count to within a few
    boundary pixels. (Without unwrapping, HealSparse happens to fill the small
    region too for a convex CCD-sized polygon, but the guard makes the intent
    explicit and is robust to degenerate orderings.)
    """
    # CCD straddling RA=0 (0.2 x 0.2 deg) and a reference CCD at RA=10.
    infile = tmp_path / "corners.txt"
    infile.write_text(
        "100-0 359.9   0.1   0.1 359.9 0.0 0.0 0.2 0.2\n"
        "200-0   9.9  10.1  10.1   9.9 0.0 0.0 0.2 0.2\n"
    )
    out = tmp_path / "cov.hsp"

    builder = CoverageMapBuilder()
    builder.run(
        args=[
            "-i", str(infile),
            "-o", str(out),
            "-c", "32",
            "-n", "1024",
        ]
    )

    import healsparse as hsp

    m = hsp.HealSparseMap.read(str(out))
    n_valid = len(m.valid_pixels)

    # Two ~0.04 deg^2 patches at nside=1024 (~3.4 arcmin pixels): a small,
    # bounded number of pixels — nowhere near a hemisphere.
    assert 0 < n_valid < 200


def test_build_map_pole_guard_skips_polygon(tmp_path, capsys):
    """A polygon with |dec| near 90 deg is skipped with a warning."""
    infile = tmp_path / "corners.txt"
    infile.write_text(
        "100-0  10.0  10.2  10.2  10.0 89.5 89.5 89.7 89.7\n"
        "200-0  10.0  10.2  10.2  10.0 20.0 20.0 20.2 20.2\n"
    )
    out = tmp_path / "cov.hsp"

    builder = CoverageMapBuilder()
    builder.run(
        args=["-i", str(infile), "-o", str(out), "-c", "32", "-n", "1024"]
    )

    captured = capsys.readouterr()
    assert "skipping CCD 100-0" in captured.out
    assert "Added 1 polygons" in captured.out


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
