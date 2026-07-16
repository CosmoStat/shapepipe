"""UNIT / PROPERTY TESTS FOR THE COVERAGE-MASK FEATURE.

Covers the pure and lightly-fixtured logic behind the per-CCD coverage nexp
masks: exposure-number parsing, the CCD-list -> unique-exposure reduction, the
handler's missing-CCD subtraction (pinned against the real ID format), per-CCD
corner extraction from plain and fpack-compressed multi-HDU headers,
``--ccd_list`` filtering, the per-CCD resume path, the builder's per-CCD row
parsing, the accumulated exposure-count (nexp) contract, the RA-wrap and pole
guards, the power-of-two ``nside`` validation, the atomic-download rename, and
the ``-h``/argv handling of the console runners. The heavy paths (real VOSpace
download, production-resolution map building, plotting, multiprocessing) are
exercised end to end by the pipeline, not here.
"""

import sys
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pytest
from astropy import wcs
from hypothesis import given
from hypothesis import strategies as st

from shapepipe.utilities import summary
from shapepipe.utilities.ccd_psf_handler import CcdPsfHandler
from shapepipe.utilities.coverage_map_builder import (
    CoverageMapBuilder,
    unwrap_ra,
)
from shapepipe.utilities.field_corners_extractor import (
    FieldCornersExtractor,
    _ccd_corners,
    _expnum_from_path,
    _image_shape,
    _parse_header_to_wcs,
)
from shapepipe.utilities.header_downloader import HeaderDownloader


def _tan_wcs(crval_ra, crval_dec=0.0, nx=2080, ny=4612):
    """Build a minimal 2-D TAN WCS with a populated pixel shape.

    ``nx``/``ny`` set ``pixel_shape`` so ``_header_text`` can emit matching
    ``NAXIS1``/``NAXIS2`` keywords.
    """
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [crval_ra, crval_dec]
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.cd = [[-1e-5, 0.0], [0.0, 1e-5]]
    w.pixel_shape = (nx, ny)
    return w


def _header_text(wcs_list):
    """Render a multi-HDU header text file: a primary HDU then one plain image
    HDU per WCS, each carrying ``NAXIS1``/``NAXIS2``.
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


def _compressed_header_text(w):
    """Render a one-CCD header mimicking an fpack tile-compressed HDU.

    ``NAXIS1``/``NAXIS2`` describe the compressed binary table (byte width, row
    count); the true image dimensions live in ``ZNAXIS1``/``ZNAXIS2``. This is
    the shape of headers fetched from '<expnum>p.fits.fz' with ``head=True``.
    """
    nx, ny = w.pixel_shape
    wcs_cards = "\n".join(str(card) for card in w.to_header().cards)
    ccd = (
        "XTENSION= 'BINTABLE'\n"
        "BITPIX  =                    8\n"
        "NAXIS   =                    2\n"
        "NAXIS1  =                    8\n"
        f"NAXIS2  =                 {ny}\n"
        "PCOUNT  =              1000000\n"
        "GCOUNT  =                    1\n"
        "TFIELDS =                    1\n"
        "ZIMAGE  =                    T\n"
        "ZBITPIX =                  -32\n"
        "ZNAXIS  =                    2\n"
        f"ZNAXIS1 =                 {nx}\n"
        f"ZNAXIS2 =                 {ny}\n"
        f"{wcs_cards}"
    )
    return f"SIMPLE  =                    T\nEND     \n{ccd}\nEND     \n"


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


def test_get_exposures_single_line(tmp_path):
    """A one-line CCD list (0-d loadtxt array) still yields one exposure."""
    ccd_list = tmp_path / "ccds.txt"
    ccd_list.write_text("2143523-0\n")

    exps = HeaderDownloader().get_exposures(str(ccd_list))

    npt.assert_array_equal(exps, np.array([2143523]))


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
    """Valid CCDs are all exposure single-HDUs minus the missing set.

    The real ``summary.get_all_shdus`` is used so the cross-component
    ``<exp>-<ccd>`` ID format is pinned end to end.
    """
    handler = CcdPsfHandler()

    # Two exposures, 3 CCDs each -> 6 candidate CCDs; two are missing.
    monkeypatch.setattr(handler, "get_exp", lambda patches: {"100", "200"})
    monkeypatch.setattr(
        handler,
        "get_exp_shdu_missing",
        lambda patches: {"100-1", "200-2"},
    )

    result = handler.get_ccds_with_psf(["P1"], n_CCD=3)

    # get_all_shdus yields "<exp>-<ccd>" for ccd in range(n_CCD).
    assert result == {"100-0", "100-2", "200-0", "200-1"}
    # Guard the assumption that the missing IDs share the produced format.
    assert set(summary.get_all_shdus({"100"}, 3)) == {"100-0", "100-1", "100-2"}


@pytest.mark.parametrize(
    ("version", "n_patch"),
    [("v1.3", 7), ("v1.4", 7), ("v1.5", 8), ("v1.6", 9), ("v2.0", 10)],
)
def test_version_to_patch_count(version, n_patch):
    """Each catalogue version maps to its patch count; v2.0 mirrors v1.6."""
    handler = CcdPsfHandler()
    handler._params["version_cat"] = version
    handler.update_params()
    assert handler._params["n_patch"] == n_patch
    assert len(handler._params["patches"]) == n_patch


def test_n_patch_option_overrides_version_default():
    """An explicit -p patch count wins over the version mapping."""
    handler = CcdPsfHandler()
    handler._params["version_cat"] = "v2.0"
    handler._params["n_patch"] = 4
    handler.update_params()
    assert handler._params["patches"] == ["P1", "P2", "P3", "P4"]


def test_invalid_version_raises():
    """An unknown catalogue version fails loudly."""
    handler = CcdPsfHandler()
    handler._params["version_cat"] = "v9.9"
    with pytest.raises(ValueError, match="v9.9"):
        handler.update_params()


# ---------------------------------------------------------------------------
# image-shape resolution (fpack ZNAXIS vs plain NAXIS)
# ---------------------------------------------------------------------------

def test_image_shape_prefers_znaxis_for_compressed_header(tmp_path):
    """A compressed HDU reports ZNAXIS dims, not the binary-table NAXIS."""
    w = _tan_wcs(100.0, 20.0, nx=2080, ny=4612)
    path = tmp_path / "1234567.txt"
    path.write_text(_compressed_header_text(w))

    (parsed_w, shape), = _parse_header_to_wcs(str(path))

    # WCS pixel_shape would wrongly report the compressed byte width (8).
    assert parsed_w.pixel_shape == (8, 4612)
    # _image_shape recovers the true image dimensions.
    assert shape == (2080, 4612)


def test_image_shape_falls_back_to_naxis(tmp_path):
    """A plain image HDU (no ZIMAGE) uses NAXIS1/NAXIS2."""
    w = _tan_wcs(100.0, 20.0, nx=2080, ny=4612)
    path = tmp_path / "1234567.txt"
    path.write_text(_header_text([w]))

    (_, shape), = _parse_header_to_wcs(str(path))

    assert shape == (2080, 4612)


def test_image_shape_raises_without_dimensions():
    """A header with no image dimensions raises a clear ``ValueError``."""
    from astropy.io.fits import Header

    header = Header()
    header["CTYPE1"] = "RA---TAN"

    with pytest.raises(ValueError, match="image dimensions"):
        _image_shape(header)


def test_compressed_header_corners_are_full_width(tmp_path):
    """Corners from a compressed header span the true CCD width, not 8 px."""
    w = _tan_wcs(100.0, 20.0, nx=2080, ny=4612)
    path = tmp_path / "1234567.txt"
    path.write_text(_compressed_header_text(w))

    (parsed_w, shape), = _parse_header_to_wcs(str(path))
    ra, dec = _ccd_corners(parsed_w, shape)

    # RA extent must reflect ~2080 px * 1e-5 deg/px * cos(dec), not 8 px.
    ra_extent = max(ra) - min(ra)
    assert ra_extent > 0.01


# ---------------------------------------------------------------------------
# per-CCD corner geometry
# ---------------------------------------------------------------------------

def test_ccd_corners_returns_four_corners_around_centre():
    """The 4 corners bracket the CCD centre and cover its full pixel extent."""
    nx, ny = 2080, 4612
    w = _tan_wcs(100.0, 20.0, nx=nx, ny=ny)

    ra, dec = _ccd_corners(w, (nx, ny))

    assert len(ra) == 4
    assert len(dec) == 4
    # The footprint straddles the CRVAL centre in both coordinates.
    assert min(ra) < 100.0 < max(ra)
    assert min(dec) < 20.0 < max(dec)
    # Extent matches the CD-scaled pixel span (with edge padding).
    npt.assert_allclose(max(dec) - min(dec), ny * 1e-5, rtol=1e-3)


def test_parse_header_to_wcs_returns_one_wcs_per_extension(tmp_path):
    """One (wcs, shape) pair is returned per CCD HDU; the primary is skipped."""
    ccd_wcs = [_tan_wcs(10.0), _tan_wcs(20.0), _tan_wcs(30.0)]

    path = tmp_path / "1234567.txt"
    path.write_text(_header_text(ccd_wcs))

    result = _parse_header_to_wcs(str(path))

    assert len(result) == len(ccd_wcs)
    assert all(shape == (2080, 4612) for _, shape in result)


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


# ---------------------------------------------------------------------------
# FieldCornersExtractor resume path (per-CCD done-set)
# ---------------------------------------------------------------------------

def _write_headers(header_dir, expnums, n_ccd=3):
    """Write one plain multi-HDU header per exposure into ``header_dir``."""
    header_dir.mkdir(exist_ok=True)
    for j, expnum in enumerate(expnums):
        ccd_wcs = [_tan_wcs(10.0 + j + i) for i in range(n_ccd)]
        (header_dir / f"{expnum}.txt").write_text(_header_text(ccd_wcs))


def test_get_done_ccds_reads_present_ids(tmp_path):
    """The done-set is the exact set of CCD IDs already in the output."""
    out = tmp_path / "corners.txt"
    out.write_text(
        "1234567-0 10 10 10 10 20 20 20 20\n"
        "1234567-2 10 10 10 10 20 20 20 20\n"
    )

    extractor = FieldCornersExtractor()
    extractor._params["output_file"] = str(out)

    assert extractor.get_done_ccds() == {"1234567-0", "1234567-2"}


def test_resume_adds_new_exposure_without_duplicating(tmp_path):
    """Resume appends a new exposure and leaves existing rows untouched."""
    header_dir = tmp_path / "headers"
    _write_headers(header_dir, [1000001, 1000002])

    out = tmp_path / "corners.txt"
    # Pre-populate with exposure 1's three CCD rows.
    pre = FieldCornersExtractor().process_single_header(
        (str(header_dir / "1000001.txt"), False, None)
    )
    with open(out, "w") as f:
        for ccd_id, ra, dec in pre:
            f.write(f"{ccd_id} " + " ".join(f"{v:.5f}" for v in ra + dec) + "\n")

    FieldCornersExtractor().run(
        args=["-i", str(header_dir), "-o", str(out), "-r"]
    )

    ids = [line.split()[0] for line in out.read_text().splitlines()]
    # No duplicate exposure-1 rows; exposure 2's three CCDs added.
    assert ids.count("1000001-0") == 1
    assert sorted(ids) == [
        "1000001-0", "1000001-1", "1000001-2",
        "1000002-0", "1000002-1", "1000002-2",
    ]


def test_resume_completes_partial_exposure(tmp_path):
    """An exposure interrupted mid-write is completed, not skipped."""
    header_dir = tmp_path / "headers"
    _write_headers(header_dir, [1000001])

    out = tmp_path / "corners.txt"
    # Simulate an interrupt: only the first of exposure 1's CCDs was written.
    pre = FieldCornersExtractor().process_single_header(
        (str(header_dir / "1000001.txt"), False, None)
    )
    ccd_id, ra, dec = pre[0]
    with open(out, "w") as f:
        f.write(f"{ccd_id} " + " ".join(f"{v:.5f}" for v in ra + dec) + "\n")

    FieldCornersExtractor().run(
        args=["-i", str(header_dir), "-o", str(out), "-r"]
    )

    ids = [line.split()[0] for line in out.read_text().splitlines()]
    # The missing CCDs are filled in; the present one is not duplicated.
    assert sorted(ids) == ["1000001-0", "1000001-1", "1000001-2"]
    assert ids.count("1000001-0") == 1


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
# CoverageMapBuilder — parsing, nexp contract, RA-wrap and pole guards
# ---------------------------------------------------------------------------

def _read_map(path):
    import healsparse as hsp

    return hsp.HealSparseMap.read(str(path))


def test_build_map_single_row(tmp_path):
    """A one-row corners file exercises the atleast_1d/2d parsing guards."""
    infile = tmp_path / "corners.txt"
    infile.write_text("100-0 9.9 10.1 10.1 9.9 0.0 0.0 0.2 0.2\n")
    out = tmp_path / "cov.hsp"

    CoverageMapBuilder().run(
        args=["-i", str(infile), "-o", str(out), "-c", "32", "-n", "1024"]
    )

    m = _read_map(out)
    assert 0 < len(m.valid_pixels) < 200
    assert m[m.valid_pixels].max() == 1


def test_build_map_nexp_counts_overlapping_exposures(tmp_path):
    """Two overlapping CCDs from different exposures give value 2 in overlap."""
    # Two 0.4x0.4 deg CCDs offset by 0.2 deg in RA -> a central overlap strip.
    infile = tmp_path / "corners.txt"
    infile.write_text(
        "100-0  9.8 10.2 10.2  9.8 19.8 19.8 20.2 20.2\n"
        "200-0 10.0 10.4 10.4 10.0 19.8 19.8 20.2 20.2\n"
    )
    out = tmp_path / "cov.hsp"

    CoverageMapBuilder().run(
        args=["-i", str(infile), "-o", str(out), "-c", "32", "-n", "1024"]
    )

    m = _read_map(out)
    values = m[m.valid_pixels]
    # The overlap is covered by both exposures (value 2); the union edges by
    # one (value 1). Both must be present; nothing exceeds 2.
    assert values.max() == 2
    assert values.min() == 1
    assert (values == 2).sum() > 0
    assert (values == 1).sum() > 0


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


def test_build_map_invariant_under_ra_shift(tmp_path):
    """A seam CCD's footprint matches an identical CCD shifted +10 deg in RA.

    Rotating both the seam CCD (via +360/unwrap) and a reference CCD onto the
    same RA and comparing pixel counts fails if the seam polygon were filling
    the ~360 deg complement. This is the real RA-wrap regression guard.
    """
    # Seam CCD straddling RA=0, and the same CCD translated to RA~10.
    seam = tmp_path / "seam.txt"
    seam.write_text("100-0 359.9 0.1 0.1 359.9 20.0 20.0 20.2 20.2\n")
    ref = tmp_path / "ref.txt"
    ref.write_text("200-0 9.9 10.1 10.1 9.9 20.0 20.0 20.2 20.2\n")

    m_seam = tmp_path / "seam.hsp"
    m_ref = tmp_path / "ref.hsp"
    CoverageMapBuilder().run(
        args=["-i", str(seam), "-o", str(m_seam), "-c", "32", "-n", "1024"]
    )
    CoverageMapBuilder().run(
        args=["-i", str(ref), "-o", str(m_ref), "-c", "32", "-n", "1024"]
    )

    n_seam = len(_read_map(m_seam).valid_pixels)
    n_ref = len(_read_map(m_ref).valid_pixels)

    # Same-size footprints at the same declination: pixel counts agree to
    # within a few boundary pixels, and are nowhere near a hemisphere.
    assert abs(n_seam - n_ref) <= max(3, int(0.1 * n_ref))
    assert 0 < n_seam < 200


def test_build_map_pole_guard_skips_polygon(tmp_path, capsys):
    """A polygon with |dec| near 90 deg is skipped with a warning."""
    infile = tmp_path / "corners.txt"
    infile.write_text(
        "100-0  10.0  10.2  10.2  10.0 89.5 89.5 89.7 89.7\n"
        "200-0  10.0  10.2  10.2  10.0 20.0 20.0 20.2 20.2\n"
    )
    out = tmp_path / "cov.hsp"

    CoverageMapBuilder().run(
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
