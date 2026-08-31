"""UNIT TESTS FOR MODULE PACKAGE: MASK_QUERY.

Exercises the exposure-side half of the query-everything mask design (PR #847):
``mask_query`` reads each SExtractor detection's windowed world position out of
the ``LDAC_OBJECTS`` extension, looks it up in the configured healsparse maps,
and writes a single integer ``FLAG_EXT`` column into a NEW catalogue beside the
input.

What is locked in here: (1) boolean maps flag with ``1``, (2) integer maps
contribute their value and ``MASK_BITS`` restricts which bits do, (3) the
off-coverage sentinel (``-1``) never flags — the one place this differs from
``make_cat``'s verbatim pass-through, argued in
:mod:`shapepipe.utilities.mask_query` — (4) several maps OR together, and
(5) the input catalogue is left untouched while the LDAC structure survives.
"""

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits

healsparse = pytest.importorskip("healsparse")

from shapepipe.modules.mask_query_package.mask_query import MaskQuery
from shapepipe.pipeline import file_io
from shapepipe.utilities import mask_query as mask_query_util

NSIDE_COVERAGE = 32
NSIDE_SPARSE = 4096

# Detection world positions (RA, Dec in degrees). The last one sits far from
# every map's coverage, so it exercises the off-coverage path.
RA = np.array([10.0, 10.1, 10.2, 200.0])
DEC = np.array([20.0, 20.1, 20.2, -40.0])


class _NullLogger:
    """Captures what the module logs, so coverage warnings can be asserted."""

    def __init__(self):
        self.info_msgs = []
        self.warning_msgs = []

    def info(self, msg, *_a, **_kw):
        self.info_msgs.append(str(msg))

    def warning(self, msg, *_a, **_kw):
        self.warning_msgs.append(str(msg))


def _write_map(path, value, dtype=np.int16, n_covered=2):
    """Build an integer map carrying ``value`` at the first ``n_covered``.

    Everything else reads the ``-1`` sentinel, i.e. off coverage.
    """
    smap = healsparse.HealSparseMap.make_empty(
        NSIDE_COVERAGE, NSIDE_SPARSE, dtype, sentinel=-1
    )
    if n_covered:
        smap.update_values_pos(
            RA[:n_covered],
            DEC[:n_covered],
            np.full(n_covered, value, dtype=dtype),
            lonlat=True,
        )
    smap.write(str(path))
    return str(path)


def _write_bool_map(path, n_covered=2):
    """Build a boolean map, ``True`` at the first ``n_covered``."""
    smap = healsparse.HealSparseMap.make_empty(
        NSIDE_COVERAGE, NSIDE_SPARSE, np.bool_
    )
    smap.update_values_pos(
        RA[:n_covered],
        DEC[:n_covered],
        np.ones(n_covered, dtype=np.bool_),
        lonlat=True,
    )
    smap.write(str(path))
    return str(path)


def _write_sexcat(path, n=None):
    """Write a synthetic LDAC SExtractor catalogue of known positions.

    Written with astropy rather than ``FITSCatalogue.save_as_fits``, because
    creating an LDAC catalogue through that API requires an existing LDAC file
    to copy the ``LDAC_IMHEAD`` HDU from. The three-HDU layout below is what
    SExtractor writes and what ``SEx_catalogue=True`` (``hdu_no=2``) indexes.
    """
    n = len(RA) if n is None else n
    imhead = fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="Field Header Card", format="1A", array=np.array(["x"])
            )
        ],
        name="LDAC_IMHEAD",
    )
    objects = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="NUMBER", format="J", array=np.arange(n)),
            fits.Column(name="XWIN_WORLD", format="D", array=RA[:n]),
            fits.Column(name="YWIN_WORLD", format="D", array=DEC[:n]),
            fits.Column(
                name="IMAFLAGS_ISO", format="J", array=np.zeros(n, dtype="i4")
            ),
        ],
        name="LDAC_OBJECTS",
    )
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(
        str(path), overwrite=True
    )
    return str(path)


def _read(path):
    cat = file_io.FITSCatalogue(str(path), SEx_catalogue=True)
    cat.open()
    data = cat.get_data()
    flag = np.copy(data["FLAG_EXT"])
    names = set(data.dtype.names)
    cat.close()
    return flag, names


def test_parse_map_paths():
    """Comma-separated paths parse, whitespace-tolerant, empties dropped."""
    assert mask_query_util.parse_map_paths(
        " /a/star.hsp, /b/maximask.hsp ,, "
    ) == ["/a/star.hsp", "/b/maximask.hsp"]


def test_flag_positions_integer_map(tmp_path):
    """An integer map contributes its value; off-coverage stays clean."""
    path = _write_map(tmp_path / "m.hsp", 4, n_covered=3)
    npt.assert_array_equal(
        mask_query_util.flag_positions([path], RA, DEC), [4, 4, 4, 0]
    )


def test_flag_positions_bits_restrict(tmp_path):
    """MASK_BITS selects which bits of an integer map flag."""
    path = _write_map(tmp_path / "m.hsp", 1028, n_covered=3)
    npt.assert_array_equal(
        mask_query_util.flag_positions([path], RA, DEC, bits=4), [4, 4, 4, 0]
    )
    # A bit the map does not carry leaves everything clean.
    npt.assert_array_equal(
        mask_query_util.flag_positions([path], RA, DEC, bits=2), [0, 0, 0, 0]
    )


def test_flag_positions_ors_maps(tmp_path):
    """Several maps combine with a bitwise OR."""
    a = _write_map(tmp_path / "a.hsp", 4, n_covered=1)
    b = _write_map(tmp_path / "b.hsp", 1024, n_covered=3)
    npt.assert_array_equal(
        mask_query_util.flag_positions([a, b], RA, DEC), [1028, 1024, 1024, 0]
    )


def test_flag_positions_boolean_map(tmp_path):
    """A boolean map flags with 1; its False sentinel stays clean."""
    path = _write_bool_map(tmp_path / "bool.hsp", n_covered=2)
    npt.assert_array_equal(
        mask_query_util.flag_positions([path], RA, DEC), [1, 1, 0, 0]
    )


def test_mask_query_writes_flag_ext(tmp_path):
    """The module writes FLAG_EXT into a new catalogue and counts the hits."""
    map_path = _write_map(tmp_path / "star.hsp", 4, n_covered=2)
    in_path = _write_sexcat(tmp_path / "sexcat-000-0.fits")
    out_path = tmp_path / "sexcat_ext-000-0.fits"

    n_flagged = MaskQuery(
        in_path, str(out_path), [map_path], w_log=_NullLogger()
    ).process()

    assert n_flagged == 2
    flag, names = _read(out_path)
    npt.assert_array_equal(flag, [4, 4, 0, 0])
    assert np.issubdtype(flag.dtype, np.integer)
    # The columns SExtractor wrote survive alongside the new one.
    assert {"NUMBER", "XWIN_WORLD", "YWIN_WORLD", "IMAFLAGS_ISO"} <= names
    # The LDAC structure survives: setools and psfex read HDU 2 by index.
    with fits.open(str(out_path)) as hdus:
        assert [hdu.name for hdu in hdus] == [
            "PRIMARY",
            "LDAC_IMHEAD",
            "LDAC_OBJECTS",
        ]

    # The input is not mutated: this module publishes a new file.
    with fits.open(in_path) as hdus:
        assert "FLAG_EXT" not in hdus[2].data.dtype.names


def test_mask_query_bits_and_all_clean(tmp_path):
    """MASK_BITS reaches the module, and a miss leaves every object clean."""
    map_path = _write_map(tmp_path / "m.hsp", 1024, n_covered=3)
    in_path = _write_sexcat(tmp_path / "sexcat-000-1.fits")
    out_path = tmp_path / "sexcat_ext-000-1.fits"

    n_flagged = MaskQuery(
        in_path, str(out_path), [map_path], bits=4, w_log=_NullLogger()
    ).process()

    assert n_flagged == 0
    flag, _ = _read(out_path)
    npt.assert_array_equal(flag, [0, 0, 0, 0])


def test_partial_read_matches_full(tmp_path):
    """The partial read returns exactly what a full read of the map returns.

    This is the guarantee that makes the optimisation safe: query_map loads
    only the coverage pixels the positions touch, and must be indistinguishable
    from HealSparseMap.read(path) at every queried position.
    """
    for dtype, writer in ((np.int16, _write_map), (np.bool_, None)):
        path = (
            _write_map(tmp_path / f"full_{dtype.__name__}.hsp", 7, n_covered=3)
            if writer is not None
            else _write_bool_map(tmp_path / "full_bool.hsp", n_covered=3)
        )
        full = healsparse.HealSparseMap.read(path)
        expected = np.asarray(full.get_values_pos(RA, DEC, lonlat=True))
        npt.assert_array_equal(
            mask_query_util.query_map(path, RA, DEC), expected
        )


def test_coverage_reported_for_both_map_kinds(tmp_path):
    """in_coverage is the coverage mask, not valid_mask.

    For a boolean map healsparse stores only the True pixels, so valid_mask
    would equal the value and could not distinguish an unmasked object from one
    the map does not reach. The distant object must read as off-coverage while
    the near, unflagged ones read as covered.
    """
    for path in (
        _write_map(tmp_path / "int.hsp", 4, n_covered=2),
        _write_bool_map(tmp_path / "bool.hsp", n_covered=2),
    ):
        _, in_coverage = mask_query_util.query_map_coverage(path, RA, DEC)
        # The 4th position is 190 deg away; the first two are on the map.
        assert in_coverage[0] and in_coverage[1]
        assert not in_coverage[3]


def test_all_off_coverage_warns(tmp_path):
    """A map that reaches nothing warns, instead of logging a silent zero."""
    path = _write_map(tmp_path / "elsewhere.hsp", 4, n_covered=2)
    far_ra = np.array([200.0, 201.0])
    far_dec = np.array([-40.0, -41.0])
    log = _NullLogger()

    flag = mask_query_util.flag_positions([path], far_ra, far_dec, w_log=log)

    npt.assert_array_equal(flag, [0, 0])
    assert any("NO object is inside" in m for m in log.warning_msgs)
    assert any("2 outside coverage" in m for m in log.info_msgs)
    # A map that DOES reach the objects must not warn.
    log2 = _NullLogger()
    mask_query_util.flag_positions([path], RA, DEC, w_log=log2)
    assert not log2.warning_msgs


def test_flag_positions_empty_input(tmp_path):
    """Zero positions is not an error and reads no map."""
    flag = mask_query_util.flag_positions(
        [str(tmp_path / "does-not-exist.hsp")], np.zeros(0), np.zeros(0)
    )
    assert flag.shape == (0,)


def test_mask_query_empty_ccd(tmp_path):
    """A CCD with no detections still publishes a catalogue, not an error.

    setools tolerates sparse-CCD attrition and psfex_interp's completeness
    floor is 0/warn, so an empty sexcat must flow through rather than raise.
    """
    map_path = _write_map(tmp_path / "star.hsp", 4, n_covered=2)
    in_path = _write_sexcat(tmp_path / "sexcat-000-2.fits", n=0)
    out_path = tmp_path / "sexcat_ext-000-2.fits"
    log = _NullLogger()

    n_flagged = MaskQuery(
        in_path, str(out_path), [map_path], w_log=log
    ).process()

    assert n_flagged == 0
    assert out_path.exists()
    flag, names = _read(out_path)
    assert flag.shape == (0,)
    assert "FLAG_EXT" in names
    assert any("No detections" in m for m in log.info_msgs)


def test_empty_coverage_raises_clearly(tmp_path):
    """A map with no coverage at all fails with a message naming the file.

    The probe read in query_map_coverage indexes the first covered pixel; with
    nothing covered that would be an IndexError from deep inside the utility.
    """
    path = _write_map(tmp_path / "empty.hsp", 4, n_covered=0)
    with pytest.raises(ValueError):
        mask_query_util.query_map_coverage(path, RA, DEC)
