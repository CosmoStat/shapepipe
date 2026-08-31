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
    def info(self, *_args, **_kwargs):
        pass


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


def _write_sexcat(path):
    """Write a synthetic LDAC SExtractor catalogue of known positions.

    Written with astropy rather than ``FITSCatalogue.save_as_fits``, because
    creating an LDAC catalogue through that API requires an existing LDAC file
    to copy the ``LDAC_IMHEAD`` HDU from. The three-HDU layout below is what
    SExtractor writes and what ``SEx_catalogue=True`` (``hdu_no=2``) indexes.
    """
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
            fits.Column(name="NUMBER", format="J", array=np.arange(len(RA))),
            fits.Column(name="XWIN_WORLD", format="D", array=RA),
            fits.Column(name="YWIN_WORLD", format="D", array=DEC),
            fits.Column(
                name="IMAFLAGS_ISO",
                format="J",
                array=np.zeros(len(RA), dtype="i4"),
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
