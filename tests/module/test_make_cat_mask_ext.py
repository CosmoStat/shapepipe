"""UNIT TESTS FOR MODULE PACKAGE: MAKE_CAT external-mask columns.

Exercises the optional per-band external healsparse mask lookup added to
``make_cat`` (PR #847 §4, the ShapePipe end of UNIONS-WL/spherex#38). A small
synthetic final catalogue carrying known ``XWIN_WORLD`` / ``YWIN_WORLD``
object positions is queried against synthetic healsparse maps of known value,
locking in: (1) the ``MASK_<BAND>`` column name and per-object values, (2) the
off-map sentinel (``-1`` for integer maps) written verbatim for objects outside
coverage, and (3) multi-band handling.
"""

import numpy as np
import numpy.testing as npt
import pytest

healsparse = pytest.importorskip("healsparse")

from shapepipe.modules.make_cat_package import make_cat


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


NSIDE_COVERAGE = 32
NSIDE_SPARSE = 4096

# Object world positions (RA, Dec in degrees). The last object sits far from
# the map coverage so it exercises the off-map sentinel path.
RA = np.array([10.0, 10.1, 10.2, 200.0])
DEC = np.array([20.0, 20.1, 20.2, -40.0])


def _make_map(value, dtype=np.int16, sentinel=-1):
    """Build a healsparse map covering the first three RA/Dec positions.

    All covered pixels carry ``value``; everything else reads the sentinel.
    """
    smap = healsparse.HealSparseMap.make_empty(
        NSIDE_COVERAGE, NSIDE_SPARSE, dtype, sentinel=sentinel
    )
    smap.update_values_pos(
        RA[:3], DEC[:3], np.full(3, value, dtype=dtype), lonlat=True
    )
    return smap


def _final_cat():
    """Synthetic final-catalogue columns with known positions."""
    return {
        "NUMBER": np.arange(len(RA)),
        "XWIN_WORLD": RA,
        "YWIN_WORLD": DEC,
    }


def test_parse_mask_ext_paths():
    """band:path pairs parse into a stripped mapping, whitespace-tolerant."""
    parsed = make_cat.parse_mask_ext_paths(
        "u:/a/mask_u.hsp, g:/b/mask_g.hsp,r:/c/mask_r.hsp"
    )
    assert parsed == {
        "u": "/a/mask_u.hsp",
        "g": "/b/mask_g.hsp",
        "r": "/c/mask_r.hsp",
    }


def test_mask_ext_columns(tmp_path):
    """Per-band columns carry the map value on-map and the sentinel off-map."""
    u_map = _make_map(16)
    g_map = _make_map(32)
    u_path = tmp_path / "mask_u.hsp"
    g_path = tmp_path / "mask_g.hsp"
    u_map.write(str(u_path))
    g_map.write(str(g_path))

    cat = _final_cat()
    make_cat.save_mask_ext_data(
        cat,
        {"u": str(u_path), "g": str(g_path)},
        _NullLogger(),
    )

    # On-map objects (first three) carry the map value; the off-map object
    # (last) carries the map's -1 sentinel.
    npt.assert_array_equal(cat["MASK_u"], [16, 16, 16, -1])
    npt.assert_array_equal(cat["MASK_g"], [32, 32, 32, -1])
    # Integer dtype preserved from the map.
    assert np.issubdtype(cat["MASK_u"].dtype, np.integer)
    assert list(cat) == [
        "NUMBER", "XWIN_WORLD", "YWIN_WORLD", "MASK_u", "MASK_g"
    ]
