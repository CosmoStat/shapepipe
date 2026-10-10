"""The per-tile final-catalogue reader reads HDF5 and FITS alike.

``shapepipe.utilities.final_cat.read_final_cat`` is what every final-catalogue
helper reads through (``create_final_cat.py``, the workflow's
``merge_final_cat.py`` by way of it, and ``scripts/python/merge_final_cat.py``).
One catalogue written both ways — HDF5 by make_cat's own writer, FITS as a
binary table at HDU 1 — must read to the same structured array, with the format
decided by content rather than file name.
"""

import shutil

import numpy as np
import pytest
from astropy.io import fits

from shapepipe.modules.make_cat_package.make_cat import write_final_cat
from shapepipe.utilities.final_cat import is_hdf5, read_final_cat


def _catalogue(rows=25, seed=3):
    """A catalogue with scalar, vector, integer and boolean columns."""
    rng = np.random.default_rng(seed)
    cat = np.zeros(rows, dtype=[
        ("NUMBER", "i4"),
        ("XWIN_WORLD", "f8"),
        ("MAG_AUTO", "f4"),
        ("FLUX_APER", "f4", (3,)),
        ("TILE_UNIQUE_ID", "i8"),
        ("NGMIX_FLAGS", "bool"),
    ])
    cat["NUMBER"] = np.arange(1, rows + 1)
    cat["XWIN_WORLD"] = rng.uniform(0, 360, rows)
    cat["MAG_AUTO"] = rng.uniform(18, 25, rows)
    cat["FLUX_APER"] = rng.uniform(0, 100, (rows, 3))
    cat["TILE_UNIQUE_ID"] = rng.integers(0, 2**40, rows)
    cat["NGMIX_FLAGS"] = rng.integers(0, 2, rows).astype(bool)
    return cat


@pytest.fixture
def both(tmp_path):
    cat = _catalogue()
    h5 = tmp_path / "final_cat-210-282.hdf5"
    write_final_cat(str(h5), {n: cat[n] for n in cat.dtype.names})
    fts = tmp_path / "final_cat-210-282.fits"
    fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU(cat, name="RESULTS")]
    ).writeto(fts)
    return cat, h5, fts


def _assert_same(got, want):
    assert got.dtype.names == want.dtype.names
    for name in want.dtype.names:
        assert got.dtype[name] == want.dtype[name], name
        assert np.array_equal(got[name], want[name]), name


def test_every_column_reads_alike_from_either_format(both):
    cat, h5, fts = both
    from_h5, from_fits = read_final_cat(str(h5)), read_final_cat(str(fts))
    _assert_same(from_h5, cat)
    _assert_same(from_fits, cat)
    assert from_fits.dtype["FLUX_APER"].shape == (3,)
    assert all(from_fits.dtype[n].base.isnative for n in cat.dtype.names)


def test_requested_columns_come_back_in_the_requested_order(both):
    cat, h5, fts = both
    wanted = ["FLUX_APER", "NUMBER", "NGMIX_FLAGS"]
    for path in (h5, fts):
        got = read_final_cat(str(path), wanted)
        assert got.dtype.names == tuple(wanted), path
        _assert_same(got, cat[wanted])


def test_format_is_decided_by_content_not_name(both, tmp_path):
    cat, h5, fts = both
    h5_named_fits = tmp_path / "hdf5_content.fits"
    fits_named_hdf5 = tmp_path / "fits_content.hdf5"
    shutil.copy(h5, h5_named_fits)
    shutil.copy(fts, fits_named_hdf5)
    assert is_hdf5(str(h5_named_fits)) and not is_hdf5(str(fits_named_hdf5))
    _assert_same(read_final_cat(str(h5_named_fits)), cat)
    _assert_same(read_final_cat(str(fits_named_hdf5)), cat)


@pytest.mark.parametrize("fmt", ["hdf5", "fits"])
def test_a_missing_column_raises_naming_it(both, fmt):
    _, h5, fts = both
    path = h5 if fmt == "hdf5" else fts
    with pytest.raises(KeyError, match="NOT_A_COLUMN") as err:
        read_final_cat(str(path), ["NUMBER", "NOT_A_COLUMN"])
    assert str(path) in str(err.value)
