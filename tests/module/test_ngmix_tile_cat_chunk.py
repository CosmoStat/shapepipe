"""A chunk's Tile_cat holds the stamps of its own rows, and only those.

Each ngmix chunk reads the tile catalogue for rows ``chunk_rows(n_obj,
row_min, row_max)``. The per-object columns stay full length; the stamp
columns (``VIGNET`` and ``SEG_VIGNET``) are held for the chunk's rows only,
indexed by tile-catalogue row exactly like the full column.
"""

import numpy as np
import pytest
from astropy.io import fits

from shapepipe.modules.ngmix_package.ngmix import (
    ChunkStamps,
    Tile_cat,
    chunk_rows,
)

N_OBJ = 23
STAMP = 7


def _ldac(path, columns):
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="1A", array=["x"])],
        name="LDAC_IMHEAD",
    )
    objects = fits.BinTableHDU.from_columns(columns, name="LDAC_OBJECTS")
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(path)
    return str(path)


@pytest.fixture
def catalogue(tmp_path):
    """A tile catalogue with VIGNET and SEG_VIGNET stamp columns."""
    rng = np.random.default_rng(3)
    number = rng.permutation(np.arange(1, 10 * N_OBJ, 10))[:N_OBJ]
    vign = rng.normal(size=(N_OBJ, STAMP, STAMP)).astype(np.float32)
    vign[::4, 0, :] = -1e30
    vign[1, 2, 2] = np.nan
    seg = rng.integers(0, 5, size=(N_OBJ, STAMP, STAMP)).astype(np.int32)
    dim = f"({STAMP}, {STAMP})"
    cat = _ldac(tmp_path / "sexcat.fits", [
        fits.Column(name="NUMBER", format="J", array=number),
        fits.Column(name="XWIN_WORLD", format="D", array=rng.uniform(size=N_OBJ)),
        fits.Column(name="YWIN_WORLD", format="D", array=rng.uniform(size=N_OBJ)),
        fits.Column(name="FLUX_AUTO", format="E", array=rng.uniform(size=N_OBJ)),
        fits.Column(
            name="VIGNET", format=f"{STAMP * STAMP}E", dim=dim, array=vign
        ),
        fits.Column(
            name="SEG_VIGNET", format=f"{STAMP * STAMP}J", dim=dim, array=seg
        ),
    ])
    return cat


def _full_columns(path):
    """Every column as an in-memory (not memory-mapped) read gives it."""
    with fits.open(path, memmap=False) as hdul:
        data = hdul[2].data
        return {name: np.copy(data[name]) for name in data.dtype.names}


@pytest.mark.parametrize(
    "row_min, row_max",
    [(-1, -1), (1, 6), (7, 15), (16, N_OBJ), (20, -1), (N_OBJ + 1, N_OBJ)],
)
def test_chunk_holds_exactly_its_rows(catalogue, row_min, row_max):
    full = _full_columns(catalogue)
    rows = chunk_rows(N_OBJ, row_min, row_max)

    tile = Tile_cat(catalogue, row_min=row_min, row_max=row_max)

    assert tile.rows == rows
    # Per-object columns: full length, every row.
    np.testing.assert_array_equal(tile.obj_id, full["NUMBER"])
    np.testing.assert_array_equal(tile.ra, full["XWIN_WORLD"])
    np.testing.assert_array_equal(tile.dec, full["YWIN_WORLD"])
    np.testing.assert_array_equal(tile.flux, full["FLUX_AUTO"])

    for stamps, column in ((tile.vign, full["VIGNET"]),
                           (tile.seg, full["SEG_VIGNET"])):
        # Same stamps, bit for bit, at the same tile-catalogue rows.
        for i_tile in rows:
            assert stamps[i_tile].dtype == column[i_tile].dtype
            assert stamps[i_tile].tobytes() == column[i_tile].tobytes()
        # Only the chunk's stamps are held, in memory of their own.
        assert stamps._stamps.shape == (len(rows), STAMP, STAMP)
        assert stamps._stamps.flags.owndata
        # The other rows are absent, and touching one says so.
        for i_tile in set(range(N_OBJ)) - set(rows):
            with pytest.raises(IndexError, match="outside this chunk"):
                stamps[i_tile]


def test_stamps_outlive_the_catalogue_file(catalogue):
    """The held stamps are copies, not views into the memory-mapped file."""
    full = _full_columns(catalogue)
    tile = Tile_cat(catalogue, row_min=5, row_max=9)
    # Overwrite the file in place with different stamps.
    with fits.open(catalogue, mode="update") as hdul:
        hdul[2].data["VIGNET"][:] = 0
        hdul[2].data["SEG_VIGNET"][:] = 0
    assert tile.vign._stamps.tobytes() == full["VIGNET"][4:9].tobytes()
    assert tile.seg._stamps.tobytes() == full["SEG_VIGNET"][4:9].tobytes()


def test_chunk_stamps_index_by_tile_row():
    column = np.arange(10 * 4).reshape(10, 2, 2)
    stamps = ChunkStamps(column, range(3, 6))
    np.testing.assert_array_equal(stamps[np.int64(3)], column[3])
    np.testing.assert_array_equal(stamps[5], column[5])
    for i_tile in (2, 6, -1):
        with pytest.raises(IndexError):
            stamps[i_tile]
