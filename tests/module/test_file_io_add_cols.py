"""UNIT TESTS FOR PIPELINE: file_io.FITSCatalogue.add_cols.

``add_cols`` appends many columns with one rebuild and one write of the
file. Its contract is that the file it writes is the file the same columns
added one ``add_col`` call at a time would produce -- which writes and
reopens the file between every column. These tests build both from the same
base catalogue and compare them HDU by HDU: names, column formats and
``TDIM``, headers, data, and finally the bytes.
"""

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits

from shapepipe.pipeline import file_io

N_OBJ = 7


def _new_columns():
    """Columns of every kind make_cat appends, in a fixed order."""
    rng = np.random.default_rng(1)
    return {
        "TILE_ID": np.full(N_OBJ, 210.282),
        "TILE_UNIQUE_ID": np.arange(N_OBJ, dtype=np.int64) + 10**9,
        "FLAGS_I16": rng.integers(-5, 5, N_OBJ).astype(np.int16),
        "N_I32": rng.integers(0, 100, N_OBJ).astype(np.int32),
        "G1_F32": rng.normal(size=N_OBJ).astype(np.float32),
        "MASK_BOOL": rng.integers(0, 2, N_OBJ).astype(bool),
        # A per-epoch slot block as a vector column, and a 2-D cell.
        "PSF_SLOTS": rng.normal(size=(N_OBJ, 12)),
        "COV": rng.normal(size=(N_OBJ, 2, 3)),
        "EXP_NAME": np.array([f"2113864-{i}" for i in range(N_OBJ)]),
    }


def _write_base(path, sex_catalogue):
    """A base catalogue; with ``sex_catalogue``, the table is HDU 2 of 4."""
    base = np.empty(N_OBJ, dtype=[("NUMBER", "i4"), ("X", "f8")])
    base["NUMBER"] = np.arange(1, N_OBJ + 1)
    base["X"] = np.linspace(0.0, 1.0, N_OBJ)
    if not sex_catalogue:
        cat = file_io.FITSCatalogue(
            str(path),
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
        )
        cat.save_as_fits(base, ext_name="RESULTS")
        return

    header = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="10A", array=["x"])],
        name="LDAC_IMHEAD",
    )
    table = fits.BinTableHDU(base, name="LDAC_OBJECTS")
    epoch = fits.BinTableHDU(
        np.array([(1, 3)], dtype=[("NUMBER", "i4"), ("CCD_N", "i4")]),
        name="EPOCH_0",
    )
    fits.HDUList([fits.PrimaryHDU(), header, table, epoch]).writeto(path)


def _open(path, sex_catalogue):
    cat = file_io.FITSCatalogue(
        str(path),
        open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
        SEx_catalogue=sex_catalogue,
    )
    cat.open()
    return cat


@pytest.mark.parametrize("sex_catalogue", [False, True])
def test_add_cols_matches_successive_add_col(tmp_path, sex_catalogue):
    """One add_cols call writes the file N add_col calls would."""
    columns = _new_columns()
    one_path = tmp_path / "one_by_one.fits"
    batch_path = tmp_path / "batched.fits"

    _write_base(one_path, sex_catalogue)
    cat = _open(one_path, sex_catalogue)
    for name, data in columns.items():
        cat.add_col(name, data)
    cat.close()

    _write_base(batch_path, sex_catalogue)
    cat = _open(batch_path, sex_catalogue)
    cat.add_cols(columns)
    cat.close()

    with fits.open(one_path) as ref, fits.open(batch_path) as new:
        assert [h.name for h in new] == [h.name for h in ref]
        for h_ref, h_new in zip(ref, new):
            assert h_new.header == h_ref.header, h_ref.name
            if not isinstance(h_ref, fits.BinTableHDU):
                continue
            assert h_new.columns.names == h_ref.columns.names
            assert h_new.columns.formats == h_ref.columns.formats
            assert h_new.columns.dims == h_ref.columns.dims
            for name in h_ref.columns.names:
                npt.assert_array_equal(
                    h_new.data[name], h_ref.data[name], err_msg=name
                )

        table = new[2 if sex_catalogue else 1]
        assert table.columns.names[-len(columns):] == list(columns)
        assert table.data["PSF_SLOTS"].shape == (N_OBJ, 12)
        npt.assert_array_equal(table.data["FLAGS_I16"], columns["FLAGS_I16"])

    assert batch_path.read_bytes() == one_path.read_bytes()


def test_add_cols_leaves_catalogue_open_and_readable(tmp_path):
    """After add_cols the instance is reopened on the written file."""
    path = tmp_path / "cat.fits"
    _write_base(path, sex_catalogue=False)
    cat = _open(path, sex_catalogue=False)
    cat.add_cols({"A": np.arange(N_OBJ), "B": np.ones(N_OBJ)})
    assert cat.get_col_names() == ["NUMBER", "X", "A", "B"]
    npt.assert_array_equal(cat.get_data()["A"], np.arange(N_OBJ))
    cat.close()


def test_add_cols_rejects_non_array_before_writing(tmp_path):
    """A non-array column raises and leaves the file as it was."""
    path = tmp_path / "cat.fits"
    _write_base(path, sex_catalogue=False)
    before = path.read_bytes()
    cat = _open(path, sex_catalogue=False)
    with pytest.raises(TypeError):
        cat.add_cols({"A": np.arange(N_OBJ), "B": list(range(N_OBJ))})
    cat.close()
    assert path.read_bytes() == before
