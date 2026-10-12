"""ngmix chunks select catalogue rows, so they cover every object once.

The workflow partitions a tile into chunks with ``ngmix_range.row_ranges``
(1-based closed row ranges) and each ngmix chunk keeps the rows
``ngmix.chunk_rows`` returns. The pair must visit every row exactly once
whatever the ``NUMBER`` column holds: an external detection catalogue keeps
its own, possibly gapped and unsorted, ``NUMBER``.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from hypothesis import given, settings
from hypothesis import strategies as st

from shapepipe.modules.ngmix_package.ngmix import chunk_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "ngmix_range.py"


def _load_ngmix_range():
    spec = importlib.util.spec_from_file_location("_ngmix_range", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ngmix_range = _load_ngmix_range()


@settings(max_examples=200, deadline=None)
@given(
    st.integers(min_value=1, max_value=300),
    st.integers(min_value=1, max_value=12),
    st.data(),
)
def test_chunks_cover_every_row_once(n_obj, n_chunks, data):
    epochs = data.draw(
        st.lists(st.integers(0, 20), min_size=n_obj, max_size=n_obj)
    )
    # Gapped, unsorted NUMBER: distinct values drawn from a wide range.
    number = np.array(
        data.draw(
            st.lists(
                st.integers(1, 10**6 - 1),
                min_size=n_obj,
                max_size=n_obj,
                unique=True,
            )
        )
    )

    selected = [
        number[i]
        for lo, hi in ngmix_range.row_ranges(epochs, n_chunks)
        for i in chunk_rows(n_obj, lo, hi)
    ]
    assert selected == list(number)


def test_unbounded_and_empty_chunks():
    assert chunk_rows(5, -1, -1) == range(0, 5)
    assert chunk_rows(5, 3, -1) == range(2, 5)
    assert chunk_rows(5, -1, 2) == range(0, 2)
    # The splitter's canonical empty range (n_obj + 1, n_obj).
    assert len(chunk_rows(5, 6, 5)) == 0


def _write_sexcat(run_dir, number, epoch_numbers, ccd_n):
    out = run_dir / "output/run_sp_tile_Sx/sextractor_runner/output"
    out.mkdir(parents=True)
    hdus = [
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns(
            [fits.Column(name="NUMBER", format="K", array=number)],
            name="LDAC_OBJECTS",
        ),
    ]
    for k, (enum, ccd) in enumerate(zip(epoch_numbers, ccd_n)):
        hdus.append(
            fits.BinTableHDU.from_columns(
                [
                    fits.Column(name="NUMBER", format="K", array=enum),
                    fits.Column(name="CCD_N", format="K", array=ccd),
                ],
                name=f"EPOCH_{k}",
            )
        )
    fits.HDUList(hdus).writeto(out / "sexcat-301-279.fits")


def test_object_epochs_accepts_gapped_unsorted_number(tmp_path):
    number = np.array([40, 7, 1000, 3])
    ccd_n = [np.array([0, -1, 5, 2]), np.array([1, 1, -1, -1])]
    _write_sexcat(tmp_path, number, [number, number], ccd_n)
    np.testing.assert_array_equal(
        ngmix_range.object_epochs(tmp_path), [2, 1, 1, 1]
    )


def test_object_epochs_refuses_misaligned_epoch_rows(tmp_path):
    number = np.array([40, 7, 1000, 3])
    ccd_n = [np.zeros(4, dtype=int), np.zeros(4, dtype=int)]
    _write_sexcat(tmp_path, number, [number, number[::-1]], ccd_n)
    with pytest.raises(SystemExit, match="aligned"):
        ngmix_range.object_epochs(tmp_path)
