"""Chunk-path derivation in the merge_sep_cats module."""

import pytest

from shapepipe.modules.merge_sep_cats_package.merge_sep_cats import chunk_path


REL = "./output/run_sp_tile_ngmix_Ng1u/ngmix_runner/output/ngmix-210-282.fits"
ABS = (
    "/scratch/run/tiles/21/210.282/output/run_sp_tile_ngmix_Ng1u"
    "/ngmix_runner/output/ngmix-210-282.fits"
)


def test_relative_path_unchanged_behaviour():
    """The bash pipeline's relative INPUT_DIR still resolves as before."""
    assert chunk_path(REL, 3) == (
        "./output/run_sp_tile_ngmix_Ng3u/ngmix_runner/output/ngmix-210-282.fits"
    )


def test_sharded_absolute_path():
    """Digits in the sharded parent dirs and in the file number are untouched."""
    assert chunk_path(ABS, 2) == (
        "/scratch/run/tiles/21/210.282/output/run_sp_tile_ngmix_Ng2u"
        "/ngmix_runner/output/ngmix-210-282.fits"
    )


def test_double_digit_chunk():
    assert "Ng12u" in chunk_path(ABS, 12)


def test_chunk_one_is_identity():
    assert chunk_path(ABS, 1) == ABS


def test_no_run_directory_raises():
    with pytest.raises(ValueError, match="run_"):
        chunk_path("/scratch/run/tiles/21/210.282/ngmix-210-282.fits", 2)
