"""UNIT TESTS FOR FIND EXPOSURES.

Covers the epoch-provenance path in ``FindExposures.get_exposure_list``:
a tile whose HISTORY header cannot be read must fail deliberately, and the
exposure filename extension must be stripped in full, not just up to the
first extension.

"""

import logging

import astropy.io.fits as fits
import pytest

from shapepipe.modules.find_exposures_package.find_exposures import (
    FindExposures,
)


def _write_tile(path, history_lines=None):
    """Minimal tile-like FITS file, optionally with HISTORY cards."""
    header = fits.Header()
    if history_lines is not None:
        for line in history_lines:
            header["HISTORY"] = line
    fits.PrimaryHDU(header=header).writeto(path)


def _make_find_exposures(tmp_path, colnum, prefix=""):
    return FindExposures(
        img_tile_path=str(tmp_path / "tile.fits"),
        output_path=str(tmp_path / "exp_numbers.txt"),
        w_log=logging.getLogger("test_find_exposures"),
        colnum=colnum,
        prefix=prefix,
    )


def test_get_exposure_list_raises_when_history_missing(tmp_path):
    """A tile with no readable HISTORY fails deliberately with IOError."""
    _write_tile(tmp_path / "tile.fits", history_lines=None)
    find_exp = _make_find_exposures(tmp_path, colnum=2)

    with pytest.raises(IOError, match="tile.fits"):
        find_exp.get_exposure_list()


def test_get_exposure_list_strips_full_extension(tmp_path):
    """A multi-extension exposure name has every extension stripped."""
    _write_tile(
        tmp_path / "tile.fits",
        history_lines=["input image 2243881p.fits.fz 6 extension(s)"],
    )
    find_exp = _make_find_exposures(tmp_path, colnum=2)

    assert find_exp.get_exposure_list() == ["2243881p"]


def test_get_exposure_list_matches_committed_cfis_colnum(tmp_path):
    """The committed CFIS config's COLNUM=3 lands on the exposure filename.

    Production HISTORY values carry a leading space after the ``HISTORY``
    keyword (the FITS comment-field convention), so splitting on space
    shifts every column by one relative to a value with no leading space;
    ``config_tile_Fe.ini`` sets ``COLNUM = 3`` to land on the filename
    under that convention.
    """
    _write_tile(
        tmp_path / "tile.fits",
        history_lines=[" input image 2243881p.fits 6 extension(s)"],
    )
    find_exp = _make_find_exposures(tmp_path, colnum=3)

    assert find_exp.get_exposure_list() == ["2243881p"]


def test_get_exposure_list_keeps_trailing_epoch_letter(tmp_path):
    """CFIS exposure names keep their trailing ``p`` (no prefix configured)."""
    _write_tile(
        tmp_path / "tile.fits",
        history_lines=["input image 2243881p.fits 6 extension(s)"],
    )
    find_exp = _make_find_exposures(tmp_path, colnum=2, prefix="")

    assert find_exp.get_exposure_list() == ["2243881p"]


def test_get_exposure_list_strips_configured_prefix(tmp_path):
    """A real prefix (e.g. simulated exposures) is still stripped."""
    _write_tile(
        tmp_path / "tile.fits",
        history_lines=["input image simu_image-000-000.fits"],
    )
    find_exp = _make_find_exposures(
        tmp_path, colnum=2, prefix="simu_image-"
    )

    assert find_exp.get_exposure_list() == ["000-000"]
