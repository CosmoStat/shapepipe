"""UNIT TESTS FOR PSFEX_INTERP."""

from types import SimpleNamespace

from astropy.io import fits
import numpy as np
import numpy.testing as npt
import pytest

from shapepipe.modules.psfex_interp_package.psfex_interp import (
    PSFExInterpolator,
)


def _write_galaxy_catalogue(path, **columns):

    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(),
            fits.BinTableHDU.from_columns(
                [
                    fits.Column(
                        name=name,
                        format="D",
                        array=np.asarray(values),
                    )
                    for name, values in columns.items()
                ]
            ),
        ]
    ).writeto(path)


def test_get_galaxy_positions_reads_requested_columns(tmp_path):

    cat_path = tmp_path / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, X_IMAGE=[1.5, 2.5], Y_IMAGE=[3.5, 4.5])
    interp = PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(tmp_path),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["X_IMAGE", "Y_IMAGE"],
    )

    interp._get_galaxy_positions()

    assert interp.gal_pos.shape == (2, 2)
    npt.assert_allclose(interp.gal_pos, [[1.5, 3.5], [2.5, 4.5]])


def test_get_galaxy_positions_raises_for_missing_column(tmp_path):

    cat_path = tmp_path / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, X_IMAGE=[1.5, 2.5])
    interp = PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(tmp_path),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["X_IMAGE", "Y_IMAGE"],
    )

    with pytest.raises(KeyError, match="Y_IMAGE"):
        interp._get_galaxy_positions()
