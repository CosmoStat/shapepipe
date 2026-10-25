"""The committed PSF star selection bounds sizes in arcsec via FWHM_WORLD."""

from pathlib import Path

import numpy as np

from shapepipe.modules.setools_package.setools import SETools

CONFIG = (
    Path(__file__).resolve().parents[2]
    / "workflow/config/cfis/star_selection.setools"
)


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


def test_star_selection_preselects_on_angular_fwhm(tmp_path):
    """0.75" stars are selected; more numerous 2" galaxies are not.

    Without a working angular upper bound in the preselection, the FWHM mode
    locks onto the galaxies and none of the stars are selected.
    """
    rng = np.random.default_rng(1)
    n_star, n_gal = 30, 40
    fwhm_arcsec = np.r_[np.full(n_star, 0.75), np.full(n_gal, 2.0)]
    pixel_scale = rng.uniform(0.185, 0.188, fwhm_arcsec.size)
    fwhm_image = fwhm_arcsec / pixel_scale + rng.uniform(-0.05, 0.05, fwhm_arcsec.size)

    catalogue = np.zeros(
        fwhm_image.size,
        dtype=[
            ("FWHM_IMAGE", "f8"),
            ("FWHM_WORLD", "f8"),
            ("MAG_AUTO", "f8"),
            ("FLAGS", "i4"),
            ("IMAFLAGS_ISO", "i4"),
            ("X_IMAGE", "f8"),
            ("Y_IMAGE", "f8"),
        ],
    )
    catalogue["FWHM_IMAGE"] = fwhm_image
    catalogue["FWHM_WORLD"] = fwhm_image * pixel_scale / 3600.0
    catalogue["MAG_AUTO"] = rng.uniform(19.0, 20.5, fwhm_image.size)
    catalogue["X_IMAGE"] = rng.uniform(0, 2048, fwhm_image.size)
    catalogue["Y_IMAGE"] = rng.uniform(0, 4612, fwhm_image.size)

    tools = SETools(catalogue, str(tmp_path), "-test", str(CONFIG), cat_file=False)
    try:
        tools.process(_NullLogger())
    finally:
        tools._config_file.close()

    np.testing.assert_array_equal(
        tools.mask["star_selection"], np.arange(fwhm_image.size) < n_star
    )
    stats = dict(
        line.split(" = ")
        for line in (tmp_path / "stat/star_stat.txt").read_text().splitlines()
        if " = " in line
    )
    assert float(stats["Mean star fwhm selected (arcsec)"]) == np.float64(
        np.mean(fwhm_image[:n_star] * pixel_scale[:n_star])
    ).round(6) or abs(
        float(stats["Mean star fwhm selected (arcsec)"]) - 0.75
    ) < 0.02
