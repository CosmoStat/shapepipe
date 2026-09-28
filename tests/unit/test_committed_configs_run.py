"""Exercise committed SExtractor/PSFEx configs against the real tools."""

from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS


ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "workflow" / "config" / "cfis"
SEX = shutil.which("source-extractor")
PSFEX = shutil.which("psfex")


def _synthetic_images(directory):
    """Write a WCS image with isolated stars and the maps used by both runners."""
    size = 512
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [size / 2, size / 2]
    wcs.wcs.cdelt = [-5.16e-5, 5.16e-5]
    wcs.wcs.crval = [180.0, 30.0]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["GAIN"] = 1.0
    header["SATURATE"] = 60000.0
    header["PHOTZP"] = 30.0

    rng = np.random.default_rng(42)
    yy, xx = np.mgrid[:size, :size]
    image = rng.normal(100.0, 3.0, (size, size)).astype(np.float32)
    for x in np.linspace(55, size - 55, 4):
        for y in np.linspace(55, size - 55, 4):
            image += 5000.0 * np.exp(
                -((xx - x) ** 2 + (yy - y) ** 2) / 8.0
            )

    fits.PrimaryHDU(image, header).writeto(directory / "image.fits")
    fits.PrimaryHDU(np.ones_like(image)).writeto(directory / "weight.fits")
    fits.PrimaryHDU(np.zeros_like(image, dtype=np.int16)).writeto(
        directory / "flag.fits"
    )


def _run_sextractor(directory, *, tile):
    config = CONFIG_DIR / ("default_tile.sex" if tile else "default_exp.sex")
    parameters = CONFIG_DIR / (
        "default_noimaflags.param" if tile else "default.param"
    )
    convolution = CONFIG_DIR / (
        "gauss_3.0_7x7.conv" if tile else "default.conv"
    )
    catalogue = directory / ("tile.cat" if tile else "exposure.cat")
    command = [
        SEX,
        "image.fits",
        "-c", str(config),
        "-PARAMETERS_NAME", str(parameters),
        "-FILTER_NAME", str(convolution),
        "-CATALOG_NAME", str(catalogue),
        "-WEIGHT_IMAGE", "weight.fits",
        "-FLAG_IMAGE", "NONE" if tile else "flag.fits",
        "-VERBOSE_TYPE", "QUIET",
        "-WRITE_XML", "N",
    ]
    result = subprocess.run(
        command,
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
        timeout=8,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with fits.open(catalogue) as hdus:
        objects = next(hdu for hdu in hdus if hdu.name == "LDAC_OBJECTS")
        assert objects.data is not None and len(objects.data) > 0
    return catalogue


def test_committed_convolution_filters_start_with_sextractor_directive():
    """SExtractor requires CONV on line 1 of every convolution file."""
    conv_files = sorted(CONFIG_DIR.rglob("*.conv"))
    assert conv_files
    for path in conv_files:
        assert path.read_text(encoding="utf-8").splitlines()[0].startswith("CONV "), path


@pytest.mark.skipif(not (SEX and PSFEX), reason="source-extractor and psfex are required")
def test_committed_sextractor_and_psfex_configs_run(tmp_path):
    """Both detection configs emit catalogues and the exposure catalogue fits a PSF."""
    _synthetic_images(tmp_path)
    _run_sextractor(tmp_path, tile=True)
    exposure_catalogue = _run_sextractor(tmp_path, tile=False)

    subprocess.run(
        [
            PSFEX,
            str(exposure_catalogue),
            "-c", str(CONFIG_DIR / "default.psfex"),
            # Smaller fixture stamp keeps this committed-config smoke test fast.
            "-PSF_SIZE", "21,21",
            "-PSF_SUFFIX", ".psf",
            "-WRITE_XML", "N",
            "-CHECKIMAGE_TYPE", "NONE",
            "-CHECKPLOT_TYPE", "NONE",
            "-VERBOSE_TYPE", "QUIET",
        ],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
        timeout=8,
    )
    psf_file = exposure_catalogue.with_suffix(".psf")
    assert psf_file.is_file()
    assert psf_file.stat().st_size > 0
