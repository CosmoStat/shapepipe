"""INTEGRATION TEST: TPV distortion through external tools.

This test verifies that the `split_exp` → SExtractor → PSFEx pipeline
chain handles TPV distortion keywords natively, without shapepipe
converting to SIP at the keyword level. Introduced alongside the
removal of the `sip_tpv` dependency.

Requires `source-extractor` and `psfex` on PATH; skipped otherwise.

"""

import os
import shutil
import subprocess
import tempfile
from unittest import TestCase, skipUnless

import numpy as np
import numpy.testing as npt
from astropy.io import fits
from astropy.wcs import WCS

HAS_SEX = shutil.which("source-extractor") is not None
HAS_PSFEX = shutil.which("psfex") is not None


def _tpv_header(nx, ny, crval1=180.0, crval2=30.0):
    """Synthetic CFIS-scale TPV header with a non-trivial radial term."""
    h = fits.Header()
    h["NAXIS"], h["NAXIS1"], h["NAXIS2"] = 2, nx, ny
    h["CTYPE1"], h["CTYPE2"] = "RA---TPV", "DEC--TPV"
    h["CUNIT1"], h["CUNIT2"] = "deg", "deg"
    h["CRPIX1"], h["CRPIX2"] = nx / 2, ny / 2
    h["CRVAL1"], h["CRVAL2"] = crval1, crval2
    scale = 5.16e-5
    h["CD1_1"], h["CD1_2"] = -scale, 0.0
    h["CD2_1"], h["CD2_2"] = 0.0, scale
    for axis in (1, 2):
        h[f"PV{axis}_0"] = 0.0
        h[f"PV{axis}_1"] = 1.0
        h[f"PV{axis}_2"] = 0.0
        h[f"PV{axis}_4"] = 2.0e-4
        h[f"PV{axis}_5"] = 0.0
        h[f"PV{axis}_6"] = 2.0e-4
    return h


@skipUnless(HAS_SEX, "source-extractor not on PATH")
class SExtractorTPVTestCase(TestCase):
    """SExtractor 2.25 parses TPV headers consistently with astropy."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        nx = ny = 512
        self.header = _tpv_header(nx, ny)
        self.src_xy = np.array(
            [
                (80, 80),
                (80, 256),
                (80, 432),
                (256, 80),
                (256, 256),
                (256, 432),
                (432, 80),
                (432, 256),
                (432, 432),
            ],
            dtype=float,
        )
        rng = np.random.default_rng(42)
        img = rng.normal(100.0, 3.0, (ny, nx)).astype(np.float32)
        yy, xx = np.mgrid[0:ny, 0:nx]
        for xi, yi in self.src_xy:
            img += 5000.0 * np.exp(-((xx - xi) ** 2 + (yy - yi) ** 2) / 8.0)
        self.image_path = os.path.join(self.tmpdir, "image.fits")
        fits.PrimaryHDU(data=img, header=self.header).writeto(self.image_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_sextractor_world_coords_match_astropy(self):
        cat_path = os.path.join(self.tmpdir, "sex.cat")
        param_path = os.path.join(self.tmpdir, "test.param")
        with open(param_path, "w") as f:
            f.write("NUMBER\nX_IMAGE\nY_IMAGE\nALPHA_J2000\nDELTA_J2000\n")

        subprocess.run(
            [
                "source-extractor",
                self.image_path,
                "-CATALOG_NAME",
                cat_path,
                "-CATALOG_TYPE",
                "ASCII_HEAD",
                "-PARAMETERS_NAME",
                param_path,
                "-DETECT_MINAREA",
                "3",
                "-DETECT_THRESH",
                "5.0",
                "-FILTER",
                "N",
                "-VERBOSE_TYPE",
                "QUIET",
            ],
            check=True,
            cwd=self.tmpdir,
        )

        rows = [
            line.split()
            for line in open(cat_path)
            if not line.startswith("#") and line.strip()
        ]
        sex = np.array(
            [
                [float(r[1]), float(r[2]), float(r[3]), float(r[4])]
                for r in rows
            ]
        )
        self.assertEqual(len(sex), len(self.src_xy))

        truth_fits_xy = self.src_xy + 1.0  # SExtractor uses 1-indexed pixels
        truth_world = WCS(self.header).wcs_pix2world(truth_fits_xy, 1)

        for tx, ty, tra, tdec in zip(
            truth_fits_xy[:, 0],
            truth_fits_xy[:, 1],
            truth_world[:, 0],
            truth_world[:, 1],
        ):
            i = int(np.argmin((sex[:, 0] - tx) ** 2 + (sex[:, 1] - ty) ** 2))
            dra_mas = (
                (sex[i, 2] - tra) * 3_600_000.0 * np.cos(np.radians(tdec))
            )
            ddec_mas = (sex[i, 3] - tdec) * 3_600_000.0
            sep_mas = np.hypot(dra_mas, ddec_mas)
            # Pixel scale ~186 mas/pix; require <10 mas (~0.05 px).
            self.assertLess(
                sep_mas,
                10.0,
                f"SExtractor world coord at ({tx}, {ty}) diverges from "
                f"astropy by {sep_mas:.3f} mas",
            )


@skipUnless(HAS_SEX and HAS_PSFEX, "source-extractor and psfex required")
class PSFExTPVTestCase(TestCase):
    """PSFEx consumes a TPV-headed SExtractor LDAC without error."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        nx = ny = 512
        self.header = _tpv_header(nx, ny)
        xs, ys = np.meshgrid(
            np.linspace(50, nx - 50, 6), np.linspace(50, ny - 50, 6)
        )
        self.src_xy = np.column_stack([xs.ravel(), ys.ravel()])
        rng = np.random.default_rng(7)
        img = rng.normal(100.0, 3.0, (ny, nx)).astype(np.float32)
        yy, xx = np.mgrid[0:ny, 0:nx]
        for xi, yi in self.src_xy:
            img += 5000.0 * np.exp(-((xx - xi) ** 2 + (yy - yi) ** 2) / 8.0)
        self.image_path = os.path.join(self.tmpdir, "image.fits")
        fits.PrimaryHDU(data=img, header=self.header).writeto(self.image_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_psfex_runs_on_tpv_ldac(self):
        param_path = os.path.join(self.tmpdir, "psfex.param")
        with open(param_path, "w") as f:
            f.write(
                "NUMBER\nX_IMAGE\nY_IMAGE\nXWIN_IMAGE\nYWIN_IMAGE\n"
                "FLUX_RADIUS\nFLUX_APER(1)\nFLUXERR_APER(1)\n"
                "FLUX_AUTO\nFLUXERR_AUTO\nSNR_WIN\nVIGNET(35,35)\n"
                "FLAGS\nFLAGS_WEIGHT\nELONGATION\n"
                "ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\n"
            )
        ldac_path = os.path.join(self.tmpdir, "sex.ldac")
        subprocess.run(
            [
                "source-extractor",
                self.image_path,
                "-CATALOG_NAME",
                ldac_path,
                "-CATALOG_TYPE",
                "FITS_LDAC",
                "-PARAMETERS_NAME",
                param_path,
                "-DETECT_MINAREA",
                "3",
                "-DETECT_THRESH",
                "5.0",
                "-FILTER",
                "N",
                "-PHOT_APERTURES",
                "8",
                "-VERBOSE_TYPE",
                "QUIET",
            ],
            check=True,
            cwd=self.tmpdir,
        )

        # LDAC_IMHEAD must preserve TPV keywords verbatim.
        with fits.open(ldac_path) as h:
            cards = [str(c) for c in h["LDAC_IMHEAD"].data[0][0]]
        preserved = {"CTYPE1", "CTYPE2", "PV1_4", "PV2_4", "PV1_6", "PV2_6"}
        for key in preserved:
            self.assertTrue(
                any(key in c for c in cards),
                f"LDAC_IMHEAD missing {key}; SExtractor dropped TPV keys",
            )

        subprocess.run(
            [
                "psfex",
                ldac_path,
                "-PSF_SUFFIX",
                ".psf",
                "-WRITE_XML",
                "N",
                "-CHECKIMAGE_TYPE",
                "NONE",
                "-CHECKPLOT_TYPE",
                "NONE",
                "-SAMPLE_MINSN",
                "5",
                "-SAMPLE_AUTOSELECT",
                "N",
                "-SAMPLE_FWHMRANGE",
                "2.0,8.0",
                "-VERBOSE_TYPE",
                "QUIET",
            ],
            check=True,
            cwd=self.tmpdir,
        )

        psf_path = os.path.join(self.tmpdir, "sex.psf")
        self.assertTrue(os.path.exists(psf_path))
        with fits.open(psf_path) as h:
            hdr = h[1].header
            self.assertGreater(hdr["ACCEPTED"], 0)
            self.assertLess(hdr["CHI2"], 5.0)

        # Independently: astropy must read the LDAC-embedded header and
        # produce a WCS that round-trips to the source header.
        with fits.open(ldac_path) as h:
            cards = [str(c) for c in h["LDAC_IMHEAD"].data[0][0]]
        embedded = fits.Header.fromstring("\n".join(cards), sep="\n")
        w_ldac = WCS(embedded)
        w_src = WCS(self.header)
        pix = np.array([[100.0, 100.0], [256.0, 256.0], [400.0, 400.0]])
        npt.assert_allclose(
            w_ldac.wcs_pix2world(pix, 1),
            w_src.wcs_pix2world(pix, 1),
            atol=1e-9,
            err_msg="LDAC-embedded WCS diverges from source header",
        )
