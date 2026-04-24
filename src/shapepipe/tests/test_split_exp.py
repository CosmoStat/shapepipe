"""UNIT TESTS FOR SPLIT_EXP.

Verify that SplitExposures writes per-CCD FITS files whose WCS (including
TPV distortion) round-trips correctly via astropy without requiring the
sip_tpv package.

"""

import os
import tempfile
from unittest import TestCase

import numpy as np
import numpy.testing as npt
from astropy.io import fits
from astropy.wcs import WCS

from shapepipe.modules.split_exp_package.split_exp import SplitExposures


def _make_tpv_header(crval1, crval2):
    """Build a synthetic MegaCam-style header with TPV distortion."""
    h = fits.Header()
    h["NAXIS"] = 2
    h["NAXIS1"] = 64
    h["NAXIS2"] = 64
    h["CTYPE1"] = "RA---TPV"
    h["CTYPE2"] = "DEC--TPV"
    h["CRPIX1"] = 32.0
    h["CRPIX2"] = 32.0
    h["CRVAL1"] = crval1
    h["CRVAL2"] = crval2
    h["CD1_1"] = -5.16e-5
    h["CD1_2"] = 0.0
    h["CD2_1"] = 0.0
    h["CD2_2"] = 5.16e-5
    h["PV1_0"] = 0.0
    h["PV1_1"] = 1.0
    h["PV1_2"] = 0.0
    h["PV1_4"] = 1.0e-4
    h["PV2_0"] = 0.0
    h["PV2_1"] = 1.0
    h["PV2_2"] = 0.0
    h["PV2_4"] = 1.0e-4
    return h


class SplitExpWCSRoundTripTestCase(TestCase):
    """Split a synthetic multi-CCD exposure and verify WCS fidelity."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.n_hdu = 3
        self.file_number_string = "-0000001"
        self.crvals = [
            (180.0, 30.0),
            (180.1, 30.0),
            (180.0, 30.1),
        ]
        self.test_pixels = np.array([[1.0, 1.0], [32.0, 32.0], [60.0, 50.0]])

        self.exp_path = os.path.join(self.tmpdir, "exp.fits")
        primary = fits.PrimaryHDU()
        ccds = [
            fits.ImageHDU(
                data=np.zeros((64, 64), dtype=np.float32),
                header=_make_tpv_header(*self.crvals[i]),
            )
            for i in range(self.n_hdu)
        ]
        fits.HDUList([primary] + ccds).writeto(self.exp_path)

    def tearDown(self):
        for name in os.listdir(self.tmpdir):
            os.remove(os.path.join(self.tmpdir, name))
        os.rmdir(self.tmpdir)

    def test_per_ccd_wcs_matches_source(self):
        splitter = SplitExposures(
            input_file_list=[self.exp_path],
            output_dir=self.tmpdir,
            file_number_string=self.file_number_string,
            output_suffix=["image"],
            n_hdu=self.n_hdu,
        )
        splitter.process()

        with fits.open(self.exp_path) as src:
            for idx in range(self.n_hdu):
                src_world = WCS(src[idx + 1].header).wcs_pix2world(
                    self.test_pixels, 1
                )
                out_path = os.path.join(
                    self.tmpdir, f"image{self.file_number_string}-{idx}.fits"
                )
                with fits.open(out_path) as out:
                    out_wcs = WCS(out[0].header)
                    out_world = out_wcs.wcs_pix2world(self.test_pixels, 1)
                    npt.assert_allclose(
                        out_world,
                        src_world,
                        err_msg=f"CCD {idx}: world coords diverge from source",
                    )
                    npt.assert_allclose(
                        out_wcs.wcs_world2pix(out_world, 1),
                        self.test_pixels,
                        atol=1e-6,
                        err_msg=f"CCD {idx}: WCS does not round-trip",
                    )

    def test_headers_npy_roundtrips(self):
        splitter = SplitExposures(
            input_file_list=[self.exp_path],
            output_dir=self.tmpdir,
            file_number_string=self.file_number_string,
            output_suffix=["image"],
            n_hdu=self.n_hdu,
        )
        splitter.process()

        headers = np.load(
            os.path.join(self.tmpdir, f"headers{self.file_number_string}.npy"),
            allow_pickle=True,
        )
        self.assertEqual(len(headers), self.n_hdu)
        for idx, entry in enumerate(headers):
            w = entry["WCS"]
            world = w.wcs_pix2world(self.test_pixels, 1)
            npt.assert_allclose(
                w.wcs_world2pix(world, 1),
                self.test_pixels,
                atol=1e-6,
                err_msg=f"Stored WCS for CCD {idx} does not round-trip",
            )

    def test_flag_output_uses_int16(self):
        flag_exp = os.path.join(self.tmpdir, "flag.fits")
        fits.HDUList(
            [fits.PrimaryHDU()]
            + [
                fits.ImageHDU(
                    data=np.ones((64, 64), dtype=np.float32) * 7,
                    header=_make_tpv_header(*self.crvals[i]),
                )
                for i in range(self.n_hdu)
            ]
        ).writeto(flag_exp)

        splitter = SplitExposures(
            input_file_list=[flag_exp],
            output_dir=self.tmpdir,
            file_number_string=self.file_number_string,
            output_suffix=["flag"],
            n_hdu=self.n_hdu,
        )
        splitter.process()

        for idx in range(self.n_hdu):
            out_path = os.path.join(
                self.tmpdir, f"flag{self.file_number_string}-{idx}.fits"
            )
            with fits.open(out_path) as out:
                self.assertEqual(out[0].data.dtype.kind, "i")
                self.assertEqual(out[0].data.dtype.itemsize, 2)
