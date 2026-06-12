"""UNIT TESTS FOR SEXTRACTOR POST-PROCESS.

Regression test for the tile-mode header log: ``merge_headers`` writes a
``TILE_ID`` metadata entry as the *first* key of
``log_exp_headers<tile>.sqlite``, and ``make_post_process`` derives the
number of CCDs from the first key's value. With the metadata key first,
``n_hdu`` becomes ``len(tile_id_string)`` instead of the number of CCDs,
silently dropping every epoch on the unscanned CCDs.

"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from shapepipe.modules.merge_headers_package.merge_headers import merge_headers
from shapepipe.modules.sextractor_package import sextractor_script

N_CCD = 3
CCD_NPIX = 100


def _make_ccd_wcs(ccd):
    """TAN WCS per CCD, each covering a disjoint sky patch."""
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = CCD_NPIX
    header["NAXIS2"] = CCD_NPIX
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CRPIX1"] = 50.0
    header["CRPIX2"] = 50.0
    header["CRVAL1"] = 180.0 + 0.5 * ccd
    header["CRVAL2"] = 30.0
    header["CD1_1"] = -1.0e-3
    header["CD1_2"] = 0.0
    header["CD2_1"] = 0.0
    header["CD2_2"] = 1.0e-3
    return WCS(header), header


def _write_exposure_headers(path):
    """Mimic split_exp: object array of {WCS, header} dicts, one per CCD."""
    entries = np.zeros(N_CCD, dtype="O")
    for ccd in range(N_CCD):
        wcs, header = _make_ccd_wcs(ccd)
        entries[ccd] = {"WCS": wcs, "header": header.tostring()}
    np.save(path, entries)


def _write_sex_ldac(path, exp_names, world_positions):
    """Minimal SExtractor LDAC catalogue with HISTORY exposure cards."""
    cards = np.array(
        [[f"HISTORY input image {name}p.fits" for name in exp_names]],
        dtype="U80",
    )
    imhead = fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="Field Header Card",
                format=f"{80 * len(exp_names)}A",
                dim=f"(80,{len(exp_names)})",
                array=cards,
            )
        ],
        name="LDAC_IMHEAD",
    )
    objects = fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="NUMBER",
                format="J",
                array=np.arange(1, len(world_positions) + 1),
            ),
            fits.Column(
                name="XWIN_WORLD",
                format="D",
                array=world_positions[:, 0],
            ),
            fits.Column(
                name="YWIN_WORLD",
                format="D",
                array=world_positions[:, 1],
            ),
        ],
        name="LDAC_OBJECTS",
    )
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(path)


def test_post_process_scans_all_ccds_despite_tile_id_key(tmp_path):
    """All CCDs are scanned even though TILE_ID is the sqlite's first key.

    The tile number ("51", length 2) is shorter than the CCD count (3): if
    n_hdu is derived from the TILE_ID entry, CCD 2 is never scanned and the
    object there loses both its epochs.

    """
    exp_names = ["123456", "654321"]
    header_files = []
    for name in exp_names:
        npy_path = tmp_path / f"headers-{name}.npy"
        _write_exposure_headers(npy_path)
        header_files.append([str(npy_path)])

    tile_number = "51"
    merge_headers(header_files, str(tmp_path), tile_number=tile_number)
    sqlite_path = tmp_path / f"log_exp_headers{tile_number}.sqlite"
    assert sqlite_path.is_file()

    # One object at the centre of CCD 0, one at the centre of CCD 2.
    positions = np.array(
        [
            _make_ccd_wcs(ccd)[0].all_pix2world([[50.0, 50.0]], 0)[0]
            for ccd in (0, 2)
        ]
    )
    cat_path = tmp_path / "sexcat.fits"
    _write_sex_ldac(cat_path, exp_names, positions)

    sextractor_script.make_post_process(
        str(cat_path),
        str(sqlite_path),
        ["XWIN_WORLD", "YWIN_WORLD"],
        ["0", str(CCD_NPIX), "0", str(CCD_NPIX)],
    )

    with fits.open(cat_path) as hdu_list:
        n_epoch = hdu_list["LDAC_OBJECTS"].data["N_EPOCH"]
        ccd_n = {
            idx: hdu_list[f"EPOCH_{idx}"].data["CCD_N"]
            for idx in range(len(exp_names))
        }

    # Both objects appear in both exposures, on their true CCDs.
    np.testing.assert_array_equal(n_epoch, [2, 2])
    for idx in range(len(exp_names)):
        np.testing.assert_array_equal(ccd_n[idx], [0, 2])
