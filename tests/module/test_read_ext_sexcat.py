"""UNIT TESTS FOR MODULE PACKAGE: READ_EXT_SEXCAT.

Drives ``make_ldac_from_ascii`` on a synthetic ASCII SExtractor-format
catalogue and a synthetic tile image, and checks the FITS-LDAC it writes is
what the tile chain downstream of ``tile_detect`` reads: the LDAC_IMHEAD
extension carrying the tile header, the SExtractor column aliases, one
``VIGNET`` stamp per object cut from the image, and ``TILE_UNIQUE_ID``. The
last test follows that column through ``make_cat.save_sextractor_data`` into
the final catalogue.
"""

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits

from shapepipe.modules.make_cat_package import make_cat
from shapepipe.modules.read_ext_sexcat_package import read_ext_sexcat as rs

NX, NY = 40, 30
STAMP = 5
# (NUMBER, X_IMAGE, Y_IMAGE): an interior object, one on the left edge, one
# in the top-right corner.
OBJECTS = [(1, 10.0, 12.0), (2, 1.0, 20.0), (7, 40.0, 30.0)]


def _write_ascii_cat(path):
    lines = [
        "#   1 NUMBER          Running object number",
        "#   2 X_IMAGE         Object position along x        [pixel]",
        "#   3 Y_IMAGE         Object position along y        [pixel]",
        "#   4 ALPHA_J2000     Right ascension of barycenter  [deg]",
        "#   5 DELTA_J2000     Declination of barycenter      [deg]",
        "#   6 MAG_AUTO        Kron-like elliptical aperture magnitude [mag]",
    ]
    for num, x, y in OBJECTS:
        lines.append(f"{num} {x} {y} {150.0 + num} {30.0 + num} {20.0 + num}")
    path.write_text("\n".join(lines) + "\n")


def _write_image(path):
    # pixel value = 1000*row + column (0-based), so every stamp pixel names
    # where it came from.
    data = (np.arange(NY)[:, None] * 1000 + np.arange(NX)[None, :]).astype(
        np.float32
    )
    hdu = fits.PrimaryHDU(data)
    hdu.header["HISTORY"] = "input image 2605805p.fits"
    hdu.header["TILEKEY"] = "kept"
    hdu.writeto(path, overwrite=True)


@pytest.fixture
def ldac(tmp_path):
    cat = tmp_path / "CFIS_cat-301-279.cat"
    img = tmp_path / "CFIS_image-301-279.fits"
    out = tmp_path / "sexcat-301-279.fits"
    _write_ascii_cat(cat)
    _write_image(img)
    rs.make_ldac_from_ascii(
        str(cat), str(img), str(out), "-301-279", stamp_size=STAMP
    )
    return out


def test_ldac_layout_and_header(ldac):
    with fits.open(ldac) as hdul:
        assert [h.name for h in hdul] == ["PRIMARY", "LDAC_IMHEAD", "LDAC_OBJECTS"]
        cards = hdul["LDAC_IMHEAD"].data[0][0]
        assert isinstance(cards, str) or cards.ndim == 1
        text = "".join(cards) if not isinstance(cards, str) else cards
        assert "TILEKEY" in text and "2605805p" in text


def test_tile_unique_id_and_aliases(ldac):
    with fits.open(ldac) as hdul:
        data = hdul["LDAC_OBJECTS"].data
    numbers = np.array([o[0] for o in OBJECTS])
    # NUMBER is renumbered to a contiguous 1..n_obj running index in output
    # row order; the original NUMBER survives only inside TILE_UNIQUE_ID.
    npt.assert_array_equal(data["NUMBER"], np.arange(1, len(OBJECTS) + 1))
    npt.assert_array_equal(data["TILE_UNIQUE_ID"], 301279 * 10**6 + numbers)
    # A FITS 'K' column reads back as big-endian ('>i8'), not np.int64's
    # native byte order, so compare kind and width rather than dtype.
    assert data["TILE_UNIQUE_ID"].dtype.kind == "i"
    assert data["TILE_UNIQUE_ID"].dtype.itemsize == 8
    npt.assert_array_equal(data["XWIN_IMAGE"], data["X_IMAGE"])
    npt.assert_array_equal(data["YWIN_IMAGE"], data["Y_IMAGE"])
    npt.assert_array_equal(data["XWIN_WORLD"], data["ALPHA_J2000"])
    npt.assert_array_equal(data["YWIN_WORLD"], data["DELTA_J2000"])


def test_vignets_are_cut_from_the_image_and_zero_padded(ldac):
    with fits.open(ldac) as hdul:
        vignets = hdul["LDAC_OBJECTS"].data["VIGNET"]
    assert vignets.shape == (len(OBJECTS), STAMP, STAMP)

    # Interior object at (10, 12), 1-based: centre pixel is (row 11, col 9).
    assert vignets[0, STAMP // 2, STAMP // 2] == 11 * 1000 + 9
    assert vignets[0, 0, 0] == 9 * 1000 + 7

    # Left edge, x = 1: the two columns left of the image are zero.
    assert (vignets[1, :, :2] == 0).all()
    assert vignets[1, STAMP // 2, STAMP // 2] == 19 * 1000 + 0

    # Top-right corner: only the lower-left quadrant of the stamp is in the
    # image.
    assert (vignets[2, STAMP // 2 + 1:, :] == 0).all()
    assert (vignets[2, :, STAMP // 2 + 1:] == 0).all()
    assert vignets[2, STAMP // 2, STAMP // 2] == 29 * 1000 + 39


@pytest.mark.parametrize("number", ["-1301-279", "-301-1279"])
def test_four_digit_grid_index_is_rejected(number):
    with pytest.raises(ValueError):
        rs._tile_id_from_file_number_string(number)


def test_tile_unique_id_reaches_the_final_catalogue(ldac, tmp_path):
    """make_cat copies every sexcat column, so the ID needs no extra wiring."""
    final = make_cat.prepare_final_cat_file(str(tmp_path), "-301-279")
    n_obj = make_cat.save_sextractor_data(final, str(ldac))
    assert n_obj == len(OBJECTS)
    with fits.open(tmp_path / "final_cat-301-279.fits") as hdul:
        data = hdul["RESULTS"].data
    assert "VIGNET" not in data.names
    npt.assert_array_equal(
        data["TILE_UNIQUE_ID"], 301279 * 10**6 + np.array([1, 2, 7])
    )
    npt.assert_allclose(data["TILE_ID"], 301.279)
