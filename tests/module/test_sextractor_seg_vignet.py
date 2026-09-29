"""UNIT TESTS FOR SEXTRACTOR MODE'S SEG_VIGNET.

``sextractor_script.add_seg_vignet`` cuts the SEGMENTATION check image on the
grid SExtractor cut each object's VIGNET on, and adds it to the sexcat as the
int32 ``SEG_VIGNET`` column. SExtractor centres VIGNET on the pixel nearest
its double-precision barycentre; the catalogue stores X_IMAGE / Y_IMAGE as
float32, which ties a few positions to an exact half pixel. Those objects
take whichever neighbouring centre reproduces their VIGNET from the image
minus the BACKGROUND check image, bit for bit.
"""

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits

from shapepipe.modules.sextractor_package import sextractor_script as ss

NX, NY = 60, 50
STAMP = 7
BIG = np.float32(-1e30)

# (NUMBER, x, y): 1-based double-precision barycentres. Object 3's x and
# object 4's y sit a float32 ulp either side of a half pixel, so the float32
# catalogue shows an exact .5 for both but SExtractor rounded them apart.
OBJECTS = [
    (1, 20.2, 15.7),
    (2, 1.3, 48.9),
    (3, 30.5 + 5e-7, 30.1),
    (4, 44.8, 22.5 - 5e-7),
    (5, 59.9, 1.2),
]


def _scene(tmp_path):
    rng = np.random.default_rng(3)
    image = rng.normal(100, 5, (NY, NX)).astype(np.float32)
    background = np.linspace(99, 101, NX * NY, dtype=np.float32).reshape(
        NY, NX
    )
    seg = np.zeros((NY, NX), np.int32)
    for number, x, y in OBJECTS:
        col, row = int(np.rint(x - 1)), int(np.rint(y - 1))
        seg[max(row - 2, 0):row + 3, max(col - 2, 0):col + 3] = number
    number = np.array([o[0] for o in OBJECTS])
    x = np.array([o[1] for o in OBJECTS])
    y = np.array([o[2] for o in OBJECTS])
    col, row = np.rint(x - 1).astype(int), np.rint(y - 1).astype(int)

    # VIGNET as SExtractor writes it: image minus background, -1e30 off the
    # image and on some neighbours' pixels.
    vignets = ss.cut_stamps(image - background, col, row, STAMP, BIG)
    seg_true = ss.cut_stamps(seg, col, row, STAMP, 0)
    vignets[(seg_true != 0) & (seg_true != number[:, None, None])] = BIG

    paths = {name: str(tmp_path / f"{name}-001-001.fits")
             for name in ("image", "background", "segmentation", "sexcat")}
    fits.PrimaryHDU(image).writeto(paths["image"])
    fits.PrimaryHDU(background).writeto(paths["background"])
    fits.PrimaryHDU(seg).writeto(paths["segmentation"])
    objects = fits.BinTableHDU.from_columns([
        fits.Column(name="NUMBER", format="J", array=number),
        fits.Column(name="X_IMAGE", format="E", array=x.astype(np.float32)),
        fits.Column(name="Y_IMAGE", format="E", array=y.astype(np.float32)),
        fits.Column(name="VIGNET", format=f"{STAMP * STAMP}E",
                    array=vignets.reshape(len(number), -1),
                    dim=f"({STAMP},{STAMP})"),
    ], name="LDAC_OBJECTS")
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="80A",
                     array=np.array(["HISTORY x"]))], name="LDAC_IMHEAD")
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(
        paths["sexcat"])
    return paths, seg, seg_true


def _add(paths):
    ss.add_seg_vignet(paths["sexcat"], paths["segmentation"], paths["image"],
                      paths["background"])
    with fits.open(paths["sexcat"]) as hdul:
        return [h.name for h in hdul], hdul["LDAC_OBJECTS"].data.copy()


def test_float32_catalogue_positions_tie_on_half_pixels():
    """The premise: the float32 columns lose which way objects 3 and 4
    rounded, so no rounding rule on the catalogue alone recovers both."""
    x3 = np.float32(OBJECTS[2][1])
    y4 = np.float32(OBJECTS[3][2])
    assert x3 % 1 == 0.5 and y4 % 1 == 0.5
    assert np.rint(OBJECTS[2][1] - 1) == np.floor(x3 - 1) + 1
    assert np.rint(OBJECTS[3][2] - 1) == np.floor(y4 - 1)


def test_seg_vignet_is_registered_with_vignet(tmp_path):
    """SEG_VIGNET is the check image on VIGNET's grid for every object, the
    half-pixel ties included: int32, VIGNET's shape, 0 off the image."""
    paths, _, seg_true = _scene(tmp_path)
    names, data = _add(paths)
    assert names == ["PRIMARY", "LDAC_IMHEAD", "LDAC_OBJECTS"]
    seg_vignets = data["SEG_VIGNET"]
    assert seg_vignets.dtype.kind == "i" and seg_vignets.dtype.itemsize == 4
    assert seg_vignets.shape == data["VIGNET"].shape
    npt.assert_array_equal(seg_vignets, seg_true)
    centre = STAMP // 2
    npt.assert_array_equal(seg_vignets[:, centre, centre], data["NUMBER"])


def test_seg_vignet_leaves_the_other_columns_alone(tmp_path):
    paths, _, _ = _scene(tmp_path)
    with fits.open(paths["sexcat"]) as hdul:
        before = hdul["LDAC_OBJECTS"].data.copy()
        imhead = hdul["LDAC_IMHEAD"].data.copy()
    _, after = _add(paths)
    assert after.names == before.names + ["SEG_VIGNET"]
    for name in before.names:
        npt.assert_array_equal(after[name], before[name])
    with fits.open(paths["sexcat"]) as hdul:
        npt.assert_array_equal(hdul["LDAC_IMHEAD"].data, imhead)


def test_sextractor_caller_names_its_check_images(tmp_path):
    """The runner finds the check images add_seg_vignet reads by type."""
    caller = ss.SExtractorCaller(
        [str(tmp_path / "image-001-001.fits"),
         str(tmp_path / "weight-001-001.fits")],
        str(tmp_path), "-001-001", "d.sex", "d.param", "d.conv",
        True, False, False, False, False, False, False,
        check_image=["BACKGROUND", "SEGMENTATION"],
    )
    assert caller.check_paths == {
        "BACKGROUND": f"{tmp_path}/background-001-001.fits",
        "SEGMENTATION": f"{tmp_path}/segmentation-001-001.fits",
    }


def test_seg_vignet_refuses_a_map_on_another_grid(tmp_path):
    paths, seg, _ = _scene(tmp_path)
    fits.PrimaryHDU(seg[:-1]).writeto(paths["segmentation"], overwrite=True)
    with pytest.raises(ValueError, match="grid"):
        _add(paths)
