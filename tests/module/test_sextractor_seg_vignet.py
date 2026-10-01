"""UNIT TESTS FOR SEXTRACTOR MODE'S SEG_VIGNET.

``sextractor_script.add_seg_vignet`` cuts the SEGMENTATION check image on the
grid SExtractor cut each object's VIGNET on and adds it to the sexcat as the
int32 ``SEG_VIGNET`` column. SExtractor centres VIGNET on the pixel nearest
its double-precision barycentre. The catalogue's float32 X_IMAGE / Y_IMAGE
tie a few positions to an exact half pixel, so the runner asks SExtractor
for X_IMAGE_DBL / Y_IMAGE_DBL too (``seg_vignet_param_file``), and
``add_seg_vignet`` centres on them and then drops them.
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
# catalogue shows an exact .5 for both but SExtractor rounded them apart, and
# each the other way from rounding the float32 value half to even.
OBJECTS = [
    (1, 20.2, 15.7),
    (2, 1.3, 48.9),
    (3, 31.5 + 5e-7, 30.1),
    (4, 44.8, 22.5 - 5e-7),
    (5, 59.9, 1.2),
]
DBL = ["X_IMAGE_DBL", "Y_IMAGE_DBL"]


def _scene(tmp_path, flat=False):
    """A sexcat whose VIGNETs are cut as SExtractor cuts them.

    ``flat`` makes the image uniform, so no pixel comparison can tell a tie
    object's two candidate centres apart.
    """
    rng = np.random.default_rng(3)
    image = (np.full((NY, NX), 100, np.float32) if flat
             else rng.normal(100, 5, (NY, NX)).astype(np.float32))
    seg = np.zeros((NY, NX), np.int32)
    for number, x, y in OBJECTS:
        col, row = int(np.rint(x - 1)), int(np.rint(y - 1))
        seg[max(row - 2, 0):row + 3, max(col - 2, 0):col + 3] = number
    number = np.array([o[0] for o in OBJECTS])
    x = np.array([o[1] for o in OBJECTS])
    y = np.array([o[2] for o in OBJECTS])
    col, row = np.rint(x - 1).astype(int), np.rint(y - 1).astype(int)

    vignets = ss.cut_stamps(image, col, row, STAMP, BIG)
    seg_true = ss.cut_stamps(seg, col, row, STAMP, 0)
    vignets[(seg_true != 0) & (seg_true != number[:, None, None])] = BIG

    paths = {name: str(tmp_path / f"{name}-001-001.fits")
             for name in ("segmentation", "sexcat")}
    fits.PrimaryHDU(seg).writeto(paths["segmentation"])
    objects = fits.BinTableHDU.from_columns([
        fits.Column(name="NUMBER", format="J", array=number),
        fits.Column(name="X_IMAGE", format="E", array=x.astype(np.float32)),
        fits.Column(name="Y_IMAGE", format="E", array=y.astype(np.float32)),
        fits.Column(name="VIGNET", format=f"{STAMP * STAMP}E",
                    array=vignets.reshape(len(number), -1),
                    dim=f"({STAMP},{STAMP})"),
        fits.Column(name="THETA_J2000", format="E", array=np.zeros(5)),
        fits.Column(name="X_IMAGE_DBL", format="D", array=x),
        fits.Column(name="Y_IMAGE_DBL", format="D", array=y),
    ], name="LDAC_OBJECTS")
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="80A",
                     array=np.array(["HISTORY x"]))], name="LDAC_IMHEAD")
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(
        paths["sexcat"])
    return paths, seg_true


def _add(paths):
    ss.add_seg_vignet(paths["sexcat"], paths["segmentation"])
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
    assert np.rint(OBJECTS[2][1] - 1) != np.rint(np.float64(x3) - 1)
    assert np.rint(OBJECTS[3][2] - 1) != np.rint(np.float64(y4) - 1)


@pytest.mark.parametrize("flat", [False, True], ids=["noisy", "flat"])
def test_seg_vignet_is_registered_with_vignet(tmp_path, flat):
    """SEG_VIGNET is the check image on VIGNET's grid for every object, the
    half-pixel ties included, even where the image cannot tell the two
    candidate centres apart: int32, VIGNET's shape, 0 off the image."""
    paths, seg_true = _scene(tmp_path, flat=flat)
    names, data = _add(paths)
    assert names == ["PRIMARY", "LDAC_IMHEAD", "LDAC_OBJECTS"]
    seg_vignets = data["SEG_VIGNET"]
    assert seg_vignets.dtype.kind == "i" and seg_vignets.dtype.itemsize == 4
    assert seg_vignets.shape == data["VIGNET"].shape
    npt.assert_array_equal(seg_vignets, seg_true)
    centre = STAMP // 2
    npt.assert_array_equal(seg_vignets[:, centre, centre], data["NUMBER"])


def test_seg_vignet_replaces_the_double_positions(tmp_path):
    """The written catalogue has the columns SExtractor writes without the
    double positions, plus SEG_VIGNET, and every other HDU and column as
    they were."""
    paths, _ = _scene(tmp_path)
    with fits.open(paths["sexcat"]) as hdul:
        before = hdul["LDAC_OBJECTS"].data.copy()
        imhead = hdul["LDAC_IMHEAD"].data.copy()
    _, after = _add(paths)
    kept = [name for name in before.names if name not in DBL]
    assert after.names == kept + ["SEG_VIGNET"]
    for name in kept:
        npt.assert_array_equal(after[name], before[name])
    with fits.open(paths["sexcat"]) as hdul:
        npt.assert_array_equal(hdul["LDAC_IMHEAD"].data, imhead)


def test_seg_vignet_needs_the_double_positions(tmp_path):
    paths, _ = _scene(tmp_path)
    with fits.open(paths["sexcat"]) as hdul:
        objects = hdul["LDAC_OBJECTS"]
        cols = [c for c in objects.columns if c.name not in DBL]
        hdus = [hdul[0].copy(), hdul[1].copy(),
                fits.BinTableHDU.from_columns(cols, name="LDAC_OBJECTS")]
    fits.HDUList(hdus).writeto(paths["sexcat"], overwrite=True)
    with pytest.raises(ValueError, match="X_IMAGE_DBL"):
        _add(paths)


def test_param_file_adds_the_double_positions(tmp_path):
    """The runner's parameter file is the configured one plus the two
    double-precision positions, appended after its last parameter."""
    dot_param = tmp_path / "default.param"
    dot_param.write_text("NUMBER  # running number\nX_IMAGE\nVIGNET(51,51)")
    out = ss.seg_vignet_param_file(str(dot_param), str(tmp_path / "out.param"))
    lines = [line.split("#")[0].strip()
             for line in open(out).read().splitlines()]
    assert [line for line in lines if line] == [
        "NUMBER", "X_IMAGE", "VIGNET(51,51)"] + DBL
    assert dot_param.read_text().count("DBL") == 0


def test_sextractor_caller_names_its_check_images(tmp_path):
    """The runner finds the SEGMENTATION check image by type."""
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
