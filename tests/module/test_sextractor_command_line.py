"""SExtractor command line built by ``SExtractorCaller``.

Without a detection image the call is single-image (``img -WEIGHT_IMAGE w``);
dual-image mode (``det,img -WEIGHT_IMAGE det_w,w``) only with a detection
image, since dual-image mode on the same image twice changes the measurements.
"""

import shlex

import pytest

from shapepipe.modules.sextractor_package.sextractor_script import (
    SExtractorCaller,
)


def _command(inputs, **flags):
    options = dict(
        use_weight=False, use_flag=False, use_psf=False,
        use_detection_image=False, use_detection_weight=False,
        use_zero_point=False, use_background=False,
    )
    options.update(flags)
    caller = SExtractorCaller(
        inputs, "/out", "-000-000", "d.sex", "d.param", "d.conv",
        check_image=[""], **options,
    )
    argv = shlex.split(caller.make_command_line("source-extractor"))
    options = {argv[i]: argv[i + 1] for i in range(2, len(argv) - 1)
               if argv[i].startswith("-")}
    return argv, options


def test_tile_single_image_with_weight():
    argv, opt = _command(["img.fits", "w.fits"], use_weight=True)
    assert argv[:2] == ["source-extractor", "img.fits"]
    assert opt["-WEIGHT_IMAGE"] == "w.fits"
    assert "-FLAG_IMAGE" not in opt
    assert opt["-CATALOG_NAME"] == "/out/sexcat-000-000.fits"
    assert opt["-CHECKIMAGE_TYPE"] == "NONE"


def test_exposure_single_image_with_weight_and_flag():
    argv, opt = _command(["img.fits", "w.fits", "f.fits"],
                         use_weight=True, use_flag=True)
    assert argv[1] == "img.fits"
    assert opt["-WEIGHT_IMAGE"] == "w.fits"
    assert opt["-FLAG_IMAGE"] == "f.fits"


def test_single_image_without_weight():
    argv, opt = _command(["img.fits"])
    assert argv[1] == "img.fits"
    assert opt["-WEIGHT_TYPE"] == "None"
    assert "-WEIGHT_IMAGE" not in opt


@pytest.mark.parametrize("detection_weight", [False, True])
def test_dual_image_with_detection_image(detection_weight):
    inputs = ["img.fits", "w.fits", "det.fits"]
    if detection_weight:
        inputs.append("det_w.fits")
    argv, opt = _command(inputs, use_weight=True, use_detection_image=True,
                         use_detection_weight=detection_weight)
    assert argv[1] == "det.fits,img.fits"
    det_w = "det_w.fits" if detection_weight else "w.fits"
    assert opt["-WEIGHT_IMAGE"] == f"{det_w},w.fits"


def test_check_images_named_per_unit():
    caller = SExtractorCaller(
        ["img.fits", "w.fits"], "/out", "-001-002", "d.sex", "d.param",
        "d.conv", True, False, False, False, False, False, False,
        check_image=["BACKGROUND", "SEGMENTATION"],
    )
    argv = shlex.split(caller.make_command_line("source-extractor"))
    assert argv[1] == "img.fits"
    assert argv[argv.index("-CHECKIMAGE_TYPE") + 1] == "BACKGROUND,SEGMENTATION"
    assert argv[argv.index("-CHECKIMAGE_NAME") + 1] == (
        "/out/background-001-002.fits,/out/segmentation-001-002.fits"
    )
