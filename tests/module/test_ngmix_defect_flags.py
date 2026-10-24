"""Adaptive defect flags on the coadd object's footprint, over kept epochs."""

import numpy as np
import pytest

from shapepipe.modules.ngmix_package.ngmix import Ngmix, prepare_postage_stamps
from tests.module.test_ngmix_defect_fill import N_STAMP, _fake_inputs


def _run(kind="inside", seg_kind="valid", megacam_flip=True):
    seg = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    seg[23:28, 23:37] = 1  # reaches 11 px right of centre, no stamp edge
    if seg_kind == "empty":
        seg[:] = 0
    if seg_kind == "clipped":
        seg[25, 36:] = 1
    epochs = {}
    for i, (ccd, bit) in enumerate(((10, 4), (20, 1024), (21, 0))):
        flag = np.zeros(seg.shape, dtype=np.float32)
        weight = np.ones(seg.shape)
        pixel = (25, 35) if kind != "outside" else (35, 35)
        if kind != "clean":
            flag[pixel] = bit
            if kind == "zero_weight" and bit:
                flag[pixel] = 0
                weight[pixel] = 0
        if megacam_flip:
            flag = Ngmix.MegaCamFlip(flag, ccd)
            weight = Ngmix.MegaCamFlip(weight, ccd)
        epochs[f"210000{i}-{ccd}"] = (flag, weight)
    # A centrally vetoed epoch and a fraction-cut epoch carry otherwise
    # unique bits on the footprint. Neither may contribute.
    for i, fraction in enumerate((False, True), start=3):
        flag = np.zeros(seg.shape, dtype=np.float32)
        flag[25, 35] = 64
        if fraction:
            flag[:10] = 64
        else:
            flag[25, 25] = 64
        epochs[f"210000{i}-20"] = (flag, np.ones(seg.shape))
    vignet, tile, psf, gal = _fake_inputs(epochs)
    tile.seg = None if seg_kind == "missing" else seg[None]
    return prepare_postage_stamps(
        vignet, 1, 0, tile, bkg_sub=False, psf_obj=psf, gal_obj=gal,
        megacam_flip=megacam_flip,
    )


@pytest.mark.decision("shape_measurement.defect_flagging")
@pytest.mark.parametrize("flip", [True, False])
@pytest.mark.parametrize("kind, expected", [
    ("inside", (1028, 2)), ("outside", (0, 0)),
    ("clean", (0, 0)), ("zero_weight", (0, 2)),
])
def test_footprint_flags_only_count_kept_epochs(kind, expected, flip):
    """Keep bit identities across epochs, ignore outside/vetoed defects.

    The final clean epoch must not overwrite the OR of earlier epochs.
    A zero-weight pixel has no instrument bit but counts if interpolated.
    """
    stamp = _run(kind, megacam_flip=flip)
    assert len(stamp.flags) == 3
    diagnostics = stamp.defect_diagnostics()
    assert (diagnostics["imaflags_iso"],
            diagnostics["n_epoch_interp_iso"]) == expected
    assert diagnostics["n_epoch_interp"] == (0 if kind == "clean" else 2)


@pytest.mark.parametrize("seg_kind", ["missing", "empty", "clipped"])
def test_unknown_or_incomplete_footprint_is_not_clean(seg_kind):
    diagnostics = _run(seg_kind=seg_kind).defect_diagnostics()
    assert diagnostics["imaflags_iso"] == -1
    assert diagnostics["n_epoch_interp_iso"] == -1
