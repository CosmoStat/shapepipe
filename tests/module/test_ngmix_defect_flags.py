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
    ccds = (10, 20, 21) if megacam_flip else (20, 21, 22)
    for i, (ccd, bit) in enumerate(zip(ccds, (4, 1024, 0))):
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
    # A footprint-rejected epoch and a fraction-cut epoch carry otherwise
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
    )


@pytest.mark.decision("shape_measurement.defect_flagging")
@pytest.mark.parametrize("flip", [True, False])
@pytest.mark.parametrize("kind, expected", [
    ("inside", (0, 0)), ("outside", (0, 0)),
    ("clean", (0, 0)), ("zero_weight", (0, 0)),
])
def test_footprint_flags_only_count_kept_epochs(kind, expected, flip):
    """Reject footprint hits and retain off-footprint defect diagnostics.

    Both flipped and unflipped epochs use their own registered footprint.
    """
    stamp = _run(kind, megacam_flip=flip)
    assert len(stamp.flags) == (1 if kind in ("inside", "zero_weight") else 3)
    diagnostics = stamp.defect_diagnostics()
    assert (diagnostics["imaflags_iso"],
            diagnostics["n_epoch_interp_iso"]) == expected
    assert diagnostics["n_epoch_interp"] == (2 if kind == "outside" else 0)


@pytest.mark.parametrize("seg_kind", ["missing", "empty"])
def test_unavailable_footprint_rejects_all_epochs(seg_kind):
    stamp = _run(seg_kind=seg_kind)
    assert not stamp.flags
    assert stamp.epoch_cuts == dict(
        considered=5, masked_fraction=1, footprint_defect=4,
    )


def test_clipped_footprint_is_tested_but_diagnostics_stay_unknown():
    stamp = _run(seg_kind="clipped")
    assert len(stamp.flags) == 1
    diagnostics = stamp.defect_diagnostics()
    assert diagnostics["imaflags_iso"] == -1
    assert diagnostics["n_epoch_interp_iso"] == -1


def test_footprint_record_preserves_bits_and_interpolated_epoch_counts():
    """The record API retains OR and sentinel semantics for downstream use."""
    from shapepipe.modules.ngmix_package.ngmix import Postage_stamp

    stamp = Postage_stamp()
    seg = np.zeros((9, 9), dtype=np.int32)
    seg[4, 4] = 17
    for bit in (4, 1024, 0):
        flag = np.zeros_like(seg)
        flag[4, 4] = bit
        stamp.record_footprint_defects(flag, flag != 0, seg, 17)
    assert stamp.imaflags_iso == 1028
    assert stamp.n_epoch_interp_iso == 2


@pytest.mark.parametrize("value", [np.nan, np.inf, -1, 1.5, float(2**63)])
def test_non_integral_or_unsafe_flags_raise(value):
    """Do not silently truncate float flags or overflow the integer OR."""
    from shapepipe.modules.ngmix_package.ngmix import Postage_stamp

    stamp = Postage_stamp()
    seg = np.zeros((5, 5), dtype=int)
    seg[2, 2] = 7
    flag = np.zeros(seg.shape)
    flag[2, 2] = value
    with pytest.raises(ValueError, match="footprint flags"):
        stamp.record_footprint_defects(flag, flag != 0, seg, 7)


def test_footprint_columns_survive_ngmix_fits(tmp_path):
    """Every shear extension keeps integer diagnostics, including -1."""
    from astropy.io import fits
    from tests.module.test_ngmix import _fake_metacal_result

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    inst._output_dir = str(tmp_path)
    inst._file_number_string = "-001-001"
    rows = []
    for i, (bits, count) in enumerate(((1028, 2), (-1, -1)), start=1):
        row = _fake_metacal_result(.18, .02, .09, .001)
        row.update(obj_id=i, imaflags_iso=bits, n_epoch_interp_iso=count)
        rows.append(row)
    inst.save_results(inst.compile_results(rows))
    with fits.open(inst.get_output_path(str(tmp_path))) as hdus:
        for hdu in hdus[1:]:
            assert hdu.data["imaflags_iso"].tolist() == [1028, -1]
            assert hdu.data["n_epoch_interp_iso"].tolist() == [2, -1]
            assert hdu.data["imaflags_iso"].dtype.kind == "i"
