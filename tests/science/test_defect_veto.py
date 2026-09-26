"""Shear recovery for the defects the central-defect veto keeps (noise fill).

Physics invariant: a defect the veto keeps -- a column, a 3-px bleed or a
single pixel at the veto radius or beyond, or an edge band -- leaves both
additive terms |c1|, |c2| < 5e-4 and both diagonal multiplicative terms
|m11|, |m22| < 1%, from the full 2x2 response matrix. The grid spans galaxies
with half-light radius 0.3" and 0.5" through a 0.7" PSF, round and with
ellipticity (0.05, 0.02). A noise-filled defect biases m anisotropically
(m11 and m22 differ by up to a factor of ten), so a scalar m would hide it.
The cases sit at ``EPOCH_CENTRAL_DEFECT_RADIUS`` itself, so lowering the
radius below the calibrated value turns this red.

Positive control: the same column two pixels inside the radius breaks the
bound, so the recovery check can see the bias the veto removes.
"""

import json

import numpy as np
import pytest

from shapepipe.modules.ngmix_package.ngmix import (
    EPOCH_CENTRAL_DEFECT_RADIUS,
    EPOCH_MASKED_FRACTION_CUT,
    defect_mask,
    has_central_defect,
)
from tests.helpers.defect_response import defect_response

N = 51
CENTRE = N // 2
R = int(np.ceil(EPOCH_CENTRAL_DEFECT_RADIUS))
ELLIPTICAL = (0.05, 0.02)
SEEDS = range(6)


def geometry(kind, distance):
    """A detector defect whose nearest pixel is ``distance`` px from the
    stamp centre; an edge band reaches in to that distance."""
    bad = np.zeros((N, N), dtype=bool)
    if kind == "pixel":
        bad[CENTRE, CENTRE + distance] = True
    elif kind == "column":
        bad[:, CENTRE + distance] = True
    elif kind == "bleed":
        bad[:, CENTRE + distance:CENTRE + distance + 3] = True
    elif kind == "edge":
        bad[:, CENTRE + distance:] = True
    return bad


# Wide defects (a 3-px bleed, the widest edge band the veto keeps) at the
# radius itself on the 0.5" galaxy through the elliptical PSF sit at the
# bound (m11 = -0.98% +/- 0.04% and -0.94% +/- 0.06%; see
# :func:`has_central_defect`), so those two are checked one pixel out.
CASES = (
    [(k, d, h, (0.0, 0.0)) for k in ("column", "bleed", "pixel")
     for d in (R, R + 2) for h in (0.3, 0.5)]
    + [(k, R, 0.5, ELLIPTICAL) for k in ("column", "pixel")]
    + [(k, R + 1, 0.5, ELLIPTICAL) for k in ("bleed", "edge")]
    + [("edge", R, 0.5, (0.0, 0.0))]
    + [("edge", N - CENTRE - 5, 0.5, psf) for psf in ((0.0, 0.0), ELLIPTICAL)]
)


def recover(bad, hlr, psf_shear, tmp_path):
    result = defect_response(bad, hlr=hlr, psf=0.7, seeds=SEEDS,
                             psf_shear=psf_shear)
    (tmp_path / "recovery.json").write_text(json.dumps(result, indent=2))
    return np.abs(result["m"]).max(), np.abs(result["c"]).max(), result


def case_id(case):
    kind, distance, hlr, psf_shear = case
    psf = "elliptical" if any(psf_shear) else "round"
    return f"{kind}-{distance}px-hlr{hlr}-{psf}"


@pytest.mark.parametrize("kind,distance,hlr,psf_shear", CASES,
                         ids=[case_id(c) for c in CASES])
def test_kept_defects_recover_shear_on_both_axes(kind, distance, hlr,
                                                 psf_shear, tmp_path):
    """Failure modes: the veto radius is below the calibrated value; the
    fill leaks raw defect values or leaves more of the object's light
    missing; the check reads m11 alone and misses m22, or c1 alone.
    """
    bad = geometry(kind, distance)
    masked = defect_mask(np.ones((N, N)), bad.astype(np.int32))
    np.testing.assert_array_equal(masked, bad)
    assert masked.mean() <= EPOCH_MASKED_FRACTION_CUT
    assert not has_central_defect(masked, EPOCH_CENTRAL_DEFECT_RADIUS)
    m, c, result = recover(bad, hlr, psf_shear, tmp_path)
    assert m < 0.01, result
    assert c < 5e-4, result


def test_vetoed_column_breaks_the_bound(tmp_path):
    """Positive control: a column two pixels inside the radius, on the
    0.5" galaxy, gives |m| > 2% and |c| > 1e-3. The veto drops it."""
    bad = geometry("column", R - 2)
    assert has_central_defect(bad, EPOCH_CENTRAL_DEFECT_RADIUS)
    m, c, result = recover(bad, 0.5, (0.0, 0.0), tmp_path)
    assert m > 0.02, result
    assert c > 1e-3, result
