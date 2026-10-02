"""Shear recovery for the defects kept under ``DEFECT_FILL = interpolate``.

Physics invariant: every defect the veto keeps leaves both additive terms
|c1|, |c2| < 5e-4 and both diagonal multiplicative terms |m11|, |m22| < 1%,
from the full 2x2 response matrix. Interpolated defects (columns, full and
finite 3-px bleeds, single pixels) are checked at
``EPOCH_INTERPOLATED_DEFECT_RADIUS`` and two pixels beyond it, on galaxies
with half-light radius 0.3" and 0.5" through a 0.7" PSF, round and with
ellipticity (0.05, 0.02), and on a 0.7" galaxy through a 0.9" PSF. Defects
too wide to interpolate are noise-filled, and are checked at
``EPOCH_CENTRAL_DEFECT_RADIUS``. The cases sit at the radii themselves, so
lowering either below its calibrated value turns this red.

Positive control: a 3-px bleed three pixels inside the interpolated-defect
radius breaks the bound.
"""

import json

import numpy as np
import pytest

from shapepipe.modules.ngmix_package.defect_interpolation import (
    interpolable_defects,
)
from shapepipe.modules.ngmix_package.ngmix import (
    EPOCH_CENTRAL_DEFECT_RADIUS,
    EPOCH_INTERPOLATED_DEFECT_RADIUS,
    EPOCH_MASKED_FRACTION_CUT,
    central_defect_vetoes,
    defect_mask,
)
from tests.helpers.defect_response import defect_response

N = 51
CENTRE = N // 2
RI = int(np.ceil(EPOCH_INTERPOLATED_DEFECT_RADIUS))
RN = int(np.ceil(EPOCH_CENTRAL_DEFECT_RADIUS))
ROUND = (0.0, 0.0)
ELLIPTICAL = (0.05, 0.02)
SEEDS = range(6)
INTERPOLATE = {"defect_fill": "interpolate"}


def geometry(kind, distance):
    """A detector defect whose nearest pixel is ``distance`` px from the
    stamp centre."""
    bad = np.zeros((N, N), dtype=bool)
    near = CENTRE + distance
    if kind == "pixel":
        bad[CENTRE, near] = True
    elif kind == "column":
        bad[:, near] = True
    elif kind == "bleed":
        bad[:, near:near + 3] = True
    elif kind == "finite_bleed":
        bad[CENTRE - 5:CENTRE + 6, near:near + 3] = True
    elif kind == "wide_bleed":
        bad[:, near:near + 5] = True
    elif kind == "edge":
        bad[:, near:] = True
    return bad


INTERPOLATED = ("column", "bleed", "finite_bleed", "pixel")
CASES = (
    [(k, d, h, 0.7, ROUND) for k in INTERPOLATED
     for d in (RI, RI + 2) for h in (0.3, 0.5)]
    + [(k, RI, 0.5, 0.7, ELLIPTICAL) for k in INTERPOLATED]
    + [(k, RI, 0.7, 0.9, psf) for k in ("bleed", "finite_bleed")
       for psf in (ROUND, ELLIPTICAL)]
    + [(k, RN, h, 0.7, ROUND) for k in ("wide_bleed", "edge")
       for h in (0.3, 0.5)]
)


def case_id(case):
    kind, distance, hlr, fwhm, psf_shear = case
    psf = "elliptical" if any(psf_shear) else "round"
    return f"{kind}-{distance}px-hlr{hlr}-psf{fwhm}-{psf}"


def recover(bad, hlr, fwhm, psf_shear, tmp_path):
    result = defect_response(bad, hlr=hlr, psf=fwhm, seeds=SEEDS,
                             psf_shear=psf_shear, options=INTERPOLATE)
    (tmp_path / "recovery.json").write_text(json.dumps(result, indent=2))
    return np.abs(result["m"]).max(), np.abs(result["c"]).max(), result


@pytest.mark.parametrize("kind,distance,hlr,fwhm,psf_shear", CASES,
                         ids=[case_id(c) for c in CASES])
def test_kept_defects_recover_shear_on_both_axes(kind, distance, hlr, fwhm,
                                                 psf_shear, tmp_path):
    """Failure modes: a veto radius is below its calibrated value; the
    interpolated pixels' quarter-turn orbit keeps its weight (a one-sided
    hole in the likelihood); the fill mask is symmetrized; a wide hole is
    interpolated; raw defect values leak into metacal."""
    bad = geometry(kind, distance)
    masked = defect_mask(np.ones((N, N)), bad.astype(np.int32))
    interpolated = interpolable_defects(masked)
    assert interpolated.any() == (kind in INTERPOLATED)
    assert masked.mean() <= EPOCH_MASKED_FRACTION_CUT
    assert not central_defect_vetoes(masked, EPOCH_CENTRAL_DEFECT_RADIUS,
                                     "interpolate")
    m, c, result = recover(bad, hlr, fwhm, psf_shear, tmp_path)
    assert m < 0.01, result
    assert c < 5e-4, result


def test_vetoed_bleed_breaks_the_bound(tmp_path):
    """Positive control: a 3-px bleed three pixels inside the
    interpolated-defect radius, on the 0.3" galaxy, gives |c| > 1e-3 and
    |m| > 1.5%. The veto drops it."""
    bad = geometry("bleed", RI - 3)
    assert central_defect_vetoes(bad, EPOCH_CENTRAL_DEFECT_RADIUS,
                                 "interpolate")
    m, c, result = recover(bad, 0.3, 0.7, ROUND, tmp_path)
    assert c > 1e-3, result
    assert m > 0.015, result
