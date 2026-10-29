"""Per-fixture shear recovery for interpolated and noise-filled defects.

These fixtures test the fill and weighting at specified distances (7/9 px
for interpolation, 13 px for noise), not the D0 epoch-selection rule.
The full response matrix must give |m11|, |m22| < 1% and |c1|, |c2| < 5e-4.
A close 3-px bleed is a positive control that breaks those bounds.
D0's justification concerns prevalence-weighted catalogue bias with
multi-epoch dilution, not a shear-recovery guarantee for every kept defect.
"""

import json

import numpy as np
import pytest

from shapepipe.modules.ngmix_package.defect_interpolation import (
    interpolable_defects,
)
from shapepipe.modules.ngmix_package.ngmix import (
    defect_mask,
)
from tests.helpers.defect_response import defect_response

N = 51
CENTRE = N // 2
RI, RN = 7, 13  # fixture distances, independent of epoch selection
ROUND = (0.0, 0.0)
ELLIPTICAL = (0.05, 0.02)
SEEDS = range(6)


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
    elif kind == "cluster4":
        bad[:, near:near + 4] = True
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
    + [("edge", RN, 0.5, 0.7, ELLIPTICAL)]
    + [(k, RN, 0.7, 0.9, psf) for k in ("cluster4", "edge")
       for psf in (ROUND, ELLIPTICAL)]
    + [("edge", N - CENTRE - 5, 0.5, 0.7, psf) for psf in (ROUND, ELLIPTICAL)]
)


def case_id(case):
    kind, distance, hlr, fwhm, psf_shear = case
    psf = "elliptical" if any(psf_shear) else "round"
    return f"{kind}-{distance}px-hlr{hlr}-psf{fwhm}-{psf}"


def recover(bad, hlr, fwhm, psf_shear, tmp_path):
    result = defect_response(bad, hlr=hlr, psf=fwhm, seeds=SEEDS,
                             psf_shear=psf_shear)
    (tmp_path / "recovery.json").write_text(json.dumps(result, indent=2))
    return np.abs(result["m"]).max(), np.abs(result["c"]).max(), result


@pytest.mark.parametrize("kind,distance,hlr,fwhm,psf_shear", CASES,
                         ids=[case_id(c) for c in CASES])
def test_fill_fixtures_recover_shear_on_both_axes(kind, distance, hlr, fwhm,
                                                 psf_shear, tmp_path):
    """Failure modes: the interpolated pixels' quarter-turn orbit keeps
    its weight (a one-sided
    hole in the likelihood); the fill mask is symmetrized; a wide hole is
    interpolated; raw defect values leak into metacal; the check reads m11
    alone and misses m22, or c1 alone."""
    bad = geometry(kind, distance)
    masked = defect_mask(np.ones((N, N)), bad.astype(np.int32))
    np.testing.assert_array_equal(masked, bad)
    interpolated = interpolable_defects(masked)
    assert interpolated.any() == (kind in INTERPOLATED)
    m, c, result = recover(bad, hlr, fwhm, psf_shear, tmp_path)
    assert m < 0.01, result
    assert c < 5e-4, result


def test_close_bleed_breaks_the_bound(tmp_path):
    """Positive control: a 3-px bleed at 4 px on the 0.3" galaxy gives
    |c| > 1e-3 and |m| > 1.5%."""
    bad = geometry("bleed", RI - 3)
    m, c, result = recover(bad, 0.3, 0.7, ROUND, tmp_path)
    assert c > 1e-3, result
    assert m > 0.015, result
