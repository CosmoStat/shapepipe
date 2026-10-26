"""Shear recovery for the defects the central veto keeps.

Physics invariant: every defect the central veto keeps leaves both additive terms
|c1|, |c2| < 5e-4 and both diagonal multiplicative terms |m11|, |m22| < 1%,
from the full 2x2 response matrix. A filled defect biases m anisotropically
(m11 and m22 can differ tenfold), so a scalar m would hide it.
Interpolated defects (columns, full and finite 3-px bleeds, single pixels)
are checked at ``EPOCH_INTERPOLATED_DEFECT_RADIUS`` and two pixels beyond
it, on galaxies with half-light radius 0.3" and 0.5" through a 0.7" PSF,
round and with ellipticity (0.05, 0.02), and on a 0.7" galaxy through a
0.9" PSF. These pass only because the interpolated pixels' quarter turns
also lose their weight: an unsymmetrized column at 8 px gives
c1 = -1.3e-3. Defects too wide to interpolate (a 4-column cluster, a 5-px
bleed, edge bands) are noise-filled, and are checked at
``EPOCH_CENTRAL_DEFECT_RADIUS`` on the same three galaxies; the 0.7"
galaxy through the 0.9" PSF sets that radius (a 4-column cluster at 10 px
gives m11 = -6.6%). The cases sit at the radii themselves, so lowering
either below its calibrated value turns this red.

The bound is checked against the veto alone, not the masked-fraction cut:
the edge bands at the noise-fill radius cover about 25% of the stamp, which
``EPOCH_MASKED_FRACTION_CUT`` (10%) drops in production. That they stay
within the bound is the evidence that the fraction cut is not a bias
control, so it can be chosen for DES comparability and robustness alone.

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
def test_kept_defects_recover_shear_on_both_axes(kind, distance, hlr, fwhm,
                                                 psf_shear, tmp_path):
    """Failure modes: a veto radius is below its calibrated value; the
    interpolated pixels' quarter-turn orbit keeps its weight (a one-sided
    hole in the likelihood); the fill mask is symmetrized; a wide hole is
    interpolated; raw defect values leak into metacal; the check reads m11
    alone and misses m22, or c1 alone."""
    bad = geometry(kind, distance)
    masked = defect_mask(np.ones((N, N)), bad.astype(np.int32))
    np.testing.assert_array_equal(masked, bad)
    interpolated = interpolable_defects(masked)
    assert interpolated.any() == (kind in INTERPOLATED)
    assert not central_defect_vetoes(masked, interpolated)
    m, c, result = recover(bad, hlr, fwhm, psf_shear, tmp_path)
    assert m < 0.01, result
    assert c < 5e-4, result


def test_vetoed_bleed_breaks_the_bound(tmp_path):
    """Positive control: a 3-px bleed three pixels inside the
    interpolated-defect radius, on the 0.3" galaxy, gives |c| > 1e-3 and
    |m| > 1.5%. The veto drops it."""
    bad = geometry("bleed", RI - 3)
    assert central_defect_vetoes(bad, interpolable_defects(bad))
    m, c, result = recover(bad, 0.3, 0.7, ROUND, tmp_path)
    assert c > 1e-3, result
    assert m > 0.015, result
