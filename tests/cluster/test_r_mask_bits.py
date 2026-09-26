"""The OR of the six r-ladder bits is the published 2D-cosmic-shear r mask.

Decision ``mask_default_cut`` (astra.yaml, catalogue_assembly) makes the
catalogue's default cut ``n1|n2|n4|n8|n64|n1024``, on the claim that this OR
reproduces ``mask_r_nside131072.hsp``, the mask the 2D analysis used. This
test checks that claim pixel-exactly on three whole coverage granules
(nside_cov 128, 2^20 nside-131072 pixels each) lying fully inside the v1.6.x
footprint. The granules sit in the interior on purpose: at the footprint
edge ``mask_r`` has no coverage and reads unmasked where the ladder masks, a
difference of order half a percent of area that carries no r data.

It reads the Feb-2025 r-band ladder beside the published mask. The DR6 ladder
that MASK_EXT_PATHS points at is not what this checks. Enforces contract
``mask-ext-ladder-columns``; skips when the maps are not readable.
"""

import os

import numpy as np
import pytest


pytestmark = [pytest.mark.slow, pytest.mark.candide]

MASK_DIR = "/n17data/UNIONS/WL/masks"
R_MASK = f"{MASK_DIR}/mask_r_nside131072.hsp"
R_BITS = ["n1", "n2", "n4", "n8", "n64", "n1024"]
# Whole nside_cov=128 granules fully inside coverage_v1.6.x.
INTERIOR_GRANULES = [32131, 29770, 43619]


def _paths():
    return [R_MASK] + [
        f"{MASK_DIR}/mask_r_nside131072_{bit}.hsp" for bit in R_BITS
    ]


@pytest.mark.parametrize("granule", INTERIOR_GRANULES)
def test_six_bit_or_reproduces_mask_r(granule):
    """OR of the six r-ladder bits equals mask_r on a whole granule."""
    unreadable = [p for p in _paths() if not os.access(p, os.R_OK)]
    if unreadable:
        pytest.skip(f"mask maps not readable: {unreadable}")
    hsp = pytest.importorskip("healsparse")

    pixels = granule * (1 << 20) + np.arange(1 << 20)

    def read(path):
        # A bit map with no coverage here reads its sentinel, False: unmasked.
        if not hsp.HealSparseCoverage.read(path).coverage_mask[granule]:
            assert path != R_MASK, f"mask_r has no coverage in {granule}"
            return np.zeros(pixels.size, dtype=bool)
        sub = hsp.HealSparseMap.read(path, pixels=[granule])
        return sub.get_values_pix(pixels, nest=True).astype(bool)

    mask_r = read(R_MASK)
    combined = np.logical_or.reduce([read(p) for p in _paths()[1:]])

    assert mask_r.any(), f"granule {granule}: mask_r masks nothing"
    n_diff = int(np.count_nonzero(mask_r != combined))
    assert n_diff == 0, (
        f"mask_default_cut: granule {granule}: OR of {R_BITS} differs from "
        f"mask_r on {n_diff} of {pixels.size} pixels"
    )
