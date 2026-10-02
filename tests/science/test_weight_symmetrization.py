"""Weight symmetrization on the uberseg neighbour side.

Physics invariant: a one-sided zero-weight hole pulls the Gaussian fit, and
zeroing the weight on its three quarter-turn copies cancels the pull. A
neighbour footprint with no light 10 px from a galaxy with half-light radius
0.5" through a 0.7" PSF gives c1 = -7.8e-3 under uberseg, and 2e-6 with
``symmetrize_weights="interpolated_and_neighbours"``.

Why the default leaves the neighbour side unsymmetrized: a real neighbour
leaves light in the target's cell, which pulls the fit the other way, and
the one-sided hole partly offsets it. With a neighbour of a tenth of the
target's flux at 10 px, c1 = +1.1e-3 unsymmetrized and +1.2e-2 symmetrized.
If this second test turns red, the default deserves another look.
"""

import numpy as np

from tests.helpers.defect_response import defect_response

N = 51
SEEDS = range(6)
NO_DEFECT = np.zeros((N, N), dtype=bool)


def uberseg_response(flux_ratio, symmetrize):
    return defect_response(
        NO_DEFECT, hlr=0.5, psf=0.7, seeds=SEEDS,
        options=dict(blend_handling="uberseg", object_number=1,
                     dilate_neighbour=1, symmetrize_weights=symmetrize),
        neighbour=dict(offset=(0, 10), flux_ratio=flux_ratio, hlr=0.5,
                       seg_radius=4),
    )


def test_symmetrization_cancels_the_neighbour_side_hole():
    """Failure modes: the neighbour side is not symmetrized under
    "interpolated_and_neighbours", or the rotation is not about the stamp
    centre (the hole's pull survives). Positive control: the unsymmetrized
    hole gives |c1| > 5e-3."""
    unsymmetrized = uberseg_response(0.0, "interpolated")
    symmetrized = uberseg_response(0.0, "interpolated_and_neighbours")
    assert abs(unsymmetrized["c"][0]) > 5e-3, unsymmetrized
    assert np.abs(symmetrized["c"]).max() < 5e-5, symmetrized
    assert np.abs(symmetrized["m"]).max() < 0.01, symmetrized


def test_symmetrizing_the_neighbour_side_exposes_neighbour_light():
    """The default ("interpolated") leaves the uberseg neighbour side
    unsymmetrized because, with neighbour light in the image, symmetrizing
    it raises c1 several-fold."""
    unsymmetrized = uberseg_response(0.1, "interpolated")
    symmetrized = uberseg_response(0.1, "interpolated_and_neighbours")
    assert abs(unsymmetrized["c"][0]) < 3e-3, unsymmetrized
    assert symmetrized["c"][0] > 5e-3, symmetrized
