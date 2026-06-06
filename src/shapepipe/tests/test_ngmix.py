"""UNIT TESTS FOR MODULE PACKAGE: NGMIX."""

from hypothesis import given
from hypothesis import strategies as st
import numpy as np
import numpy.testing as npt

from shapepipe.modules.ngmix_package.ngmix import Ngmix

rotated_ccds = st.integers(max_value=17) | st.sampled_from([36, 37])
unrotated_ccds = st.integers(min_value=18, max_value=40).filter(
    lambda ccd_nb: ccd_nb not in [36, 37]
)


@given(
    st.lists(
        st.lists(
            st.floats(allow_nan=False, allow_infinity=False), min_size=1, max_size=8
        ),
        min_size=1,
        max_size=8,
    ).filter(lambda rows: len({len(row) for row in rows}) == 1),
    rotated_ccds,
)
def test_megacam_flip_preserves_shape_and_is_involution(rows, ccd_nb):

    vign = np.array(rows, dtype=float)
    flipped = Ngmix.MegaCamFlip(vign, ccd_nb)

    assert flipped.shape == vign.shape
    npt.assert_allclose(Ngmix.MegaCamFlip(flipped, ccd_nb), vign)


@given(
    st.lists(
        st.lists(
            st.floats(allow_nan=False, allow_infinity=False), min_size=1, max_size=8
        ),
        min_size=1,
        max_size=8,
    ).filter(lambda rows: len({len(row) for row in rows}) == 1),
    unrotated_ccds,
)
def test_megacam_flip_leaves_unrotated_ccds_unchanged(rows, ccd_nb):

    vign = np.array(rows, dtype=float)

    npt.assert_allclose(Ngmix.MegaCamFlip(vign, ccd_nb), vign)


def _metacal_noshear_g(seed):
    """Run do_ngmix_metacal on a fixed simulated stamp with RandomState(seed).

    The simulated data is drawn from a separate, fixed seed so the only
    seed-dependent randomness in play is the noise/metacal RNG inside the
    fitter — exactly the channel that must be reproducible.
    """
    from shapepipe.modules.ngmix_package.ngmix import (
        Postage_stamp,
        do_ngmix_metacal,
        get_prior,
    )
    from shapepipe.testing.simulate import make_data

    rng = np.random.RandomState(seed)
    prior = get_prior(0.1857, rng)
    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(123),
        shear=(0.02, 0.0),
        noise=1e-4,
        n_epochs=2,
        img_size=51,
    )
    stamp = Postage_stamp(bkg_sub=False, megacam_flip=False)
    stamp.gals, stamp.psfs, stamp.weights, stamp.flags, stamp.jacobs = (
        gals,
        psfs,
        weights,
        flags,
        jacobs,
    )
    res, _ = do_ngmix_metacal(stamp, prior, 1.0, rng)
    return np.asarray(res["noshear"]["g"])


def test_metacal_is_reproducible_with_fixed_seed():
    """Same seed -> identical metacal shear.

    The module seeds ``self._rng = RandomState(seed)`` per tile precisely so a
    rerun reproduces. This guards the noise-image and masked-pixel draws in
    ``prepare_ngmix_weights`` against silently falling back to the unseeded
    global ``numpy.random`` state, which would make shear estimates
    irreproducible from one run to the next.
    """
    npt.assert_array_equal(_metacal_noshear_g(42), _metacal_noshear_g(42))


def test_weight_map_recovers_injected_inverse_variance():
    """A supplied inverse-variance map must not be renormalized twice."""
    from shapepipe.modules.ngmix_package.ngmix import prepare_ngmix_weights
    from shapepipe.testing.simulate import make_data

    noise = 1e-3
    gals, _, _, weights, flags, _ = make_data(
        rng=np.random.RandomState(123),
        shear=(0.0, 0.0),
        noise=noise,
        n_epochs=1,
        img_size=201,
    )
    _, weight_map, _ = prepare_ngmix_weights(
        gals[0], weights[0], flags[0], np.random.RandomState(0)
    )

    recovered = np.median(weight_map[weight_map > 0])
    npt.assert_allclose(recovered, 1.0 / noise**2, rtol=0.15)
