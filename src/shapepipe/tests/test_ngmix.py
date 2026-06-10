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


def _fake_metacal_result(T, T_err, T_psf, T_psf_err):
    """Build one minimal metacal result as produced by the process loop.

    The PSF size entries mirror the conversion done at the source:
    ``r50_PSFo = sqrt(2 ln 2) * sigma`` with ``sigma = sqrt(T_psf / 2)``.
    """
    from shapepipe.modules.ngmix_package.ngmix import SIGMA_TO_R50

    sigma_psf = np.sqrt(max(T_psf, 0) / 2)
    per_type = {
        "nfev": 1,
        "g": [0.01, -0.02],
        "g_cov": np.diag([1e-4, 1e-4]),
        "T": T,
        "T_err": T_err,
        "flux": 100.0,
        "flux_err": 1.0,
        "s2n": 50.0,
        "flags": 0,
    }
    res = {
        "obj_id": 1,
        "n_epoch_model": 1,
        "moments_fail": 0,
        "g_PSFo": [0.001, 0.002],
        "g_err_PSFo": [1e-5, 1e-5],
        "T_PSFo": T_psf,
        "T_err_PSFo": T_psf_err,
        "r50_PSFo": SIGMA_TO_R50 * sigma_psf,
        "r50_err_PSFo": (
            SIGMA_TO_R50 * T_psf_err / (4 * sigma_psf)
            if sigma_psf > 0
            else np.nan
        ),
        "mcal_flags": 0,
    }
    res.update(
        {name: dict(per_type) for name in ("1m", "1p", "2m", "2p", "noshear")}
    )
    return res


def test_compile_results_size_columns_are_half_light_radii():
    """Every r50 column is a true half-light radius, on both sides.

    Galaxy ``r50 = sqrt(ln 2 * T)`` (not the raw area ``pars[4]``), PSF
    ``r50psf = sqrt(2 ln 2) * sigma_psf`` (not bare sigma), and the
    redundant ``*_psfo_ngmix`` size duplicates are gone.
    """
    from shapepipe.modules.ngmix_package.ngmix import Ngmix, SIGMA_TO_R50

    T, T_err = 0.18, 0.02
    T_psf, T_psf_err = 0.09, 0.001

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    out = inst.compile_results([_fake_metacal_result(T, T_err, T_psf, T_psf_err)])

    noshear = out["noshear"]

    # galaxy: r50 = sqrt(ln2 * T) = 1.17741 * sqrt(T / 2), error dT * r50/(2T)
    r50_expected = np.sqrt(np.log(2) * T)
    npt.assert_allclose(noshear["r50"], [r50_expected])
    npt.assert_allclose(noshear["r50_err"], [r50_expected * T_err / (2 * T)])

    # PSF: r50psf = sqrt(2 ln2) * sigma, sigma = sqrt(T_psf / 2)
    sigma_psf = np.sqrt(T_psf / 2)
    npt.assert_allclose(noshear["r50psf"], [SIGMA_TO_R50 * sigma_psf])
    npt.assert_allclose(
        noshear["r50psf_err"], [SIGMA_TO_R50 * T_psf_err / (4 * sigma_psf)]
    )

    # galaxy/PSF r50 are now commensurable: same convention on both sides
    npt.assert_allclose(
        np.array(noshear["r50"]) / np.array(noshear["r50psf"]),
        [np.sqrt(T / T_psf)],
    )

    # areas pass through untouched, with the err asymmetry fixed
    npt.assert_allclose(noshear["Tpsf"], [T_psf])
    npt.assert_allclose(noshear["Tpsf_err"], [T_psf_err])

    # the *_psfo_ngmix size duplicates are retired
    for retired in (
        "T_psfo_ngmix",
        "T_err_psfo_ngmix",
        "r50_psfo_ngmix",
        "r50_err_psfo_ngmix",
    ):
        assert retired not in noshear


def test_compile_results_nonpositive_T_gives_nan_r50():
    """A non-positive fitted area cannot yield a real half-light radius."""
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    out = inst.compile_results([_fake_metacal_result(-0.05, 0.02, 0.09, 0.001)])

    assert np.isnan(out["noshear"]["r50"]).all()
    assert np.isnan(out["noshear"]["r50_err"]).all()
