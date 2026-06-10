"""UNIT TESTS FOR MODULE PACKAGE: NGMIX."""

from hypothesis import given
from hypothesis import strategies as st
import numpy as np
import numpy.testing as npt
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package.ngmix import Ngmix

rotated_ccds = st.integers(max_value=17) | st.sampled_from([36, 37])
unrotated_ccds = st.integers(min_value=18, max_value=40).filter(
    lambda ccd_nb: ccd_nb not in [36, 37]
)


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


def _empty_sqlite(path):
    db = SqliteDict(path)
    db.close()


def test_ngmix_accepts_optional_background_rms_vignet(tmp_path):
    """The optional seventh ngmix input is the BACKGROUND_RMS vignet sqlite."""
    names = ("gal", "bkg", "psf", "weight", "flag", "headers", "bkg_rms")
    sqlite_paths = [
        tmp_path / f"{name}.sqlite" for name in names
    ]
    for sqlite_path in sqlite_paths:
        _empty_sqlite(str(sqlite_path))
    input_paths = (
        ["tile_cat.fits"]
        + [str(path) for path in sqlite_paths[:5]]
        + [str(sqlite_paths[6])]
    )

    ngmix = Ngmix(
        input_paths,
        str(tmp_path),
        "-001-001",
        30.0,
        0.186,
        str(sqlite_paths[5]),
        _NullLogger(),
    )

    assert ngmix._vignet_cat.bkg_rms_vign_cat is not None
    ngmix._vignet_cat.close()


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


def test_background_rms_builds_per_pixel_inverse_variance():
    """BACKGROUND_RMS supplies the variance, while weight and flag supply masks."""
    from shapepipe.modules.ngmix_package.ngmix import prepare_ngmix_weights

    gal = np.ones((3, 3))
    weight = np.ones((3, 3))
    flag = np.zeros((3, 3))
    bkg_rms = np.array(
        [
            [1.0, 2.0, 4.0],
            [0.5, 0.0, np.nan],
            [3.0, 2.0, 1.0],
        ]
    )
    weight[2, 0] = 0
    flag[2, 1] = 1

    _, weight_map, noise_img = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(0), bkg_rms=bkg_rms
    )

    expected = np.array(
        [
            [1.0, 0.25, 0.0625],
            [4.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    npt.assert_allclose(weight_map, expected)
    assert noise_img.shape == gal.shape


def test_rescale_epoch_fluxes_scales_background_rms_like_image_counts():
    from astropy.io import fits
    from shapepipe.modules.ngmix_package.ngmix import rescale_epoch_fluxes

    header = fits.Header()
    header["FSCALE"] = 2.0

    gal_scaled, weight_scaled, bkg_rms_scaled = rescale_epoch_fluxes(
        np.ones((2, 2)),
        np.ones((2, 2)) * 8.0,
        header,
        np.ones((2, 2)) * 3.0,
    )

    npt.assert_allclose(gal_scaled, 2.0)
    npt.assert_allclose(weight_scaled, 2.0)
    npt.assert_allclose(bkg_rms_scaled, 6.0)


def test_spatially_varying_rms_survives_rescale_to_observation():
    """A spatially-varying RMS map reaches the Observation as exact per-pixel
    inverse variance through the full rescale -> prepare -> Observation chain.

    This is the #604 integration guard: the injected heteroscedastic truth
    ``1 / (Fscale * rms)**2`` must come back exactly (the RMS branch is
    deterministic and never looks at galaxy flux), with Megapipe-masked and
    flagged pixels zeroed.
    """
    from astropy.io import fits
    from shapepipe.modules.ngmix_package.ngmix import (
        make_ngmix_observation,
        rescale_epoch_fluxes,
    )
    from shapepipe.testing.simulate import make_data

    img_size = 51
    yy, xx = np.mgrid[:img_size, :img_size]
    rms = 1e-3 * (1.0 + 1.5 * xx / (img_size - 1) + 0.5 * yy / (img_size - 1))

    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(123),
        shear=(0.0, 0.0),
        noise=rms,
        n_epochs=1,
        img_size=img_size,
    )
    weights[0][:3, :3] = 0.0  # Megapipe-masked corner
    flags[0][-3:, -3:] = 8  # flagged corner

    header = fits.Header()
    header["FSCALE"] = 2.0
    gal_scaled, weight_scaled, rms_scaled = rescale_epoch_fluxes(
        gals[0], weights[0], header, bkg_rms=rms
    )

    obs = make_ngmix_observation(
        gal_scaled,
        weight_scaled,
        flags[0],
        psfs[0],
        jacobs[0],
        np.random.RandomState(0),
        bkg_rms=rms_scaled,
    )

    good = (weights[0] != 0) & (flags[0] == 0)
    expected = np.zeros_like(rms)
    expected[good] = 1.0 / (header["FSCALE"] * rms[good]) ** 2
    npt.assert_allclose(obs.weight, expected)


def test_constant_stamp_fallback_yields_finite_zero_weights():
    """A degenerate constant stamp (sigma_mad == 0) must not emit NaN weights."""
    from shapepipe.modules.ngmix_package.ngmix import prepare_ngmix_weights

    gal = np.ones((4, 4))
    weight = np.zeros((4, 4))  # fully masked
    flag = np.zeros((4, 4))

    with np.errstate(divide="raise", invalid="raise"):
        _, weight_map, _ = prepare_ngmix_weights(
            gal, weight, flag, np.random.RandomState(0)
        )

    npt.assert_array_equal(weight_map, 0.0)
