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


def test_compile_results_nan_fills_failed_fit_types():
    """A failed fit type must be recorded with NaNs, not crash the tile.

    ngmix 2.x ``run_fitter`` does not raise on failure: after ``ntry`` it
    returns a result with ``flags != 0`` that carries none of the
    measurement keys (g, g_cov, T, T_err, flux, flux_err, s2n).
    ``compile_results`` previously indexed those keys directly, so a single
    failed object KeyError-crashed the whole tile at save time, hours in.
    """
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    res = _fake_metacal_result(0.18, 0.02, 0.09, 0.001)
    res["1p"] = {"flags": 0x8, "nfev": 5}  # failed fit: only flags/nfev
    res["moments_fail"] = 1

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    out = inst.compile_results([res])

    failed = out["1p"]
    for col in (
        "g1", "g2", "g1_err", "g2_err", "T", "T_err", "flux", "flux_err",
        "s2n", "mag", "mag_err", "r50", "r50_err",
    ):
        assert np.isnan(failed[col]).all(), col
    assert failed["flags"] == [0x8]

    # successful types are untouched
    npt.assert_allclose(out["noshear"]["flux"], [100.0])
    assert out["noshear"]["flags"] == [0]


def test_compile_results_nonpositive_T_gives_nan_r50():
    """A non-positive fitted area cannot yield a real half-light radius."""
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    out = inst.compile_results([_fake_metacal_result(-0.05, 0.02, 0.09, 0.001)])

    assert np.isnan(out["noshear"]["r50"]).all()
    assert np.isnan(out["noshear"]["r50_err"]).all()


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
