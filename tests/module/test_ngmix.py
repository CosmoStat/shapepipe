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
    res, _, _ = do_ngmix_metacal(stamp, prior, 1.0, rng)
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


# The two PSF families carry distinct ellipticity AND size (shapepipe#749):
# the original image PSF (-> *_psf_orig) and the metacal reconvolution kernel
# (-> *_psf_reconv) are independent fits, so a regression to the old aliasing
# (both from one source) would surface here. The original PSF is elliptical
# and smaller than the round, enlarged reconvolution kernel.
RECONV_PSF_G = [0.001, 0.002]
RECONV_PSF_G_ERR = [1e-5, 1e-5]
ORIG_PSF_G = [0.004, -0.003]
ORIG_PSF_G_ERR = [2e-5, 3e-5]
# original-PSF size, deliberately != the reconvolution T_psf the builder is
# handed, so the size un-aliasing is exercised end to end.
ORIG_PSF_T = 0.07
ORIG_PSF_T_ERR = 0.0008


def _psf_r50_pair(T_psf, T_psf_err, sigma_to_r50):
    """``(r50, r50_err)`` for a PSF size, mirroring the source conversion."""
    sigma = np.sqrt(max(T_psf, 0) / 2)
    return (
        sigma_to_r50 * sigma,
        sigma_to_r50 * T_psf_err / (4 * sigma) if sigma > 0 else np.nan,
    )


def _fake_metacal_result(T, T_err, T_psf, T_psf_err):
    """Build one minimal metacal result as produced by the process loop.

    Both PSF families are filled with self-named scalar keys, as the process
    loop back-fills them: ``*_psf_orig`` (original image PSF) and
    ``*_psf_reconv`` (metacal reconvolution kernel). The reconvolution-kernel
    size is the builder's ``T_psf`` argument; the original-PSF size is the
    fixed ``ORIG_PSF_T``, deliberately different so the un-aliasing is tested.
    PSF r50 = ``sqrt(2 ln 2) * sigma`` with ``sigma = sqrt(T / 2)``.
    """
    from shapepipe.modules.ngmix_package.ngmix import SIGMA_TO_R50

    r50_reconv, r50_err_reconv = _psf_r50_pair(T_psf, T_psf_err, SIGMA_TO_R50)
    r50_orig, r50_err_orig = _psf_r50_pair(
        ORIG_PSF_T, ORIG_PSF_T_ERR, SIGMA_TO_R50
    )
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
        "mcal_types_fail": 0,
        # original image PSF (psfex/mccd) family
        "g1_psf_orig": ORIG_PSF_G[0],
        "g2_psf_orig": ORIG_PSF_G[1],
        "g1_err_psf_orig": ORIG_PSF_G_ERR[0],
        "g2_err_psf_orig": ORIG_PSF_G_ERR[1],
        "T_psf_orig": ORIG_PSF_T,
        "T_err_psf_orig": ORIG_PSF_T_ERR,
        "r50_psf_orig": r50_orig,
        "r50_err_psf_orig": r50_err_orig,
        # metacal reconvolution-kernel family
        "g1_psf_reconv": RECONV_PSF_G[0],
        "g2_psf_reconv": RECONV_PSF_G[1],
        "g1_err_psf_reconv": RECONV_PSF_G_ERR[0],
        "g2_err_psf_reconv": RECONV_PSF_G_ERR[1],
        "T_psf_reconv": T_psf,
        "T_err_psf_reconv": T_psf_err,
        "r50_psf_reconv": r50_reconv,
        "r50_err_psf_reconv": r50_err_reconv,
        "mcal_flags": 0,
    }
    res.update(
        {name: dict(per_type) for name in ("1m", "1p", "2m", "2p", "noshear")}
    )
    return res


def test_compile_results_size_columns_are_half_light_radii():
    """Every r50 column is a true half-light radius, on every side.

    Galaxy ``r50 = sqrt(ln 2 * T)`` (not the raw area ``pars[4]``) and both
    PSF families ``r50 = sqrt(2 ln 2) * sigma`` (not bare sigma), in the same
    convention so galaxy / PSF radii are commensurable. The reconvolution
    kernel (``*_psf_reconv``) and the original image PSF (``*_psf_orig``) are
    independent fits with their own, different sizes (shapepipe#749).
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

    # reconvolution kernel: r50 = sqrt(2 ln2) * sigma, sigma = sqrt(T_psf / 2)
    sigma_reconv = np.sqrt(T_psf / 2)
    npt.assert_allclose(noshear["r50_psf_reconv"], [SIGMA_TO_R50 * sigma_reconv])
    npt.assert_allclose(
        noshear["r50_err_psf_reconv"],
        [SIGMA_TO_R50 * T_psf_err / (4 * sigma_reconv)],
    )
    npt.assert_allclose(noshear["T_psf_reconv"], [T_psf])
    npt.assert_allclose(noshear["T_err_psf_reconv"], [T_psf_err])

    # galaxy/reconv-PSF r50 are commensurable: same convention on both sides
    npt.assert_allclose(
        np.array(noshear["r50"]) / np.array(noshear["r50_psf_reconv"]),
        [np.sqrt(T / T_psf)],
    )

    # original image PSF carries its OWN size, un-aliased from the kernel
    sigma_orig = np.sqrt(ORIG_PSF_T / 2)
    npt.assert_allclose(noshear["T_psf_orig"], [ORIG_PSF_T])
    npt.assert_allclose(noshear["T_err_psf_orig"], [ORIG_PSF_T_ERR])
    npt.assert_allclose(noshear["r50_psf_orig"], [SIGMA_TO_R50 * sigma_orig])
    npt.assert_allclose(
        noshear["r50_err_psf_orig"],
        [SIGMA_TO_R50 * ORIG_PSF_T_ERR / (4 * sigma_orig)],
    )
    assert noshear["T_psf_orig"] != noshear["T_psf_reconv"]
    assert noshear["r50_psf_orig"] != noshear["r50_psf_reconv"]


def test_compile_results_psf_families_are_unaliased():
    """The two PSF-ellipticity families document *different* PSFs (#749).

    ``*_psf_reconv`` carries the metacal RECONVOLUTION kernel; ``*_psf_orig``
    carries the ORIGINAL image PSF, fit independently. Before the fix both
    came from one ``average_multiepoch_psf`` result and so were byte-
    identical; here they must differ, each tracing its own source. The
    companion size un-aliasing is checked in
    ``test_compile_results_size_columns_are_half_light_radii``.
    """
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    noshear = inst.compile_results(
        [_fake_metacal_result(0.18, 0.02, 0.09, 0.001)]
    )["noshear"]

    # reconvolution kernel
    npt.assert_allclose(noshear["g1_psf_reconv"], [RECONV_PSF_G[0]])
    npt.assert_allclose(noshear["g2_psf_reconv"], [RECONV_PSF_G[1]])
    npt.assert_allclose(noshear["g1_err_psf_reconv"], [RECONV_PSF_G_ERR[0]])
    npt.assert_allclose(noshear["g2_err_psf_reconv"], [RECONV_PSF_G_ERR[1]])
    # original image PSF
    npt.assert_allclose(noshear["g1_psf_orig"], [ORIG_PSF_G[0]])
    npt.assert_allclose(noshear["g2_psf_orig"], [ORIG_PSF_G[1]])
    npt.assert_allclose(noshear["g1_err_psf_orig"], [ORIG_PSF_G_ERR[0]])
    npt.assert_allclose(noshear["g2_err_psf_orig"], [ORIG_PSF_G_ERR[1]])
    # the un-aliasing: the two families are no longer the same value
    assert noshear["g1_psf_orig"] != noshear["g1_psf_reconv"]
    assert noshear["g2_psf_orig"] != noshear["g2_psf_reconv"]


def test_compile_results_nan_fills_failed_fit_types():
    """A failed fit type must be recorded with NaNs, not crash the tile.

    ngmix 2.x ``run_fitter`` does not raise on failure: after ``ntry`` it
    returns a result with ``flags != 0`` that carries none of the
    measurement keys (g, g_cov, T, T_err, flux, flux_err, s2n).
    ``compile_results`` previously indexed those keys directly, so a single
    failed object KeyError-crashed the whole tile at save time, hours in.

    The failed object must also carry nonzero ``mcal_flags`` (the OR of
    the per-type fit flags, computed as the process loop does) so the
    downstream NGMIX_MCAL_FLAGS quality cut can see it.
    """
    from shapepipe.modules.ngmix_package.ngmix import Ngmix, get_mcal_flags

    res = _fake_metacal_result(0.18, 0.02, 0.09, 0.001)
    res["1p"] = {"flags": 0x8, "nfev": 5}  # failed fit: only flags/nfev
    res["mcal_types_fail"] = 1
    res["mcal_flags"] = get_mcal_flags(res)

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
    assert failed["mcal_flags"] == [0x8] and failed["mcal_flags"][0] != 0

    # successful types are untouched, but share the object's mcal_flags
    npt.assert_allclose(out["noshear"]["flux"], [100.0])
    assert out["noshear"]["flags"] == [0]
    assert out["noshear"]["mcal_flags"] == [0x8]


def test_get_mcal_flags_ors_per_type_fit_flags():
    """mcal_flags is the bitwise OR of all per-type fit flags (v1 contract).

    The rewrite hard-coded ``res['mcal_flags'] = 0``, making the
    NGMIX_MCAL_FLAGS column constant-zero so any mcal_flags == 0 quality
    cut passed everything.
    """
    from shapepipe.modules.ngmix_package.ngmix import get_mcal_flags

    res = {name: {"flags": 0} for name in ("noshear", "1p", "1m", "2p", "2m")}
    assert get_mcal_flags(res) == 0

    res["1p"]["flags"] = 0x8
    res["2m"]["flags"] = 0x2
    assert get_mcal_flags(res) == 0xA


def test_compile_results_nonpositive_T_gives_nan_r50():
    """A non-positive fitted area cannot yield a real half-light radius."""
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    out = inst.compile_results([_fake_metacal_result(-0.05, 0.02, 0.09, 0.001)])

    assert np.isnan(out["noshear"]["r50"]).all()
    assert np.isnan(out["noshear"]["r50_err"]).all()


def test_average_multiepoch_psf_skips_failed_psf_epochs():
    """A failed-PSF epoch must be skipped, not KeyError the whole object.

    With ``ignore_failed_psf=True`` the bootstrapper drops failed-PSF
    epochs from the galaxy fit but leaves them in the returned obsdict,
    where their ``psf.meta['result']`` carries only flags/pars (no T/g).
    Averages must come from the surviving epochs only, and ``n_epoch``
    must count those survivors (v1 counted PSF-fit-successful epochs).
    """
    from types import SimpleNamespace
    from shapepipe.modules.ngmix_package.ngmix import average_multiepoch_psf

    def epoch(result, weight_value):
        return SimpleNamespace(
            psf=SimpleNamespace(meta={"result": result}),
            weight=np.full((2, 2), weight_value),
        )

    good = {
        "flags": 0,
        "T": 0.2,
        "T_err": 0.01,
        "g": np.array([0.01, 0.02]),
        "g_err": np.array([1e-3, 2e-3]),
    }
    failed = {"flags": 3, "nfev": 5, "pars": np.zeros(6)}  # no T/g keys

    psf_res = average_multiepoch_psf(
        {"noshear": [epoch(good, 1.0), epoch(failed, 4.0)]}
    )

    npt.assert_allclose(psf_res["T_psf"], 0.2)
    npt.assert_allclose(psf_res["T_psf_err"], 0.01)
    npt.assert_allclose(psf_res["g_psf"], [0.01, 0.02])
    npt.assert_allclose(psf_res["g_psf_err"], [1e-3, 2e-3])
    assert psf_res["n_epoch"] == 1


def test_average_original_psf_fits_each_gal_psf_with_galaxy_weight():
    """The original-PSF average fits ``gal_obs.psf`` and weights by gal S/N.

    Restores the independent original-image-PSF fit (shapepipe#749): each
    epoch's psfex/mccd PSF stamp (``gal_obs.psf``) is fit by the shared
    ``psf_runner``, then weight-averaged by the *galaxy* inverse-variance
    weight (matching :func:`average_multiepoch_psf`). A pre-set
    ``psf.meta['result']`` per epoch stands in for the runner's fit so the
    averaging math is checked in isolation.
    """
    from types import SimpleNamespace
    from shapepipe.modules.ngmix_package.ngmix import average_original_psf

    def gal_epoch(psf_result, weight_value):
        # gal_obs carries a .psf whose meta['result'] is what the runner
        # would set; .weight is the galaxy inverse-variance map.
        return SimpleNamespace(
            psf=SimpleNamespace(meta={"result": psf_result}),
            weight=np.full((2, 2), weight_value),
        )

    # runner stub: a no-op .go, since meta['result'] is pre-populated. It
    # must be *called* (the fit happens), so record that it was.
    calls = []
    runner = SimpleNamespace(go=lambda obs: calls.append(obs))

    good_a = {
        "flags": 0, "T": 0.30, "T_err": 0.02,
        "g": np.array([0.04, -0.03]), "g_err": np.array([1e-3, 1e-3]),
    }
    good_b = {
        "flags": 0, "T": 0.20, "T_err": 0.01,
        "g": np.array([0.06, -0.01]), "g_err": np.array([2e-3, 2e-3]),
    }
    failed = {"flags": 3, "nfev": 5, "pars": np.zeros(6)}  # no T/g keys

    gal_obs_list = [
        gal_epoch(good_a, 1.0),       # galaxy weight-sum = 4
        gal_epoch(good_b, 2.0),       # galaxy weight-sum = 8
        gal_epoch(failed, 4.0),       # dropped on flags != 0
    ]
    psfo_res = average_original_psf(gal_obs_list, runner)

    # the runner fit every epoch's PSF
    assert len(calls) == 3
    # weighted over the two surviving epochs (weights 4 and 8)
    w = np.array([4.0, 8.0])
    npt.assert_allclose(
        psfo_res["g_psf"],
        (good_a["g"] * w[0] + good_b["g"] * w[1]) / w.sum(),
    )
    npt.assert_allclose(
        psfo_res["T_psf"],
        (good_a["T"] * w[0] + good_b["T"] * w[1]) / w.sum(),
    )
    assert psfo_res["n_epoch"] == 2
    # original PSF is elliptical — not the round reconvolution kernel
    assert abs(psfo_res["g_psf"][0]) > 1e-3


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
