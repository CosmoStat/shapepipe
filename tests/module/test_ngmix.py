"""UNIT TESTS FOR MODULE PACKAGE: NGMIX."""

from astropy.io import fits
from hypothesis import given
from hypothesis import strategies as st
import numpy as np
import numpy.testing as npt
import pytest
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


def _fake_metacal_result(T, T_err, T_psf, T_psf_err):
    """Build one minimal metacal result as produced by the process loop.

    Both PSF families are filled with self-named scalar keys, as the process
    loop back-fills them: ``*_psf_orig`` (original image PSF) and
    ``*_psf_reconv`` (metacal reconvolution kernel). The reconvolution-kernel
    size is the builder's ``T_psf`` argument; the original-PSF size is the
    fixed ``ORIG_PSF_T``, deliberately different so the un-aliasing is tested.
    """
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
        # metacal reconvolution-kernel family
        "g1_psf_reconv": RECONV_PSF_G[0],
        "g2_psf_reconv": RECONV_PSF_G[1],
        "g1_err_psf_reconv": RECONV_PSF_G_ERR[0],
        "g2_err_psf_reconv": RECONV_PSF_G_ERR[1],
        "T_psf_reconv": T_psf,
        "T_err_psf_reconv": T_psf_err,
        "mcal_flags": 0,
    }
    res.update(
        {name: dict(per_type) for name in ("1m", "1p", "2m", "2p", "noshear")}
    )
    return res


def test_compile_results_size_columns_are_unaliased():
    """Each PSF family carries its OWN fitted area ``T``, un-aliased (#749).

    The galaxy area ``T`` is the fitted ``pars[4]``; the reconvolution kernel
    (``*_psf_reconv``) and the original image PSF (``*_psf_orig``) are
    independent fits with their own, different sizes. Half-light radius is no
    longer carried in the catalogue — it is derivable from ``T`` downstream
    via :func:`cs_util.size.T_to_r50` — so only the ``T`` columns are checked.
    """
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    T, T_err = 0.18, 0.02
    T_psf, T_psf_err = 0.09, 0.001

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    out = inst.compile_results([_fake_metacal_result(T, T_err, T_psf, T_psf_err)])

    noshear = out["noshear"]

    # galaxy area is the fitted T
    npt.assert_allclose(noshear["T"], [T])
    npt.assert_allclose(noshear["T_err"], [T_err])

    # reconvolution kernel
    npt.assert_allclose(noshear["T_psf_reconv"], [T_psf])
    npt.assert_allclose(noshear["T_err_psf_reconv"], [T_psf_err])

    # original image PSF carries its OWN size, un-aliased from the kernel
    npt.assert_allclose(noshear["T_psf_orig"], [ORIG_PSF_T])
    npt.assert_allclose(noshear["T_err_psf_orig"], [ORIG_PSF_T_ERR])
    assert noshear["T_psf_orig"] != noshear["T_psf_reconv"]


def test_compile_results_psf_families_are_unaliased():
    """The two PSF-ellipticity families document *different* PSFs (#749).

    ``*_psf_reconv`` carries the metacal RECONVOLUTION kernel; ``*_psf_orig``
    carries the ORIGINAL image PSF, fit independently. Before the fix both
    came from one ``average_multiepoch_psf`` result and so were byte-
    identical; here they must differ, each tracing its own source. The
    companion size un-aliasing is checked in
    ``test_compile_results_size_columns_are_unaliased``.
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
        "s2n", "mag", "mag_err",
    ):
        assert np.isnan(failed[col]).all(), col
    assert failed["flags"] == [0x8]
    assert failed["mcal_flags"] == [0x8] and failed["mcal_flags"][0] != 0

    # The object-level PSF columns survive on a failed fit type: they are copied
    # outside the flags==0 guard, so the failed type keeps a full-length PSF row
    # (gating that copy on the galaxy flags would shorten it -> save crash).
    npt.assert_allclose(failed["T_psf_orig"], [ORIG_PSF_T])
    npt.assert_allclose(failed["T_psf_reconv"], [0.09])
    npt.assert_allclose(failed["g1_psf_orig"], [ORIG_PSF_G[0]])
    npt.assert_allclose(failed["g1_psf_reconv"], [RECONV_PSF_G[0]])

    # successful types are untouched, but share the object's mcal_flags
    npt.assert_allclose(out["noshear"]["flux"], [100.0])
    assert out["noshear"]["flags"] == [0]
    assert out["noshear"]["mcal_flags"] == [0x8]


def test_save_results_batches_survive_a_later_failed_fit(tmp_path):
    """A batch-1-all-success / batch-2-has-a-failure sequence must not crash.

    With ``SAVE_BATCH`` enabled, ``save_results`` writes batch 1 fresh via
    ``save_as_fits`` (which locks column dtypes from that data) and appends
    every later batch into a structured array typed to those locked dtypes
    (shapepipe#795). ``nfev_fit`` was the only int-like column with a NaN
    fallback: if batch 1 has zero failed fits, ``nfev_fit`` locks to int64,
    and the first later batch containing a failed fit fed a NaN into that
    int64 column, raising ``ValueError: cannot convert float NaN to integer``
    at append time. The fix uses the -1 sentinel (ngmix's own
    failed/absent-nfev convention) instead of NaN, so the column stays
    int64 across every batch.
    """
    from shapepipe.modules.ngmix_package.ngmix import Ngmix

    inst = object.__new__(Ngmix)
    inst._zero_point = 30.0
    inst._output_dir = str(tmp_path)
    inst._file_number_string = "-000-000"

    # Batch 1: every fit type succeeds -> nfev_fit locks to int64 on write.
    batch1 = [_fake_metacal_result(0.18, 0.02, 0.09, 0.001)]
    inst.save_results(inst.compile_results(batch1))

    # Batch 2: one fit type fails (no "nfev" key -> hits the fallback path).
    batch2_res = _fake_metacal_result(0.18, 0.02, 0.09, 0.001)
    batch2_res["1p"] = {"flags": 0x8}  # failed fit: no nfev at all
    batch2_res["mcal_types_fail"] = 1
    batch2_res["mcal_flags"] = 0x8

    # Must not raise (pre-fix: ValueError: cannot convert float NaN to integer).
    inst.save_results(inst.compile_results([batch2_res]))

    with fits.open(inst.get_output_path(str(tmp_path))) as hdul:
        nfev = hdul["1P"].data["nfev_fit"]
        assert nfev.dtype.kind == "i"
        npt.assert_array_equal(nfev, [1, -1])


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
        # would set; .weight is the galaxy inverse-variance map. The PSF
        # observation must expose .copy() (average_original_psf fits a copy
        # so gal_obs.psf is never mutated); the copy carries the same
        # pre-populated meta so the averaging math is unchanged.
        psf_obs = SimpleNamespace(meta={"result": psf_result})
        psf_obs.copy = lambda _self=psf_obs: SimpleNamespace(
            meta=dict(_self.meta)
        )
        return SimpleNamespace(
            psf=psf_obs,
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
    psf_orig_res = average_original_psf(gal_obs_list, runner)

    # the runner fit every epoch's PSF
    assert len(calls) == 3
    # weighted over the two surviving epochs (weights 4 and 8)
    w = np.array([4.0, 8.0])
    npt.assert_allclose(
        psf_orig_res["g_psf"],
        (good_a["g"] * w[0] + good_b["g"] * w[1]) / w.sum(),
    )
    npt.assert_allclose(
        psf_orig_res["T_psf"],
        (good_a["T"] * w[0] + good_b["T"] * w[1]) / w.sum(),
    )
    assert psf_orig_res["n_epoch"] == 2
    # original PSF is elliptical — not the round reconvolution kernel
    assert abs(psf_orig_res["g_psf"][0]) > 1e-3


def _do_ngmix_metacal_on_psf(psf_shear, seed=7, n_epochs=2):
    """Run do_ngmix_metacal on a simulated stamp with a given PSF shape.

    Returns ``(resdict, psf_reconv, psf_orig, gal_obs_list)``, where the last
    is the list metacal consumed — so a caller can assert that the original-PSF
    pre-fit left it pristine.
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
        shear=(0.0, 0.0),
        noise=1e-5,
        n_epochs=n_epochs,
        img_size=51,
        psf_shear=psf_shear,
    )
    stamp = Postage_stamp(bkg_sub=False, megacam_flip=False)
    stamp.gals, stamp.psfs, stamp.weights, stamp.flags, stamp.jacobs = (
        gals,
        psfs,
        weights,
        flags,
        jacobs,
    )
    resdict, psf_reconv, psf_orig = do_ngmix_metacal(stamp, prior, 1.0, rng)
    return resdict, psf_reconv, psf_orig


def test_original_psf_prefit_leaves_gal_obs_psf_pristine():
    """The original-PSF pre-fit must not mutate the PSF metacal consumes (#749).

    ``average_original_psf`` fits each epoch's psfex/mccd PSF with the shared
    ``psf_runner``. ``PSFRunner.go`` sets ``.gmix`` (and ``.meta['result']``)
    on the observation it fits; if that were ``gal_obs.psf`` itself, the gmix
    would survive metacal's deep copy and be reused as the
    ``MetacalFitGaussPSF`` fallback when admom+ML both fail — silently
    rescuing objects the base branch dropped (``BootPSFFailure``). The fix
    fits a *copy*, so ``gal_obs.psf`` carries NO gmix afterward, matching the
    base-branch state that makes this add-column refactor bit-identical on the
    galaxy results.
    """
    from ngmix.observation import ObsList
    from shapepipe.modules.ngmix_package.ngmix import (
        average_original_psf,
        get_prior,
        make_ngmix_observation,
        make_runners,
    )
    from shapepipe.testing.simulate import make_data

    rng = np.random.RandomState(7)
    prior = get_prior(0.1857, rng)
    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(123),
        shear=(0.0, 0.0),
        noise=1e-5,
        n_epochs=2,
        img_size=51,
        psf_shear=(0.05, -0.04),
    )
    gal_obs_list = ObsList()
    for n_e in range(2):
        gal_obs_list.append(
            make_ngmix_observation(
                gals[n_e], weights[n_e], flags[n_e], psfs[n_e], jacobs[n_e],
                rng,
            )
        )

    # Pre-state: a freshly built PSF observation carries no gmix.
    assert not any(obs.psf.has_gmix() for obs in gal_obs_list)

    runner, psf_runner = make_runners(prior, 1.0, rng)
    average_original_psf(gal_obs_list, psf_runner)

    # The pre-fit fit a copy, so the originals are untouched: no gmix and no
    # leaked fit result in meta.
    assert not any(obs.psf.has_gmix() for obs in gal_obs_list), (
        "original-PSF pre-fit leaked a gmix into gal_obs.psf"
    )
    assert not any("result" in obs.psf.meta for obs in gal_obs_list), (
        "original-PSF pre-fit leaked a fit result into gal_obs.psf.meta"
    )


def test_reconv_psf_is_rounder_and_larger_than_orig_on_elliptical_psf():
    """Reconvolution kernel is round + dilated relative to the original PSF.

    On an ELLIPTICAL PSF stamp the metacal reconvolution kernel is round and
    enlarged by construction, so ``T_psf_reconv > T_psf_orig`` and
    ``|g_psf_reconv| < |g_psf_orig|``. This end-to-end guard catches a
    ``psf_res`` / ``psf_orig_res`` transposition in :func:`do_ngmix_metacal`
    (the two families share keys, so a swap is otherwise silent).
    """
    _, psf_reconv, psf_orig = _do_ngmix_metacal_on_psf((0.05, -0.04))

    assert psf_reconv["T_psf"] > psf_orig["T_psf"]
    assert np.hypot(*psf_reconv["g_psf"]) < np.hypot(*psf_orig["g_psf"])


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


def _one_epoch_stamp():
    """One-epoch (gal, weight, flag, psf, galsim-jacobian) for centroid tests."""
    from shapepipe.testing.simulate import make_data

    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(1),
        shear=(0.0, 0.0),
        noise=1e-5,
        n_epochs=1,
        img_size=51,
        psf_shear=(0.0, 0.0),
    )
    return gals[0], weights[0], flags[0], psfs[0], jacobs[0]


def test_wcs_centroid_uses_stored_offset():
    """``centroid_source='wcs'`` places the galaxy Jacobian origin at the stamp
    center plus the propagated coadd-centroid offset — it does not re-derive
    the centroid by re-projecting through the WCS (the #767 fix)."""
    from shapepipe.modules.ngmix_package.ngmix import make_ngmix_observation

    gal, weight, flag, psf, jacob = _one_epoch_stamp()
    offset = np.array([0.3, -0.2])  # [row, col]

    obs = make_ngmix_observation(
        gal, weight, flag, psf, jacob, np.random.RandomState(7),
        centroid_source="wcs", offset=offset,
    )

    row0, col0 = obs.jacobian.get_cen()
    center = (gal.shape[0] - 1) / 2
    npt.assert_allclose(
        [row0, col0], [center + offset[0], center + offset[1]]
    )


def test_wcs_centroid_requires_offset():
    """``centroid_source='wcs'`` with no propagated offset fails fast."""
    from shapepipe.modules.ngmix_package.ngmix import make_ngmix_observation

    gal, weight, flag, psf, jacob = _one_epoch_stamp()

    with pytest.raises(ValueError, match="requires the coadd-centroid offset"):
        make_ngmix_observation(
            gal, weight, flag, psf, jacob, np.random.RandomState(7),
            centroid_source="wcs", offset=None,
        )
