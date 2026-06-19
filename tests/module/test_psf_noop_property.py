"""PROPERTY-BASED TESTS: the original-PSF pre-fit is a science no-op.

``average_original_psf`` fits each epoch's psfex/mccd PSF stamp
(``gal_obs.psf``) with the shared ``psf_runner`` to populate the independent
original-PSF columns (``NGMIX_*_PSF_ORIG``; shapepipe#749). It fits a *copy* of
each PSF observation, so the observation metacal later deep-copies and consumes
(``boot.go(gal_obs_list)``) is left untouched.

This is the invariant that makes the add-column refactor bit-identical to a
no-original-fit branch on the galaxy/shear results: if the pre-fit cannot alter
``gal_obs.psf``, it cannot alter what metacal consumes, so the galaxy results
are unchanged. ``PSFRunner.go`` sets ``.gmix`` (and ``.meta['result']``) on the
observation it fits; were that ``gal_obs.psf`` itself, a stray gmix would
survive metacal's deep copy and be reused as the ``MetacalFitGaussPSF`` fallback
when admom+ML both fail, silently rescuing objects the base branch dropped
(``BootPSFFailure``).

The unit guard ``test_original_psf_prefit_leaves_gal_obs_psf_pristine`` in
``tests/module/test_ngmix.py`` pins this on one fixed input; here it is
generalised over varied PSF shapes, noise, epoch counts, stamp sizes, and RNG
seeds with hypothesis, so the no-op holds across the input space rather than at
a single point. The harness (``make_data`` -> ``make_ngmix_observation`` ->
``make_runners`` -> ``average_original_psf``) mirrors that unit test exactly.
"""

import copy

from hypothesis import given, settings
from hypothesis import strategies as st
import numpy as np
import numpy.testing as npt
from ngmix.observation import ObsList

from shapepipe.modules.ngmix_package.ngmix import (
    average_original_psf,
    get_prior,
    make_ngmix_observation,
    make_runners,
)
from shapepipe.testing.simulate import make_data

PIXEL_SCALE = 0.1857


def _jacobian_snapshot(jacobian):
    """The six floats fully defining an ngmix Jacobian (affine + origin)."""
    return (
        jacobian.dvdrow, jacobian.dvdcol,
        jacobian.dudrow, jacobian.dudcol,
        jacobian.row0, jacobian.col0,
    )


def _build_gal_obs_list(psf_shear, noise, n_epochs, img_size, data_seed, obs_seed):
    """Build the per-epoch galaxy ObsList metacal would consume.

    Mirrors the construction in
    ``test_original_psf_prefit_leaves_gal_obs_psf_pristine``: simulate epochs
    with ``make_data`` (here over hypothesis-varied PSF shape, noise, epoch
    count, stamp size, and seed), then wrap each through
    ``make_ngmix_observation`` so each ``gal_obs.psf`` is a real, freshly built
    PSF observation carrying no gmix.
    """
    rng = np.random.RandomState(obs_seed)
    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(data_seed),
        shear=(0.0, 0.0),
        noise=noise,
        n_epochs=n_epochs,
        img_size=img_size,
        psf_shear=psf_shear,
    )
    gal_obs_list = ObsList()
    for n_e in range(n_epochs):
        gal_obs_list.append(
            make_ngmix_observation(
                gals[n_e], weights[n_e], flags[n_e], psfs[n_e], jacobs[n_e],
                rng,
            )
        )
    return gal_obs_list


# Physically sensible strategies: |g| well below 1 (Moffat shear must stay
# realisable), positive noise, >=1 epoch, odd stamp large enough to host the
# 0.55" FWHM PSF, and bounded RNG seeds.
psf_g_component = st.floats(
    min_value=-0.3, max_value=0.3, allow_nan=False, allow_infinity=False
)
noise_strategy = st.floats(
    min_value=1e-6, max_value=1e-3, allow_nan=False, allow_infinity=False
)
n_epochs_strategy = st.integers(min_value=1, max_value=3)
img_size_strategy = st.sampled_from([41, 51, 65])
seed_strategy = st.integers(min_value=0, max_value=2**31 - 1)


@given(
    g1=psf_g_component,
    g2=psf_g_component,
    noise=noise_strategy,
    n_epochs=n_epochs_strategy,
    img_size=img_size_strategy,
    data_seed=seed_strategy,
    obs_seed=seed_strategy,
)
@settings(deadline=None, max_examples=25)
def test_average_original_psf_leaves_gal_obs_psf_pristine_property(
    g1, g2, noise, n_epochs, img_size, data_seed, obs_seed
):
    """For arbitrary inputs, the original-PSF pre-fit leaves gal_obs.psf pristine.

    The no-op invariant that guarantees bit-identical galaxy/shear results:
    after ``average_original_psf``, every epoch's ``gal_obs.psf`` is unchanged
    in every channel metacal reads —

    * ``has_gmix()`` is still False (no fitted gmix leaked in), and no fit
      ``result`` leaked into ``.meta``;
    * the PSF ``image``, ``weight``, and ``jacobian`` are bit-for-bit what they
      were before the call.

    Each assertion is load-bearing: were the pre-fit to fit ``gal_obs.psf``
    itself instead of a copy, ``PSFRunner.go`` would set ``.gmix`` /
    ``.meta['result']`` (first two assertions catch that), and centroid/weight
    handling inside the fit could perturb the arrays (the array/jacobian
    snapshots catch that). The fit still *runs* on every epoch — this pins that
    it runs harmlessly on a copy, not that it is skipped.
    """
    psf_shear = (g1, g2)
    rng = np.random.RandomState(obs_seed)
    prior = get_prior(PIXEL_SCALE, rng)
    gal_obs_list = _build_gal_obs_list(
        psf_shear, noise, n_epochs, img_size, data_seed, obs_seed
    )

    # Pre-state: a freshly built PSF observation carries no gmix or result, and
    # we snapshot every array channel metacal will consume.
    assert not any(obs.psf.has_gmix() for obs in gal_obs_list)
    pre = [
        {
            "image": obs.psf.image.copy(),
            "weight": obs.psf.weight.copy(),
            "jacobian": _jacobian_snapshot(obs.psf.jacobian),
            "meta": copy.deepcopy(dict(obs.psf.meta)),
        }
        for obs in gal_obs_list
    ]

    runner, psf_runner = make_runners(prior, 1.0, rng)
    psf_orig_res = average_original_psf(gal_obs_list, psf_runner)

    # The pre-fit actually did work (so the no-op is non-vacuous): at least one
    # epoch's PSF was fit and averaged in.
    assert psf_orig_res["n_epoch"] >= 1

    for obs, snapshot in zip(gal_obs_list, pre):
        psf = obs.psf
        # No fitted gmix leaked into the observation metacal consumes.
        assert not psf.has_gmix(), (
            "original-PSF pre-fit leaked a gmix into gal_obs.psf"
        )
        # No fit result leaked into its metadata.
        assert "result" not in psf.meta, (
            "original-PSF pre-fit leaked a fit result into gal_obs.psf.meta"
        )
        assert dict(psf.meta) == snapshot["meta"], (
            "original-PSF pre-fit mutated gal_obs.psf.meta"
        )
        # The image, weight, and jacobian metacal reads are bit-for-bit intact.
        npt.assert_array_equal(
            psf.image, snapshot["image"],
            err_msg="original-PSF pre-fit mutated gal_obs.psf.image",
        )
        npt.assert_array_equal(
            psf.weight, snapshot["weight"],
            err_msg="original-PSF pre-fit mutated gal_obs.psf.weight",
        )
        assert _jacobian_snapshot(psf.jacobian) == snapshot["jacobian"], (
            "original-PSF pre-fit mutated gal_obs.psf.jacobian"
        )
