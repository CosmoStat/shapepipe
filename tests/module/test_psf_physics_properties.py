"""PROPERTY-BASED TESTS FOR NGMIX METACAL PSF PHYSICS.

Hypothesis tests over *random elliptical PSF stamps*. The metacal
reconvolution kernel is round and dilated by construction, so for ANY input
PSF ellipticity / size it must come back rounder and larger than the original
image-PSF fit:

    T_reconv >= T_orig    and    |g_reconv| <= |g_orig|.

This is a genuine metacal physics invariant (shapepipe#749) and, because the
two PSF families (``reconv`` -> ``average_multiepoch_psf``, ``orig`` ->
``average_original_psf``) share key names, it also catches any orig/reconv
routing transposition in :func:`do_ngmix_metacal` — a swap that is otherwise
silent. The single-example sibling lives in
``test_ngmix.test_reconv_psf_is_rounder_and_larger_than_orig_on_elliptical_psf``;
this widens it to a swept family of PSFs.
"""

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st


def _do_ngmix_metacal_on_psf(psf_shear, psf_fwhm=0.55, img_size=51, seed=7,
                             n_epochs=2):
    """Run do_ngmix_metacal on a simulated stamp with a given PSF shape.

    Mirrors ``test_ngmix._do_ngmix_metacal_on_psf`` (the harness this file is
    asked to reuse), exposing ``psf_fwhm`` so the PSF *size* can be swept too.
    The galaxy stamp is drawn from a fixed seed (123); only the metacal RNG
    varies with ``seed``. Returns ``(resdict, reconv, orig)`` — the
    :class:`MetacalResult` NamedTuple, still unpacking positionally.
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
        img_size=img_size,
        psf_fwhm=psf_fwhm,
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
    return do_ngmix_metacal(stamp, prior, 1.0, rng)


# Physically-sensible PSF strategies: ellipticity well inside |g| < 1 (real
# PSFs are mildly elliptical; a magnitude up to ~0.18 keeps the Moffat shear
# valid and the fit well-behaved), and a Moffat FWHM in arcsec spanning the
# realistic ground-based range. NB: the original-PSF fit uses the default galaxy
# prior (psf_fit_prior="galaxy"), which shrinks the *recovered* |g_orig| to
# O(1e-4) even for input |g|~0.1 (the shapepipe#749 prior-domination) -- so the
# dilation (T) inequality, not the rounding (|g|) one, is the primary guard.
_g_complex = st.complex_numbers(
    min_magnitude=0.02, max_magnitude=0.18, allow_nan=False, allow_infinity=False
)
_psf_fwhm = st.floats(min_value=0.5, max_value=0.7, allow_nan=False,
                      allow_infinity=False)


@settings(deadline=None, max_examples=25)
@given(g=_g_complex, psf_fwhm=_psf_fwhm)
def test_reconv_psf_rounder_and_larger_than_orig_over_random_psfs(g, psf_fwhm):
    """For ANY elliptical PSF, the reconv kernel is rounder and larger.

    Sweeps the input PSF ellipticity ``(g1, g2)`` (magnitude 0.02-0.18, any
    orientation) and the Moffat FWHM, runs the full :func:`do_ngmix_metacal`,
    and asserts the metacal reconvolution kernel is round + dilated relative to
    the independently-fit original image PSF:

        ``T_reconv >= T_orig``  and  ``|g_reconv| <= |g_orig|``.

    The ``T`` inequality is the robust transposition catcher: the kernel is
    dilated, so ``T_reconv`` exceeds ``T_orig`` by a clear margin (~0.01-0.03)
    and a strict ``>`` is asserted below -- a transposition makes ``T_reconv``
    the smaller value and fails it hard. The ``|g|`` inequality corroborates but
    with a thin margin: under the default galaxy prior the recovered ``|g_orig|``
    is only O(1e-4) (shapepipe#749 prior-domination), so it is a secondary check.
    """
    _, reconv, orig = _do_ngmix_metacal_on_psf(
        (g.real, g.imag), psf_fwhm=psf_fwhm
    )

    g_reconv = np.hypot(*reconv["g_psf"])
    g_orig = np.hypot(*orig["g_psf"])

    # Round: the reconvolution kernel carries less ellipticity than the
    # original PSF it was built to dominate (tiny tolerance for fit noise).
    assert g_reconv <= g_orig + 1e-6, (
        f"reconv |g|={g_reconv:.5f} not <= orig |g|={g_orig:.5f} "
        f"(input PSF g=({g.real:.3f},{g.imag:.3f}), fwhm={psf_fwhm:.3f})"
    )

    # Larger: the reconvolution kernel is strictly dilated relative to the
    # original PSF (observed dT ~0.01-0.03, far above fit noise); a transposition
    # makes T_reconv the smaller value and fails this hard.
    assert reconv["T_psf"] > orig["T_psf"], (
        f"reconv T={reconv['T_psf']:.5f} not dilated above orig T={orig['T_psf']:.5f} "
        f"(input PSF g=({g.real:.3f},{g.imag:.3f}), fwhm={psf_fwhm:.3f})"
    )
