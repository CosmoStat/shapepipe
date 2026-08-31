"""Guardrail tests for the PSF-model-error extension of ``make_data``.

``make_data`` renders the galaxy with the *true* PSF and, by default, hands ngmix
that same PSF as the model stamp (zero model error). The ``psf_model_fwhm_ratio``
and ``psf_model_shear`` knobs let the model stamp differ from the true PSF, so a
controlled deconvolution error can be injected. These tests pin the two
properties the digital twin's PSF-model-error sweep relies on:

  1. with no error injected the output is byte-for-byte the legacy result, and
  2. an injected error perturbs ONLY the PSF model stamp, never the galaxy image
     (the galaxy is always convolved with the true PSF).
"""

import numpy as np

from shapepipe.testing.simulate import make_data

_KW = dict(shear=(0.02, 0.0), noise=1e-4, n_epochs=2, img_size=51,
           pixel_scale=0.1857)


def _run(seed=123, **extra):
    gals, psfs, *_ = make_data(rng=np.random.RandomState(seed), **_KW, **extra)
    return np.array(gals), np.array(psfs)


def test_no_error_is_byte_for_byte():
    """Explicit no-error kwargs reproduce the implicit default exactly."""
    g0, p0 = _run()
    g1, p1 = _run(psf_model_fwhm_ratio=1.0, psf_model_shear=None)
    np.testing.assert_array_equal(g0, g1)
    np.testing.assert_array_equal(p0, p1)


def test_size_error_touches_only_the_psf_stamp():
    """A model size error changes the PSF stamp but leaves the galaxy intact."""
    g0, p0 = _run()
    g1, p1 = _run(psf_model_fwhm_ratio=1.05)
    np.testing.assert_array_equal(g0, g1)          # galaxy uses the TRUE PSF
    assert not np.allclose(p0, p1)                  # model stamp is perturbed


def test_ellipticity_error_touches_only_the_psf_stamp():
    """A model ellipticity error changes the PSF stamp, not the galaxy."""
    g0, p0 = _run()
    g1, p1 = _run(psf_model_shear=(0.03, 0.0))
    np.testing.assert_array_equal(g0, g1)
    assert not np.allclose(p0, p1)
