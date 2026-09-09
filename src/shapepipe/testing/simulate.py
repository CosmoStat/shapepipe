"""Simulation utilities for ShapePipe validation tests.

Stable across shapepipe branches and ngmix versions — contains only galsim
and numpy, no shapepipe processing or ngmix fitting code.

Author: Axel Guinot
"""

import numpy as np
import galsim


def make_data(
    rng,
    shear,
    noise=1e-5,
    n_epochs=1,
    share_shift=False,
    gal_hlr=0.3,
    gal_flux=1000.0,
    psf_fwhm=0.55,
    pixel_scale=0.1857,
    img_size=201,
    psf_shear=(0.0, 0.0),
    wcs=None,
    return_centers=False,
    psf_model_fwhm_ratio=1.0,
    psf_model_shear=None,
):
    """Simulate an exponential galaxy with Moffat PSF.

    Parameters
    ----------
    rng : numpy.random.RandomState
        Random number generator.
    shear : tuple of float
        True shear (g1, g2).
    noise : float or numpy.ndarray, optional
        Noise sigma — a scalar, or a per-pixel map of shape
        ``(img_size, img_size)`` for spatially-varying noise. Default 1e-5.
    n_epochs : int, optional
        Number of epochs. Default 1.
    share_shift : bool, optional
        If True, all epochs share the same sub-pixel shift.
        If False, each epoch draws an independent shift.
    gal_hlr : float, optional
        Galaxy half-light radius in arcsec. Default 0.3.
    gal_flux : float, optional
        Galaxy flux. Default 1000.
    psf_fwhm : float, optional
        PSF FWHM in arcsec. Default 0.55.
    pixel_scale : float, optional
        Pixel scale in arcsec/pixel. Default 0.1857.
    img_size : int, optional
        Stamp size in pixels (square). Default 201.
    psf_shear : tuple of float, optional
        Ellipticity (g1, g2) applied to the Moffat PSF. Default (0, 0), a
        round PSF. A non-round PSF is what the PSF-leakage / additive-bias
        guardrail needs: the PSF carries a shape that an unbiased deconvolution
        must remove from a round galaxy.
    wcs : galsim.BaseWCS, optional
        Uniform drawing WCS. Default None means ``galsim.PixelScale
        (pixel_scale)`` with the legacy world-space sub-pixel shift draw —
        bit-identical rng consumption and values to earlier versions of this
        function. A non-trivial jacobian (rotation/shear — the ngmix#72
        sensitivity axis) switches the shift draw to PIXEL space (uniform
        ±0.5 pixel per axis), guaranteeing the object lands within half a
        pixel of the stamp centre — the contract that round-to-nearest-pixel
        centroid logic (e.g. ngmix ``centroid_source="wcs"``) depends on.
        Note ``psfs_sigmas`` are HSM pixel-unit sigmas of the *WCS-drawn*
        PSF stamp; under a sheared WCS they are not circularly symmetric
        measures.
    return_centers : bool, optional
        If True, also return per-epoch true object centres as
        ``galsim.PositionD`` in galsim image coordinates (1-based; stamp
        centre is ``(img_size + 1) / 2``). Default False keeps the legacy
        6-tuple return.
    psf_model_fwhm_ratio : float, optional
        Multiplicative size error injected into the PSF *model* stamp handed to
        ngmix, relative to the *true* PSF that convolves the galaxy. Default 1.0
        (no error: the same PSF object renders the galaxy and the model stamp,
        byte-for-byte identical to the legacy behaviour). A value != 1.0 makes
        the model PSF the wrong size — over-sized (>1) over-deconvolves, under-
        sized (<1) leaves residual smoothing — so metacal deconvolves by a PSF
        that does not match the one in the data. This is the controlled-PSF-
        model-error knob: it tests whether metacal's response correction (which
        is computed with the same wrong model) absorbs a deconvolution error,
        the one shape-measurement systematic the zero-error sim could not probe.
    psf_model_shear : tuple of float or None, optional
        Extra ellipticity (Δg1, Δg2) added to the PSF *model* stamp on top of
        ``psf_shear``, relative to the true PSF. Default None (no error). A
        non-zero value mismodels the PSF shape, injecting a controlled additive
        (PSF-leakage) error into the recovered shear. When both this is None and
        ``psf_model_fwhm_ratio`` is 1.0 the model PSF *is* the true PSF object,
        so the output is byte-for-byte the legacy result.

    Returns
    -------
    gals : list of numpy.ndarray
    psfs : list of numpy.ndarray
    psfs_sigmas : list of float
    weights : list of numpy.ndarray
    flags : list of numpy.ndarray
    jacob_lists : list of galsim.BaseWCS
    centers : list of galsim.PositionD, only if ``return_centers``
    """
    psf_noise = 1.0e-6
    scale = pixel_scale

    if wcs is None:
        wcs = galsim.PixelScale(scale)

        def draw_shift():
            # Legacy world-space draw — these exact lines (incl. dy-then-dx
            # order) keep rng consumption and values bit-identical to the
            # pre-wcs version of this function.
            dy, dx = rng.uniform(low=-scale / 2, high=scale / 2, size=2)
            return dx, dy, dx / scale, dy / scale

    else:
        jac = wcs.jacobian()
        transform = np.array([[jac.dudx, jac.dudy], [jac.dvdx, jac.dvdy]])

        def draw_shift():
            # Pixel-space draw (same dy-then-dx order): under a sheared
            # jacobian a world-space ±scale/2 box maps to a pixel-space
            # parallelogram exceeding ±0.5 pixel, which would silently break
            # rounded-centroid logic with off-by-one-pixel centres.
            dpy, dpx = rng.uniform(low=-0.5, high=0.5, size=2)
            dx, dy = transform @ (dpx, dpy)
            return float(dx), float(dy), dpx, dpy

    if share_shift:
        dx, dy, dpx, dpy = draw_shift()

    gals, psfs, psfs_sigmas, weights, flags, jacob_lists, centers = (
        [], [], [], [], [], [], []
    )

    for epoch in range(n_epochs):
        if not share_shift:
            dx, dy, dpx, dpy = draw_shift()

        psf = galsim.Moffat(beta=2.5, fwhm=psf_fwhm).shear(
            g1=psf_shear[0], g2=psf_shear[1]
        )
        obj = galsim.Convolve(
            psf,
            galsim.Exponential(half_light_radius=gal_hlr, flux=gal_flux).shear(
                g1=shear[0], g2=shear[1]
            ),
        ).shift(dx, dy)

        # The PSF *model* stamp handed to ngmix. By default it is the very same
        # object that convolved the galaxy (zero model error, byte-for-byte). A
        # size and/or shape error makes it differ from the true PSF, so metacal
        # deconvolves by the wrong PSF — the systematic the zero-error sim could
        # not exercise.
        if psf_model_fwhm_ratio == 1.0 and psf_model_shear is None:
            psf_model = psf
        else:
            psf_model = galsim.Moffat(
                beta=2.5, fwhm=psf_fwhm * psf_model_fwhm_ratio
            ).shear(
                g1=psf_shear[0] + (psf_model_shear[0] if psf_model_shear else 0.0),
                g2=psf_shear[1] + (psf_model_shear[1] if psf_model_shear else 0.0),
            )

        psf_im_ = psf_model.drawImage(nx=img_size, ny=img_size, wcs=wcs)
        psfs_sigmas.append(galsim.hsm.FindAdaptiveMom(psf_im_).moments_sigma)
        psf_im = psf_im_.array.astype(np.float64)
        im = obj.drawImage(nx=img_size, ny=img_size, wcs=wcs).array.astype(np.float64)

        psf_im += rng.normal(scale=psf_noise, size=psf_im.shape)
        im += rng.normal(scale=noise, size=im.shape)

        gals.append(im)
        psfs.append(psf_im)
        weights.append(im * 0 + 1.0 / noise ** 2)
        flags.append(im * 0)
        jacob_lists.append(wcs.jacobian())
        centers.append(galsim.PositionD(
            (img_size + 1) / 2 + dpx, (img_size + 1) / 2 + dpy
        ))

    if return_centers:
        return gals, psfs, psfs_sigmas, weights, flags, jacob_lists, centers
    return gals, psfs, psfs_sigmas, weights, flags, jacob_lists


def make_data_pujol(
    rng,
    shear_list,
    noise=1e-5,
    n_epochs=1,
    share_shift=False,
    gal_hlr=0.3,
    gal_flux=1000.0,
    psf_fwhm=0.55,
    pixel_scale=0.1857,
    img_size=201,
):
    """Simulate galaxies at multiple shear values sharing the same noise per epoch.

    Implements the Pujol (2018) denoising approach: all shear variants receive
    the identical noise and sub-pixel shift each epoch, so noise contributions
    to the shear response R and to the m/c estimator cancel exactly.

    Parameters
    ----------
    rng : numpy.random.RandomState
    shear_list : list of tuple of float
        List of (g1, g2) values to simulate simultaneously.
    noise : float, optional
        Per-pixel noise sigma.
    n_epochs : int, optional
    share_shift : bool, optional
        If True all epochs share the same sub-pixel shift; the shift is
        always shared across shear variants regardless of this flag.
    gal_hlr, gal_flux, psf_fwhm, pixel_scale, img_size : float / int, optional
        Same meaning as in make_data.

    Returns
    -------
    list of tuple
        One ``(gals, psfs, psfs_sigmas, weights, flags, jacob_lists)`` per
        entry in ``shear_list``, all sharing identical per-epoch noise.
    """
    psf_noise_sigma = 1.0e-6
    scale = pixel_scale
    wcs = galsim.PixelScale(scale)
    n_shears = len(shear_list)

    # output accumulators
    all_gals        = [[] for _ in range(n_shears)]
    all_psfs        = [[] for _ in range(n_shears)]
    all_psfs_sigmas = [[] for _ in range(n_shears)]
    all_weights     = [[] for _ in range(n_shears)]
    all_flags       = [[] for _ in range(n_shears)]
    all_jacobs      = [[] for _ in range(n_shears)]

    if share_shift:
        dy, dx = rng.uniform(low=-scale / 2, high=scale / 2, size=2)

    for epoch in range(n_epochs):
        if not share_shift:
            dy, dx = rng.uniform(low=-scale / 2, high=scale / 2, size=2)

        # One noise draw shared by all shear variants this epoch
        noise_img     = rng.normal(scale=noise,            size=(img_size, img_size))
        psf_noise_img = rng.normal(scale=psf_noise_sigma,  size=(img_size, img_size))

        psf    = galsim.Moffat(beta=2.5, fwhm=psf_fwhm)
        psf_im_ = psf.drawImage(nx=img_size, ny=img_size, wcs=wcs)
        psf_sigma = galsim.hsm.FindAdaptiveMom(psf_im_).moments_sigma
        psf_im = psf_im_.array.astype(np.float64) + psf_noise_img

        jacob  = wcs.jacobian()
        weight = np.full((img_size, img_size), 1.0 / noise ** 2)
        flag   = np.zeros((img_size, img_size))

        for j, shear in enumerate(shear_list):
            obj = galsim.Convolve(
                psf,
                galsim.Exponential(half_light_radius=gal_hlr, flux=gal_flux).shear(
                    g1=shear[0], g2=shear[1]
                ),
            ).shift(dx, dy)
            im = obj.drawImage(nx=img_size, ny=img_size, wcs=wcs).array.astype(np.float64)
            im += noise_img   # same noise for every shear variant

            all_gals[j].append(im)
            all_psfs[j].append(psf_im)
            all_psfs_sigmas[j].append(psf_sigma)
            all_weights[j].append(weight.copy())
            all_flags[j].append(flag.copy())
            all_jacobs[j].append(jacob)

    return [
        (all_gals[j], all_psfs[j], all_psfs_sigmas[j],
         all_weights[j], all_flags[j], all_jacobs[j])
        for j in range(n_shears)
    ]
