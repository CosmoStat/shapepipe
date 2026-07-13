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
        must remove from a round galaxy. The same sheared PSF is used both to
        convolve the object and as the PSF model stamp.
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

        psf_im_ = psf.drawImage(nx=img_size, ny=img_size, wcs=wcs)
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
