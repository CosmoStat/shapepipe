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
):
    """Simulate an exponential galaxy with Moffat PSF.

    Parameters
    ----------
    rng : numpy.random.RandomState
        Random number generator.
    shear : tuple of float
        True shear (g1, g2).
    noise : float, optional
        Per-pixel noise sigma. Default 1e-5.
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

    Returns
    -------
    gals : list of numpy.ndarray
    psfs : list of numpy.ndarray
    psfs_sigmas : list of float
    weights : list of numpy.ndarray
    flags : list of numpy.ndarray
    jacob_lists : list of galsim.BaseWCS
    """
    psf_noise = 1.0e-6
    scale = pixel_scale
    wcs = galsim.PixelScale(scale)

    if share_shift:
        dy, dx = rng.uniform(low=-scale / 2, high=scale / 2, size=2)

    gals, psfs, psfs_sigmas, weights, flags, jacob_lists = [], [], [], [], [], []

    for epoch in range(n_epochs):
        if not share_shift:
            dy, dx = rng.uniform(low=-scale / 2, high=scale / 2, size=2)

        psf = galsim.Moffat(beta=2.5, fwhm=psf_fwhm)
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

    return gals, psfs, psfs_sigmas, weights, flags, jacob_lists
