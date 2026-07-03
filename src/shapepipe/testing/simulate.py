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
