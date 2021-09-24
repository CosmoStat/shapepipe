# -*- coding: utf-8 -*-

"""GALAXY TOOLS

This module defines methods to deal with galaxy images.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import numpy as np


def sigma_to_fwhm(sigma, pixel_scale=1.0):
    """Sigma to FWHM.

    Transform standard deviation of a 1D Gaussian, sigma, to FWHM (Full Width Half Maximum).

    Parameters
    ----------
    sigma : numpy.ndarray
        input standard deviation(s)
    pixel_scale : float, optional, default=1
        pixel size in arcsec, set to 1 if no scaling
        required

    Returns
    -------
    fwhm : (array of) float
        output fwhm(s)

    Raises
    ------
    TypeError
        If ``sigma`` is not of type numpy array or float
    TypeError
        If ``sigma`` array values are not of type float
    TypeError
        If ``pixel_scale`` is not of type float

    Notes
    -----
    To compute the FWHMh for a 1D Gaussian N(x), solve the equation

    ..math::

        N(x) = (\sigma \sqrt{2\pi})^{-1} \exp[x^2/2\sigma^2] = \frac 1 2 N(x)


    for :math:`x`. The FWHM is :math:`x + (-x) = 2x`. The solution is

    ..math::

        \textrm{FWHM} = 2 \sqrt(2 \ln 2) \sigma \approx 2.355 \sigma

    """

    if not isinstance(sigma, (np.ndarray, float)):
        raise TypeError(
            f'Sigma must be of type numpy array or float, not {type(sigma)}.'
        )
    elif isinstance(sigma, np.ndarray) and sigma.dtype != np.float64:
        raise TypeError(
            f'Sigma array values must be of type float, not {sigma.dtype}.'
        )

    if not isinstance(pixel_scale, float):
        raise TypeError(
            f'The pixel scale must of type float, not {type(sigma)}.'
        )

    if pixel_scale <= 0:
        raise ValueError(
            f'Invalid pixel scale {pixel_scale}, needs to be positive'
        )
    if sigma <= 0:
        raise ValueError(
            f'Invalid standard deviatoin {sigma}, needs to be positive'
        )

    cst = 2.35482004503

    return sigma * cst * pixel_scale
