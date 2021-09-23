# -*- coding: utf-8 -*-

"""GALAXY TOOLS

This module defines methods to deal with galaxy images.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import numpy as np


def sigma_to_fwhm(sigma, pixel_scale=1.0):
    """Sigma to FWHM.

    Transform from size sigma to FWHM. The conversion factor corresponds to a
    1D Gaussian profile.

    Parameters
    ----------
    sigma : numpy.ndarray
        input size(s)
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

    cst = 2.35482004503

    return sigma * cst * pixel_scale
