# -*- coding: utf-8 -*-

"""GALAXY TOOLS

This module defines methods to deal with galaxy images.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 06/2021

:Package: ShapePipe

"""


def sigma_to_fwhm(sigma, pixel_size=1):
    """sigma to fwhm

    Transform from size sigma to FWHM.
    The conversion factor corresponds to a 1D Gaussian profile.

    Parameters
    ----------
    sigma : (array of) float
        input size(s)
    pixel_size : float, optional, default=1
        pixel size in arcsec, set to 1 if no scaling
        required

    Returns
    -------
    fwhm : (array of) float
        output fwhm(s)
    """

    return sigma * 2.355 * pixel_size
