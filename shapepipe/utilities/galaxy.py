# -*- coding: utf-8 -*-

"""GALAXY TOOLS

This module defines methods to deal with galaxy images.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 06/2021

:Package: ShapePipe

"""


def sigma_to_fwhm(sigma, pixel_scale=1):
    """sigma to fwhm

    Transform from size sigma to FWHM.
    The conversion factor corresponds to a 1D Gaussian profile.

    Parameters
    ----------
    sigma : (array of) float
        input size(s)
    pixel_scale : float, optional, default=1
        pixel size in arcsec, set to 1 if no scaling
        required

    Returns
    -------
    fwhm : (array of) float
        output fwhm(s)
    """

    # cst = 2 * sqrt(2 * ln(2))
    cst = 2.35482004503

    return sigma * cst * pixel_scale
