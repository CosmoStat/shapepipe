#!/usr/bin/env python3

"""GET_CCDS_WITH_PSF

Obtain list of CCDs (single-exposure single-HDU files) for which valid PSF information
is available. This can serve to create a footprint coverage mask.

Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys

from shapepipe.utilities.ccd_psf_handler import CcdPsfHandler


def run_ccd_psf_handler(args=None):
    """Run CCD PSF Handler.

    Create instance and run the CCD PSF handler.

    Parameters
    ----------
    args : list, optional
        command line arguments

    Returns
    -------
    int
        exit code

    """
    # Create instance

    obj = CcdPsfHandler()

    return obj.run(args=args)


def main(argv=None):
    """Main.
    """
    # A scripts to call the ccd psf class is created by pyproject.toml
    return 0
