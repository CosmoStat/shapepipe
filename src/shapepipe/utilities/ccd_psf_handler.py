"""CCD_PSF_HANDLER

Obtain list of CCDs (single-exposure single-HDU files) for which valid PSF
information is available. This can serve to create a footprint coverage mask.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
import numpy as np

from cs_util import args as cs_args
from cs_util import logging

from shapepipe.utilities import summary


class CcdPsfHandler(object):
    """CCD PSF Handler Class.

    Handles extraction of CCDs with valid PSF information from shapepipe
    patches.
    """

    def __init__(self):
        """Initialize the handler."""
        self.params_default()

    def params_default(self):
        """Set default parameters and command line options."""

        self._params = {
            "version_cat": "v1.6",
            "n_CCD": 40,
            "output": None,
        }

        self._short_options = {
            "version_cat": "-V",
            "n_CCD": "-n",
            "output": "-o",
        }

        self._types = {
            "n_CCD": "int",
        }

        self._help_strings = {
            "version_cat": "catalogue major version, allowed are v1.3, v1.4, v1.5, v1.6; default is {}",
            "n_CCD": "number of CCDs per exposure; default is {}",
            "output": "output file path; default is ccds_with_psf_<version>.txt",
        }

    def set_params_from_command_line(self, args):
        """Set Params From Command line.

        Only use when calling using python from command line.
        Does not work from ipython or jupyter.

        Parameters
        ----------
        args : list
            command line arguments

        """
        # Read command line options
        options = cs_args.parse_options(
            self._params,
            self._short_options,
            self._types,
            self._help_strings,
            args=args,
        )
        self._params = options

        # Save calling command
        logging.log_command(args)

    def update_params(self):
        """Update parameters.

        Set derived parameters based on input parameters.
        """
        # Determine number of patches based on version
        version = self._params["version_cat"]
        if version in ("v1.3", "v1.4"):
            n_patch = 7
        elif version == "v1.5":
            n_patch = 8
        elif version == "v1.6":
            n_patch = 9
        else:
            raise ValueError(f"Invalid version {version}")

        self._params["n_patch"] = n_patch
        self._params["patches"] = [f"P{x}" for x in np.arange(n_patch) + 1]

        # Set output file if not specified
        if self._params["output"] is None:
            self._params["output"] = f"ccds_with_psf_{version}.txt"

    def check_params(self):
        """Check parameters for validity."""
        # Add any parameter validation logic here
        pass

    def get_lines(self, fname):
        """Get Lines.

        Return list of lines read from a text file.

        Parameters
        ----------
        fname : str
            input file name

        Returns
        -------
        list
            IDs

        """
        IDs = []
        with open(fname) as f:
            lines = f.readlines()
            for line in lines:
                IDs.append(line.rstrip())

        return IDs

    def get_exp_shdu_missing(self, patches):
        """Get Exp Shdu Missing.

        Returns set of missing CCDs (single-exposure single-HDU IDs) from a list of patches.

        Parameters
        ----------
        patches : list
            input patches

        Returns
        -------
        set
           missing CCD IDs

        """
        exp_shdu_missing_all = set()

        for patch in patches:

            path_exp_shdu_missing = f"{patch}/summary/missing_job_32_all.txt"
            exp_shdu_missing = self.get_lines(path_exp_shdu_missing)

            print(
                f"Patch {patch}: Found {len(exp_shdu_missing)} missing ccds",
                end="; ",
            )

            exp_shdu_missing_all.update(exp_shdu_missing)

            print(f"cumulative {len(exp_shdu_missing_all)} missing ccds")

        print()

        return exp_shdu_missing_all

    def get_exp(self, patches):
        """Get Exp.

        Return set of exposures from a list of patches.

        Parameters
        ----------
        patches : list
            input patches

        Returns
        -------
        set
           exposure IDs

        """
        exp_all = set()

        for patch in patches:

            path_exp = f"{patch}/exp_numbers.txt"
            exp = self.get_lines(path_exp)

            print(f"Patch {patch}: Found {len(exp)} exposures", end="; ")

            exp_all.update(exp)

            print(f"cumulative {len(exp_all)} exposures")

        print()

        return exp_all

    def get_ccds_with_psf(self, patches, n_CCD=40):
        """Get CCDs With PSF.

        Return set of CCDs from list of patches.

        Parameters
        ----------
        patches : list
            input patches
        n_CCD : int
            number of CCDs per exposure

        Returns
        -------
        set
           CCD IDs

        """
        # Get missing CCDs
        print("=== get missing CCDs ===")
        exp_shdu_missing_all = self.get_exp_shdu_missing(patches)

        # Get all exposures used in tiles
        print("=== get exposures ===")
        exp_all = self.get_exp(patches)

        # Turn exposures into exposure-single-HDU names (CCDs)
        exp_shdu_all = summary.get_all_shdus(exp_all, n_CCD)

        print(f"Found {len(exp_shdu_all)} CCDs")

        return exp_shdu_all

    def save(self, IDs, path):
        """Save.

        Save list of IDs to text file.

        Parameters
        ----------
        IDs : set
            input IDs
        path : str
            output file name

        """
        with open(path, "w") as f_out:
            for ID in IDs:
                print(ID, file=f_out)

    def run(self, args=None):
        """Run.

        Main execution method.

        Parameters
        ----------
        args : list, optional
            command line arguments

        Returns
        -------
        int
            exit code (0 for success)

        """
        if args is None:
            args = sys.argv

        # Set parameters from command line
        self.set_params_from_command_line(args)
        self.update_params()
        self.check_params()

        # Get parameters
        patches = self._params["patches"]
        version = self._params["version_cat"]
        n_CCD = self._params["n_CCD"]
        output = self._params["output"]

        print(
            f"=== get_ccds_with_psf for version {version}, patches {patches} ==="
        )

        print("=== method 1: exp_list - missing === ")
        exp_shdu_all = self.get_ccds_with_psf(patches, n_CCD)

        self.save(exp_shdu_all, output)

        print(f"Results saved to {output}")

        return 0
