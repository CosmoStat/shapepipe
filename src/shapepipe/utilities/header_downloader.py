"""HEADER_DOWNLOADER

Download FITS headers from VOSpace for exposures listed in a CCD file.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
import os
from os.path import exists

import numpy as np
import vos

from cs_util import args as cs_args
from cs_util import logging


class HeaderDownloader(object):
    """Header Downloader Class.

    Downloads FITS headers from VOSpace for exposures in a CCD list.
    """

    def __init__(self):
        """Initialize the downloader."""
        self.params_default()

    def params_default(self):
        """Set default parameters and command line options."""

        self._params = {
            "input_file": None,
            "output_dir": "header",
            "vospace_path": "vos:cfis/pitcairn", # Default UNIONS exposure path
            "overwrite": False,
            "dir_for_links": None,
        }

        self._short_options = {
            "input_file": "-i",
            "output_dir": "-o",
            "vospace_path": "-p",
            "overwrite": "-O",
            "dir_for_links": "-d",
        }

        self._types = {
            "overwrite": "bool",
        }

        self._help_strings = {
            "input_file": "input CCD list file (txt or csv); required",
            "output_dir": "output directory for headers; default is {}",
            "vospace_path": "VOSpace base path; default is {}",
            "overwrite": "overwrite existing header files; default is {}",
            "dir_for_links": "directory to check for existing headers to link instead of download; default is {}",
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
        # Ensure output directory ends without trailing slash
        if self._params["output_dir"].endswith("/"):
            self._params["output_dir"] = self._params["output_dir"][:-1]

    def check_params(self):
        """Check parameters for validity."""
        if self._params["input_file"] is None:
            raise ValueError("Input file is required (use -i or --input_file)")

        if not exists(self._params["input_file"]):
            raise FileNotFoundError(
                f"Input file not found: {self._params['input_file']}"
            )

        # Check if dir_for_links exists if specified
        if self._params["dir_for_links"] is not None:
            if not exists(self._params["dir_for_links"]):
                raise FileNotFoundError(
                    f"Directory for links not found: {self._params['dir_for_links']}"
                )
            if not os.path.isdir(self._params["dir_for_links"]):
                raise ValueError(
                    f"Path is not a directory: {self._params['dir_for_links']}"
                )

        # Create output directory if it doesn't exist
        if not exists(self._params["output_dir"]):
            os.makedirs(self._params["output_dir"])
            if self._params["verbose"]:
                print(f"Created output directory: {self._params['output_dir']}")

    def get_exposures(self, ccd_list_file):
        """Get Exposures.

        Extract unique exposure numbers from CCD list file.

        Parameters
        ----------
        ccd_list_file : str
            path to CCD list file (txt or csv)

        Returns
        -------
        np.array
            unique exposure numbers

        """
        exps = []

        # Check if CSV format
        if ccd_list_file.endswith(".csv"):
            import astropy.table

            t = astropy.table.Table.read(ccd_list_file)
            expccd = t["CCD"].data
            r = np.char.split(expccd, sep="-")

            for i, r1 in enumerate(r):
                exp = int(r1[0])
                exps.append(exp)
        else:
            # Text format
            f = np.loadtxt(ccd_list_file, dtype="str", encoding="ascii")
            r = np.char.split(f, sep="-")

            for i, r1 in enumerate(r):
                exp = int(r1[0])
                exps.append(exp)

        exps = np.array(exps)
        uniq = np.unique(exps)

        return uniq

    def get_fits_header(self, expnum, client):
        """Get FITS Header.

        Download FITS header from VOSpace, or create a symbolic link
        if the header exists in dir_for_links.

        Parameters
        ----------
        expnum : int
            exposure number
        client : vos.Client
            VOSpace client

        Returns
        -------
        bool
            True if successful, False otherwise

        """
        vospace_path = self._params["vospace_path"]
        output_dir = self._params["output_dir"]
        overwrite = self._params["overwrite"]
        dir_for_links = self._params["dir_for_links"]

        source = f"{vospace_path}/{expnum:d}p.fits.fz"
        dest = f"{output_dir}/{expnum:d}.txt"

        if exists(dest) and not overwrite:
            return True

        # Check if header exists in dir_for_links
        if dir_for_links is not None:
            link_source = os.path.abspath(f"{dir_for_links}/{expnum:d}.txt")
            if exists(link_source):
                try:
                    # Remove existing file/link if overwrite is True
                    if exists(dest):
                        os.remove(dest)
                    # Create symbolic link
                    os.symlink(link_source, dest)
                    return True
                except Exception as e:
                    print(f"Could not create symlink from {link_source}: {e}")
                    # Fall through to download if symlink fails

        # Download from VOSpace
        try:
            client.copy(source, dest, head=True)
            return True
        except Exception as e:
            print(f"Could not copy {source}: {e}")
            return False

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
        input_file = self._params["input_file"]
        output_dir = self._params["output_dir"]
        verbose = self._params["verbose"]

        if verbose:
            print(f"Reading CCD list from: {input_file}")

        # Extract unique exposures
        exps = self.get_exposures(input_file)

        print(f"Found {len(exps)} unique exposures")

        if verbose:
            print(f"Downloading headers to: {output_dir}/")

        # Initialize VOSpace client once
        client = vos.Client()

        # Download headers
        n_success = 0
        n_failed = 0

        for i, exp in enumerate(exps):
            success = self.get_fits_header(exp, client)
            if success:
                n_success += 1
            else:
                n_failed += 1

            if verbose and i % 100 == 0:
                print(f"{i:6d} / {len(exps):6d}")

        print(f"Downloaded {n_success} headers")
        if n_failed > 0:
            print(f"Failed to download {n_failed} headers")

        return 0
