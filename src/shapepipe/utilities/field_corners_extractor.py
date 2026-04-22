"""FIELD_CORNERS_EXTRACTOR

Extract field corner coordinates from FITS headers.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
import os
import glob
import re
from os.path import exists
from multiprocessing import Pool, cpu_count

import numpy as np
from astropy import wcs
from astropy.io.fits import Header

from cs_util import args as cs_args
from cs_util import logging

# Filename suffix convention for header files: <6+ digit exposure number>.txt
_EXPNUM_RE = re.compile(r'(\d+)\.txt$')


def _expnum_from_path(path):
    """Extract exposure number from a header filename like ``1234567.txt``."""
    match = _EXPNUM_RE.search(path)
    if match is None:
        raise ValueError(f"Could not extract exposure number from {path!r}")
    return int(match.group(1))


def _parse_header_to_wcs(path):
    """Parse a multi-HDU header text file and return one WCS per HDU."""
    with open(path, "r") as f:
        string = f.read()
    tokens = re.split(r"^(END\s+)", string, flags=re.MULTILINE)
    return [
        wcs.WCS(Header.fromstring(tokens[i] + tokens[i + 1], sep="\n"))
        for i in range(2, len(tokens) - 1, 2)
    ]


def _megacam_field_corners(w):
    """Return (ra, dec) of the 4 MegaCam field corners.

    The hard-coded CCD indices (0, 8, 35, 27) and pixel coords
    (2079/32, 0) are the MegaCam corner convention; CCDs are flipped so
    the nominal bottom-right CCD reports the top-left sky corner, etc.
    """
    tl = w[0].pixel_to_world(2079, 0)
    tr = w[8].pixel_to_world(32, 0)
    br = w[35].pixel_to_world(2079, 0)
    bl = w[27].pixel_to_world(32, 0)
    return (
        [tl.ra.deg, tr.ra.deg, br.ra.deg, bl.ra.deg],
        [tl.dec.deg, tr.dec.deg, br.dec.deg, bl.dec.deg],
    )


class FieldCornersExtractor(object):
    """Field Corners Extractor Class.

    Extracts RA/Dec coordinates of field corners from FITS headers.
    """

    def __init__(self):
        """Initialize the extractor."""
        self.params_default()

    def params_default(self):
        """Set default parameters and command line options."""

        self._params = {
            "input_dir": "header",
            "output_file": "exp_ra_dec.txt",
            "resume": False,
            "n_processes": 1,
            "verbose": False,
        }

        self._short_options = {
            "input_dir": "-i",
            "output_file": "-o",
            "resume": "-r",
            "n_processes": "-n",
        }

        self._types = {
            "resume": "bool",
            "n_processes": "int",
            "verbose": "bool",
        }

        self._help_strings = {
            "input_dir": "input directory containing header files; default is {}",
            "output_file": "output file for field corners; default is {}",
            "resume": "resume from existing output file; default is {}",
            "n_processes": f"number of parallel processes (1=serial, 0=auto={cpu_count()}); default is {{}}",
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
        # Ensure input directory ends without trailing slash
        if self._params["input_dir"].endswith("/"):
            self._params["input_dir"] = self._params["input_dir"][:-1]

    def check_params(self):
        """Check parameters for validity."""
        if not exists(self._params["input_dir"]):
            raise FileNotFoundError(
                f"Input directory not found: {self._params['input_dir']}"
            )

        # Set n_processes to cpu_count if 0
        if self._params["n_processes"] == 0:
            self._params["n_processes"] = cpu_count()

        if self._params["n_processes"] < 0:
            raise ValueError(
                f"n_processes must be >= 0, got {self._params['n_processes']}"
            )

    @staticmethod
    def process_single_header(path_and_verbose):
        """Process Single Header.

        Worker function to process a single header file.
        Static method so it can be pickled for multiprocessing.

        Parameters
        ----------
        path_and_verbose : tuple
            (path, verbose) where path is header file path and verbose is bool

        Returns
        -------
        tuple or None
            (expnum, ra_list, dec_list) on success, None on failure

        """
        path, verbose = path_and_verbose
        expnum = _expnum_from_path(path)

        try:
            w = _parse_header_to_wcs(path)
            ra, dec = _megacam_field_corners(w)
            return (expnum, ra, dec)
        except Exception as e:
            if verbose:
                print(f"Failed to process {expnum}: {e}")
            return None

    def get_done_exposures(self):
        """Get Done Exposures.

        Read exposures already processed from output file.

        Returns
        -------
        set
            set of exposure numbers already processed

        """
        output_file = self._params["output_file"]

        if not exists(output_file):
            return set()

        try:
            done = np.loadtxt(output_file, usecols=(0), dtype=int)
            return set(done)
        except Exception as e:
            if self._params["verbose"]:
                print(f"Could not read existing output file: {e}")
            return set()

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
        input_dir = self._params["input_dir"]
        output_file = self._params["output_file"]
        resume = self._params["resume"]
        verbose = self._params["verbose"]

        # Find all header files
        paths = glob.glob(f"{input_dir}/*.txt")
        n = len(paths)

        if n == 0:
            print(f"No header files found in {input_dir}/")
            return 1

        print(f"{n} header files found")

        # Get list of already processed exposures if resuming
        done = set()
        if resume:
            done = self.get_done_exposures()
            print(f"{len(done)} already done")

        # Build todo list
        todo = [p for p in paths if _expnum_from_path(p) not in done]

        n_todo = len(todo)
        print(f"{n_todo} to process")

        if n_todo == 0:
            print("Nothing to do")
            return 0

        # Get n_processes
        n_processes = self._params["n_processes"]

        # Process headers
        if n_processes == 1:
            # Serial processing
            results = self._process_serial(todo, verbose)
        else:
            # Parallel processing
            print(f"Using {n_processes} parallel processes")
            results = self._process_parallel(todo, n_processes, verbose)

        # Filter out failed results (None values)
        results = [r for r in results if r is not None]

        # Sort results by exposure number
        results.sort(key=lambda x: x[0])

        # Write results to file
        mode = "a" if resume else "w"
        with open(output_file, mode, buffering=1) as f:
            for expnum, ra, dec in results:
                f.write(f"{expnum:6d} ")
                np.savetxt(f, ra, fmt="%9.5f", newline=" ")
                np.savetxt(f, dec, fmt="%9.5f", newline=" ")
                f.write("\n")

        n_success = len(results)
        n_failed = n_todo - n_success

        print(f"Processed {n_success} exposures")
        if n_failed > 0:
            print(f"Failed to process {n_failed} exposures")

        print(f"Results written to {output_file}")

        return 0

    def _process_serial(self, todo, verbose):
        """Process Serial.

        Process headers serially.

        Parameters
        ----------
        todo : list
            list of header file paths to process
        verbose : bool
            verbose output

        Returns
        -------
        list
            list of results (expnum, ra, dec) tuples

        """
        results = []
        n_todo = len(todo)

        for i, p in enumerate(todo):
            result = self.process_single_header((p, verbose))
            results.append(result)

            if verbose and i % 100 == 0:
                print(f"{i:6d} / {n_todo:6d}")

        return results

    def _process_parallel(self, todo, n_processes, verbose):
        """Process Parallel.

        Process headers in parallel using multiprocessing.

        Parameters
        ----------
        todo : list
            list of header file paths to process
        n_processes : int
            number of parallel processes
        verbose : bool
            verbose output

        Returns
        -------
        list
            list of results (expnum, ra, dec) tuples

        """
        # Prepare arguments for worker function
        args = [(p, verbose) for p in todo]

        # Create pool and process
        with Pool(processes=n_processes) as pool:
            results = pool.map(self.process_single_header, args)

        return results
