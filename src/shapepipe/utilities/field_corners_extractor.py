"""FIELD_CORNERS_EXTRACTOR

Extract per-CCD sky-footprint corner coordinates from FITS headers.

For each CCD (image HDU) in an exposure's multi-HDU header, the four corners
of the CCD are projected from pixel bounds to the sky, producing one row per
CCD keyed on its ``<expnum>-<ccd_idx>`` ID.

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


def _image_shape(header):
    """Return the CCD image ``(nx, ny)`` pixel shape from a header.

    For fpack tile-compressed HDUs (``ZIMAGE = T``), ``NAXIS1``/``NAXIS2``
    describe the compressed *binary table* (row width in bytes, row count),
    not the image, and astropy's WCS never maps the true dimensions into
    ``pixel_shape``. The true image size is carried by ``ZNAXIS1``/``ZNAXIS2``,
    which are preferred here; a plain image HDU falls back to
    ``NAXIS1``/``NAXIS2``.

    Parameters
    ----------
    header : astropy.io.fits.Header
        single-HDU header

    Returns
    -------
    tuple
        ``(nx, ny)`` image pixel dimensions

    Raises
    ------
    ValueError
        if neither the ``Z*`` nor the plain ``NAXIS`` dimensions are present
    """
    if header.get("ZIMAGE", False):
        keys = ("ZNAXIS1", "ZNAXIS2")
    else:
        keys = ("NAXIS1", "NAXIS2")

    if keys[0] not in header or keys[1] not in header:
        raise ValueError(
            f"Header is missing image dimensions ({keys[0]}/{keys[1]})"
        )

    return int(header[keys[0]]), int(header[keys[1]])


def _parse_header_to_wcs(path):
    """Parse a multi-HDU header text file into ``(wcs, (nx, ny))`` per HDU.

    The primary HDU is skipped. Each remaining HDU yields its WCS together with
    the CCD image pixel shape (``ZNAXIS1/2`` for compressed HDUs, else
    ``NAXIS1/2``); the shape is read from the header because the WCS drops the
    ``Z*`` keywords.
    """
    with open(path, "r") as f:
        string = f.read()
    tokens = re.split(r"^(END\s+)", string, flags=re.MULTILINE)
    result = []
    for i in range(2, len(tokens) - 1, 2):
        header = Header.fromstring(tokens[i] + tokens[i + 1], sep="\n")
        result.append((wcs.WCS(header), _image_shape(header)))
    return result


def _ccd_corners(w, shape):
    """Return (ra, dec) of the 4 corners of a single CCD.

    The polygon is built from the CCD's pixel bounds ``shape = (nx, ny)`` using
    pixel edges ``(-0.5, n - 0.5)`` so the quadrilateral covers the full CCD
    area rather than pixel centres. Corners are returned in a consistent
    counterclockwise pixel order (bottom-left, bottom-right, top-right,
    top-left); since ``pixel_to_world`` maps this order to a simple,
    non-self-intersecting sky quadrilateral for a CCD-sized field, and
    HealSparse ``Polygon`` is orientation-agnostic, the resulting polygon is a
    valid convex footprint.

    Parameters
    ----------
    w : astropy.wcs.WCS
        WCS of a single CCD
    shape : tuple
        ``(nx, ny)`` CCD image pixel dimensions

    Returns
    -------
    tuple
        ``(ra_list, dec_list)`` each of length 4, in degrees
    """
    nx, ny = shape

    # Pixel-edge corners, counterclockwise: BL, BR, TR, TL
    px = [-0.5, nx - 0.5, nx - 0.5, -0.5]
    py = [-0.5, -0.5, ny - 0.5, ny - 0.5]

    sky = w.pixel_to_world(px, py)
    return list(sky.ra.deg), list(sky.dec.deg)


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
            "ccd_list": None,
            "resume": False,
            "n_processes": 1,
            "verbose": False,
        }

        self._short_options = {
            "input_dir": "-i",
            "output_file": "-o",
            "ccd_list": "-l",
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
            "output_file": "output file for per-CCD corners; default is {}",
            "ccd_list": "file of valid CCD IDs (output of get_ccds_with_psf); when given, only listed CCDs are written; on --resume, CCDs new to an expanded list are added without duplicating existing rows; default is all CCDs",
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

        if self._params["ccd_list"] is not None and not exists(
            self._params["ccd_list"]
        ):
            raise FileNotFoundError(
                f"CCD list file not found: {self._params['ccd_list']}"
            )

        # Set n_processes to cpu_count if 0
        if self._params["n_processes"] == 0:
            self._params["n_processes"] = cpu_count()

        if self._params["n_processes"] < 0:
            raise ValueError(
                f"n_processes must be >= 0, got {self._params['n_processes']}"
            )

    @staticmethod
    def load_ccd_list(path):
        """Load CCD List.

        Read valid CCD IDs (one ``<expnum>-<ccd_idx>`` per line) into a set.

        Parameters
        ----------
        path : str
            path to the CCD list file

        Returns
        -------
        set
            valid CCD IDs

        """
        with open(path) as f:
            return {line.strip() for line in f if line.strip()}

    @staticmethod
    def process_single_header(args):
        """Process Single Header.

        Worker function to process a single header file into per-CCD corners.
        Static method so it can be pickled for multiprocessing.

        Parameters
        ----------
        args : tuple
            ``(path, verbose, valid_ccds)`` where ``path`` is the header file
            path, ``verbose`` is a bool, and ``valid_ccds`` is a set of CCD IDs
            to keep (or ``None`` to keep all)

        Returns
        -------
        list or None
            list of ``(ccd_id, ra_list, dec_list)`` for the exposure's CCDs on
            success, ``None`` on failure

        """
        path, verbose, valid_ccds = args
        expnum = _expnum_from_path(path)

        try:
            wcs_shapes = _parse_header_to_wcs(path)
        except Exception as e:
            if verbose:
                print(f"Failed to process {expnum}: {e}")
            return None

        rows = []
        for ccd_idx, (w, shape) in enumerate(wcs_shapes):
            ccd_id = f"{expnum}-{ccd_idx}"
            if valid_ccds is not None and ccd_id not in valid_ccds:
                continue
            try:
                ra, dec = _ccd_corners(w, shape)
            except Exception as e:
                if verbose:
                    print(f"Failed to process CCD {ccd_id}: {e}")
                continue
            rows.append((ccd_id, ra, dec))

        return rows

    def get_done_ccds(self):
        """Get Done CCDs.

        Read the set of CCD IDs already present in the output file. Resume is
        keyed on individual CCD IDs, not exposure numbers: a CCD counts as done
        only if its own row is present. This keeps resume correct when a write
        was interrupted mid-exposure (the missing CCDs are filled in) and when
        a rerun uses an expanded ``--ccd_list`` (the newly requested CCDs are
        added), and it never duplicates a row.

        Returns
        -------
        set
            CCD IDs already written

        """
        output_file = self._params["output_file"]

        if not exists(output_file):
            return set()

        try:
            ids = np.atleast_1d(
                np.loadtxt(output_file, usecols=(0), dtype=str)
            )
            return set(ids.tolist())
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
            args = sys.argv[1:]

        # Set parameters from command line
        self.set_params_from_command_line(args)
        self.update_params()
        self.check_params()

        # Get parameters
        input_dir = self._params["input_dir"]
        output_file = self._params["output_file"]
        resume = self._params["resume"]
        verbose = self._params["verbose"]

        # Load the set of valid CCD IDs to keep, if given
        valid_ccds = None
        if self._params["ccd_list"] is not None:
            valid_ccds = self.load_ccd_list(self._params["ccd_list"])
            print(f"{len(valid_ccds)} valid CCDs in {self._params['ccd_list']}")

        # Find all header files
        paths = glob.glob(f"{input_dir}/*.txt")
        n = len(paths)

        if n == 0:
            print(f"No header files found in {input_dir}/")
            return 1

        print(f"{n} header files found")

        # On resume, read the CCD IDs already written so we can skip them
        # per-CCD (not per-exposure): a partially written exposure is completed
        # rather than skipped, and no row is ever duplicated.
        done_ccds = set()
        if resume:
            done_ccds = self.get_done_ccds()
            print(f"{len(done_ccds)} CCDs already done")

        # Every header still needs parsing (a header may hold both done and
        # not-yet-done CCDs); the per-CCD filter below decides what is written.
        todo = paths

        # Get n_processes
        n_processes = self._params["n_processes"]

        # Process headers
        if n_processes == 1:
            # Serial processing
            results = self._process_serial(todo, verbose, valid_ccds)
        else:
            # Parallel processing
            print(f"Using {n_processes} parallel processes")
            results = self._process_parallel(
                todo, n_processes, verbose, valid_ccds
            )

        # Flatten per-exposure CCD lists (dropping failed exposures), then drop
        # CCDs already present in the output.
        rows = [
            row
            for res in results if res is not None
            for row in res
            if row[0] not in done_ccds
        ]

        # Sort by exposure number, then CCD index
        rows.sort(key=lambda r: (int(r[0].split("-")[0]),
                                 int(r[0].split("-")[1])))

        # Write results to file: "<ccd_id> ra1 ra2 ra3 ra4 dec1 dec2 dec3 dec4"
        mode = "a" if resume else "w"
        with open(output_file, mode, buffering=1) as f:
            for ccd_id, ra, dec in rows:
                f.write(f"{ccd_id} ")
                np.savetxt(f, ra, fmt="%9.5f", newline=" ")
                np.savetxt(f, dec, fmt="%9.5f", newline=" ")
                f.write("\n")

        n_exp_success = sum(1 for res in results if res is not None)
        n_failed = len(todo) - n_exp_success

        print(f"Processed {n_exp_success} exposures, {len(rows)} new CCDs")
        if n_failed > 0:
            print(f"Failed to process {n_failed} exposures")

        print(f"Results written to {output_file}")

        return 0

    def _process_serial(self, todo, verbose, valid_ccds):
        """Process Serial.

        Process headers serially.

        Parameters
        ----------
        todo : list
            list of header file paths to process
        verbose : bool
            verbose output
        valid_ccds : set or None
            CCD IDs to keep, or ``None`` to keep all

        Returns
        -------
        list
            list of per-exposure CCD-row lists (or ``None`` for failures)

        """
        results = []
        n_todo = len(todo)

        for i, p in enumerate(todo):
            result = self.process_single_header((p, verbose, valid_ccds))
            results.append(result)

            if verbose and i % 100 == 0:
                print(f"{i:6d} / {n_todo:6d}")

        return results

    def _process_parallel(self, todo, n_processes, verbose, valid_ccds):
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
        valid_ccds : set or None
            CCD IDs to keep, or ``None`` to keep all

        Returns
        -------
        list
            list of per-exposure CCD-row lists (or ``None`` for failures)

        """
        # Prepare arguments for worker function
        args = [(p, verbose, valid_ccds) for p in todo]

        # Create pool and process
        with Pool(processes=n_processes) as pool:
            results = pool.map(self.process_single_header, args)

        return results
