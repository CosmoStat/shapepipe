#! /usr/bin/env python3

import sys
import os
import re
import glob
import warnings
from tqdm import tqdm

# Suppress astropy WCS warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings("ignore", category=FITSFixedWarning)
from joblib import Parallel, delayed

import numpy as np
from astropy.io import fits
import galsim

from cs_util import args as cs_args
from cs_util import logging


def transform_shape(mom_list, jac):
    """Transform Shape.

    Transform shape (ellipticity and size) using a Jacobian.

    Parameters
    ----------
    mom_list : list
        input moment measurements; each list element contains
        first and second ellipticity component and size
    jac : galsim.JacobianWCS
        Jacobian transformation matrix information

    Returns
    -------
    list
        transformed shape parameters, which are
        first and second ellipticity component and size

    """
    scale, shear, theta, flip = jac.getDecomposition()

    sig_tmp = mom_list[2] * scale
    shape = galsim.Shear(g1=mom_list[0], g2=mom_list[1])
    if flip:
        # The following output is not observed
        print("FLIP!")
        shape = galsim.Shear(g1=-shape.g1, g2=shape.g2)
    shape = galsim.Shear(g=shape.g, beta=shape.beta + theta)
    shape = shear + shape

    return shape.g1, shape.g2, sig_tmp


def transform_shape_vectorized(e1, e2, sigma, scales, shear_g1, shear_g2, thetas, flips):
    """Transform Shape (Vectorized).

    Transform shapes (ellipticity and size) using Jacobian decomposition arrays.

    Parameters
    ----------
    e1 : np.ndarray
        First ellipticity component for each object
    e2 : np.ndarray
        Second ellipticity component for each object
    sigma : np.ndarray
        Size for each object
    scales : np.ndarray
        Scale factor from Jacobian decomposition for each object
    shear_g1 : np.ndarray
        First shear component from Jacobian for each object
    shear_g2 : np.ndarray
        Second shear component from Jacobian for each object
    thetas : np.ndarray
        Rotation angle from Jacobian for each object
    flips : np.ndarray
        Boolean flip flag from Jacobian for each object

    Returns
    -------
    tuple of np.ndarray
        Transformed (e1, e2, sigma) arrays

    """
    n = len(e1)
    new_e1 = np.zeros(n)
    new_e2 = np.zeros(n)
    new_sigma = sigma * scales

    for i in range(n):
        shape = galsim.Shear(g1=e1[i], g2=e2[i])
        if flips[i]:
            shape = galsim.Shear(g1=-shape.g1, g2=shape.g2)
        shape = galsim.Shear(g=shape.g, beta=shape.beta + thetas[i] * galsim.radians)
        shear = galsim.Shear(g1=shear_g1[i], g2=shear_g2[i])
        shape = shear + shape
        new_e1[i] = shape.g1
        new_e2[i] = shape.g2

    return new_e1, new_e2, new_sigma


class Loc2Glob(object):
    r"""Change from local to global coordinates.

    Class to pass from local coordinates to global coordinates under
    CFIS (CFHT) MegaCam instrument. The geometrical informcation of the
    instrument is encoded in this function.

    Parameters
    ----------
    x_gap : int
        Gap between the CCDs along the horizontal direction;
        default is ``70`` (MegaCam value)
    y_gap : int
        Gap between the CCDs along the vertical direction;
        Default is ``425`` (MegaCam value)
    x_npix : int
        Number of pixels per CCD along the horizontal direction;
        default is ``2048`` (MegaCam value)
    y_npix : int
        Number of pixels per CCD along the vertical direction;
        default to ``4612`` (MegaCam value)
    ccd_tot : int
        Total number of CCDs;
        default to ``40`` (MegaCam value)

    Notes
    -----
    This is the geometry of MegaCam. Watch out with the conventions ba,ab that means where
    is the local coordinate system origin for each CCD.
    For more info check out MegaCam's instrument webpage.

    Examples
    --------
        'COMMENT (North on top, East to the left)',
        'COMMENT    --------------------------',
        'COMMENT    ba ba ba ba ba ba ba ba ba',
        'COMMENT    00 01 02 03 04 05 06 07 08',
        'COMMENT --------------------------------',
        'COMMENT ba ba ba ba ba ba ba ba ba ba ba',
        'COMMENT 36 09 10 11 12 13 14 15 16 17 37',
        'COMMENT --------------*-----------------',
        'COMMENT 38 18 19 20 21 22 23 24 25 26 39',
        'COMMENT ab ab ab ab ab ab ab ab ab ab ab',
        'COMMENT --------------------------------',
        'COMMENT    27 28 29 30 31 32 33 34 35',
        'COMMENT    ab ab ab ab ab ab ab ab ab',
        'COMMENT    __________________________'
    """

    def __init__(
        self, x_gap=70, y_gap=425, x_npix=2048, y_npix=4612, ccd_tot=40
    ):
        r"""Initialize with instrument geometry."""
        self.x_gap = x_gap
        self.y_gap = y_gap
        self.x_npix = x_npix
        self.y_npix = y_npix
        self.ccd_tot = ccd_tot

    def loc2glob_img_coord(self, ccd_n, x_coor, y_coor):
        """loc2glob Img Coord.

        Go from the local to the global img (pixel) coordinate system.

        Global system with (0,0) in the intersection of ccds [12,13,21,22].

        Parameters
        ----------
        ccd_n: int
            CCD number of the considered positions
        x_coor: float
            Local coordinate system hotizontal value
        y_coor: float
            Local coordinate system vertical value

        Returns
        -------
        glob_x_coor: float
            Horizontal position in global coordinate system
        glob_y_coor: float
            Vertical position in global coordinate system

        """
        # Flip axes
        x_coor, y_coor = self.flip_coord(ccd_n, x_coor, y_coor)

        # Calculate the shift
        x_shift, y_shift = self.shift_coord(ccd_n)

        # Return new coordinates
        return x_coor + x_shift, y_coor + y_shift

    def flip_coord(self, ccd_n, x_coor, y_coor):
        r"""Change of coordinate convention.

        So that all of them are coherent on the global coordinate system.
        So that the origin is on the south-west corner.
        Positive: South to North ; West to East.
        """
        if ccd_n < 18 or ccd_n in [36, 37]:
            x_coor = self.x_npix - x_coor + 1
            y_coor = self.y_npix - y_coor + 1
        else:
            pass

        return x_coor, y_coor

    def x_coord_range(self):
        r"""Return range of the x coordinate."""
        max_x = self.x_npix * 6 + self.x_gap * 5
        min_x = self.x_npix * (-5) + self.x_gap * (-5)
        return min_x, max_x

    def y_coord_range(self):
        r"""Return range of the y coordinate."""
        max_y = self.y_npix * 2 + self.y_gap * 1
        min_y = self.y_npix * (-2) + self.y_gap * (-2)
        return min_y, max_y

    def shift_coord(self, ccd_n):
        r"""Provide the shifting.

        It is needed to go from the local coordinate
        system origin to the global coordinate system origin.
        """
        if ccd_n < 9:
            # first row
            x_shift = (ccd_n - 4) * (self.x_gap + self.x_npix)
            y_shift = self.y_gap + self.y_npix
            return x_shift, y_shift

        elif ccd_n < 18:
            # second row, non-ears
            x_shift = (ccd_n - 13) * (self.x_gap + self.x_npix)
            y_shift = 0.0
            return x_shift, y_shift

        elif ccd_n < 27:
            # third row non-ears
            x_shift = (ccd_n - 22) * (self.x_gap + self.x_npix)
            y_shift = -1.0 * (self.y_gap + self.y_npix)
            return x_shift, y_shift

        elif ccd_n < 36:
            # fourth row
            x_shift = (ccd_n - 31) * (self.x_gap + self.x_npix)
            y_shift = -2.0 * (self.y_gap + self.y_npix)
            return x_shift, y_shift

        elif ccd_n < 37:
            # ccd= 36 ears, second row
            x_shift = (-5.0) * (self.x_gap + self.x_npix)
            y_shift = 0.0
            return x_shift, y_shift

        elif ccd_n < 38:
            # ccd= 37 ears, second row
            x_shift = 5.0 * (self.x_gap + self.x_npix)
            y_shift = 0.0
            return x_shift, y_shift

        elif ccd_n < 39:
            # ccd= 38 ears, third row
            x_shift = (-5.0) * (self.x_gap + self.x_npix)
            y_shift = -1.0 * (self.y_gap + self.y_npix)
            return x_shift, y_shift

        elif ccd_n < 40:
            # ccd= 39 ears, third row
            x_shift = 5.0 * (self.x_gap + self.x_npix)
            y_shift = -1.0 * (self.y_gap + self.y_npix)
            return x_shift, y_shift


class Glob2CCD(object):
    r"""Get the CCD ID number from the global coordinate position.

    The Loc2Glob() object as input is the one that defines the instrument's
    geometry.

    Parameters
    ----------
    loc2glob: Loc2Glob object
        Object with the desired focal plane geometry.
    with_gaps: bool
        If add the gaps to the CCD area.
    """

    def __init__(self, loc2glob, with_gaps=True):
        # Save loc2glob object
        self.loc2glob = loc2glob
        self.with_gaps = with_gaps
        self.ccd_list = np.arange(self.loc2glob.ccd_tot)
        # Init edges defininf the CCDs
        self.edge_x_list, self.edge_y_list = self.build_all_edges()

    def build_all_edges(self):
        """Build the edges for all the CCDs in the focal plane."""
        edge_xy_list = []
        for idx in (0, 1):
            edge_list = np.array(
                [self.build_edge(ccd_n)[idx] for ccd_n in self.ccd_list]
            )
            edge_xy_list.append(edge_list)

        return edge_xy_list

    def build_edge(self, ccd_n):
        """Build the edges of the `ccd_n` in global coordinates."""
        if self.with_gaps:
            corners = np.array(
                [
                    [-self.loc2glob.x_gap / 2, -self.loc2glob.y_gap / 2],
                    [
                        self.loc2glob.x_npix + self.loc2glob.x_gap / 2,
                        -self.loc2glob.y_gap / 2,
                    ],
                    [
                        -self.loc2glob.x_gap / 2,
                        self.loc2glob.y_npix + self.loc2glob.y_gap / 2,
                    ],
                    [
                        self.loc2glob.x_npix + self.loc2glob.x_gap / 2,
                        self.loc2glob.y_npix + self.loc2glob.y_gap / 2,
                    ],
                ]
            )
        else:
            corners = np.array(
                [
                    [0, 0],
                    [self.loc2glob.x_npix, 0],
                    [0, self.loc2glob.y_npix],
                    [self.loc2glob.x_npix, self.loc2glob.y_npix],
                ]
            )

        glob_corners = np.array(
            [
                self.loc2glob.loc2glob_img_coord(ccd_n, pos[0], pos[1])
                for pos in corners
            ]
        )

        edge_xy = []
        for idx in (0, 1):
            edge = np.array(
                [np.min(glob_corners[:, idx]), np.max(glob_corners[:, idx])]
            )
            edge_xy.append(edge)

        return edge_xy

    def is_inside(self, x, y, edge_x, edge_y):
        """Is the position inside the edges.

        Return True if the position is within the rectangle
        defined by the edges.

        Parameters
        ----------
        x: float
            Horizontal position in global coordinate system.
        y: float
            Vertical position in global coordinate system.
        edge_x: np.ndarray
            Edge defined as `np.array([min_x, max_x])`.
        edge_y: np.ndarray
            Edge defined as `np.array([min_y, max_y])`.
        """
        if (
            (x > edge_x[0])
            and (x < edge_x[1])
            and (y > edge_y[0])
            and (y < edge_y[1])
        ):
            return True
        else:
            return False

    def get_ccd_n(self, x, y):
        """Returns the CCD number from the position `(x, y)`.

        Returns `None` if the position is not found.
        """
        bool_list = np.array(
            [
                self.is_inside(x, y, edge_x, edge_y)
                for edge_x, edge_y in zip(self.edge_x_list, self.edge_y_list)
            ]
        )

        try:
            return self.ccd_list[bool_list][0]
        except Exception:
            return None

    def get_ccd_n_vectorized(self, x, y):
        """Returns CCD numbers for arrays of positions (x, y).

        Vectorized version for processing many positions at once.

        Parameters
        ----------
        x : np.ndarray
            Array of horizontal positions in global coordinate system.
        y : np.ndarray
            Array of vertical positions in global coordinate system.

        Returns
        -------
        np.ndarray
            Array of CCD numbers. Returns -1 for positions not found.
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

        # Shape: (n_points, n_ccds)
        in_x = (x[:, None] > self.edge_x_list[:, 0]) & (
            x[:, None] < self.edge_x_list[:, 1]
        )
        in_y = (y[:, None] > self.edge_y_list[:, 0]) & (
            y[:, None] < self.edge_y_list[:, 1]
        )
        inside = in_x & in_y

        # Return first matching CCD for each point (-1 if none found)
        ccd_ids = np.where(inside.any(axis=1), inside.argmax(axis=1), -1)
        return ccd_ids


class Convert(object):

    def __init__(self):

        self.params_default()

    def set_params_from_command_line(self, args):
        """Set Params From Command line.

        Only use when calling using python from command line.
        Does not work from ipython or jupyter.

        """
        # Read command line options
        options = cs_args.parse_options(
            self._params,
            self._short_options,
            self._types,
            self._help_strings,
        )
        self._params = options

        # Save calling command
        logging.log_command(args)

    def params_default(self):

        self._params = {
            "input_base_dir": ".",
            "output_base_dir": ".",
            "mode": "merge",
            "patches": "",
            "psf": "psfex",
            "file_pattern_psfint": "validation_psf",
            "n_jobs": 8,
        }

        self._short_options = {
            "input_base_dir": "-i",
            "mode": "-m",
            "psf": "-p",
            "patches": "-P",
            "n_jobs": "-j",
        }

        self._types = {
            "n_jobs": "int",
        }

        self._help_strings = {
            "input_base_dir": (
                "input base dir, runs are expected in"
                + " <input_base_dir>/P<patch?>/output;"
                " default is {}"
            ),
            "mode": (
                "run mode, allowed are 'merge', 'test'; default is" + " '{}'"
            ),
            "psf": "PSF model, allowed are 'psfex' and 'mccd'; default is {}",
            "patches": "(list of) input patches",
            "n_jobs": "number of parallel jobs; default is {}",
        }

        # Output column names with types
        self._dt = [
            ("X", float),
            ("Y", float),
            ("RA", float),
            ("DEC", float),
            ("E1_PSF_HSM", float),
            ("E2_PSF_HSM", float),
            ("SIGMA_PSF_HSM", float),
            ("FLAG_PSF_HSM", float),
            ("E1_STAR_HSM", float),
            ("E2_STAR_HSM", float),
            ("SIGMA_STAR_HSM", float),
            ("FLAG_STAR_HSM", float),
            ("CCD_NB", int),
        ]

        # Extra columns for MCCD:737
        self._dt_mccd = self._dt.copy()
        self._dt_mccd.append(("GLOB_X", float))
        self._dt_mccd.append(("GLOB_Y", float))

    def update_params(self):
        """Update Params.

        Update parameters.

        """
        if self._params["psf"] == "psfex":
            #self._params["sub_dir_pattern"] = "run_sp_exp_202"
            self._params["sub_dir_pattern"] = "run_sp_combined_psf"
            self._params["sub_dir_psfint"] = "psfex_interp_runner"
        elif self._params["psf"] == "mccd":
            self._params["sub_dir_pattern"] = "run_sp_exp_SxSePsf_202"
            self._params["sub_dir_psfint"] = "mccd_fit_val_runner"
            self._params["sub_dir_setools"] = "setools_runner/output/mask"
            self._params["file_pattern_setools"] = "star_selection"
        else:
            raise ValueError(f"Invalid PSF model {self._params['psf']}")
        self._params["sub_dir_psfint"] = (
            f"{self._params['sub_dir_psfint']}/output"
        )

    def run(self):
        """Run.

        Main processing function.

        """
        if self._params["mode"] == "test":
            patch_nums = ["3", "4"]
        else:
            patch_nums = cs_args.my_string_split(self._params["patches"])

        n_jobs = self._params["n_jobs"]
        do_parallel = n_jobs != 1
        print(f"n_jobs={n_jobs}, patches={patch_nums}")

        for patch in patch_nums:

            output_dir = f"{self._params['output_base_dir']}/P{patch}"
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir, exist_ok=False)

            patch_dir = f"{self._params['input_base_dir']}/P{patch}/output/"
            subdirs = f"{patch_dir}/{self._params['sub_dir_pattern']}*"
            exp_run_dirs = glob.glob(subdirs)
            exp_run_dirs.sort()
            n_exp_runs = len(exp_run_dirs)
            print(f"Patch {patch}: {n_exp_runs} exposure runs in {patch_dir}")

            if self._params["mode"] == "test":
                exp_run_dirs = exp_run_dirs[:2]
                n_exp_runs = len(exp_run_dirs)

            # Loop over exposure runs
            if not do_parallel:
                for idx_exp, exp_run_dir in tqdm(
                    enumerate(exp_run_dirs),
                    total=n_exp_runs,
                    desc=f"Patch {patch}",
                ):
                    self.transform_exposures(
                        output_dir, patch, idx_exp, exp_run_dir
                    )
            else:
                res = Parallel(n_jobs=n_jobs, backend="threading")(
                    delayed(self.transform_exposures)(
                        output_dir, patch, idx_exp, exp_run_dir
                    )
                    for idx_exp, exp_run_dir in tqdm(
                        enumerate(exp_run_dirs),
                        total=n_exp_runs,
                        desc=f"Patch {patch}",
                    )
                )

    def transform_exposures(self, output_dir, patch, idx, exp_run_dir):
        """Transform exposures.

        Transform shapes for exposure for a given run (input exp run dir).

        """
        psf_dir = f"{exp_run_dir}/{self._params['sub_dir_psfint']}"
        try:
            all_files = os.listdir(psf_dir)
        except Exception:
            return

        for file_name in all_files:
            if self._params["file_pattern_psfint"] not in file_name:
                continue

            tmp = re.findall(r"\d+", file_name)

            if self._params["psf"] == "psfex":
                exp_name, ccd_id = int(tmp[0]), int(tmp[1])

                # Output path per input file, encoding exposure and CCD
                output_path = (
                    f"{output_dir}/{self._params['file_pattern_psfint']}_conv"
                    f"-{exp_name}-{ccd_id}.fits"
                )

            elif self._params["psf"] == "mccd":
                exp_name = int(tmp[0])
                ccd_id = "all"

                # Output path per input file, encoding exposure (all CCDs in
                # one file)
                output_path = (
                    f"{output_dir}/{self._params['file_pattern_psfint']}_conv"
                    f"-{exp_name}.fits"
            )
            if os.path.exists(output_path):
                continue

            psf_file_path = f"{psf_dir}/{file_name}"

            try:
                if self._params["psf"] == "psfex":
                    psf_file_hdus = fits.open(psf_file_path, memmap=False)
                    psf_file = psf_file_hdus[2].data
                    header_file = psf_file_hdus[1].data
                    psf_file_hdus.close()
                    mod = "RA"
                else:
                    psf_file = fits.getdata(psf_file_path, 1, memmap=True)
                    mod = "RA_LIST"
            except Exception:
                continue

            new_e1_psf = np.zeros_like(psf_file[mod])
            new_e2_psf = np.zeros_like(psf_file[mod])
            new_sig_psf = np.zeros_like(psf_file[mod])
            new_e1_star = np.zeros_like(psf_file[mod])
            new_e2_star = np.zeros_like(psf_file[mod])
            new_sig_star = np.zeros_like(psf_file[mod])

            if self._params["psf"] == "psfex":
                header = fits.Header.fromstring(
                    "\n".join(header_file[0][0]), sep="\n"
                )
                wcs = galsim.AstropyWCS(header=header)

                k = 0
                for ind, obj in enumerate(psf_file):
                    try:
                        jac = wcs.jacobian(
                            image_pos=galsim.PositionD(
                                obj["X"],
                                obj["Y"],
                            )
                        )
                    except Exception:
                        continue
                    g1_psf_tmp, g2_psf_tmp, sig_psf_tmp = transform_shape(
                        [
                            obj["E1_PSF_HSM"],
                            obj["E2_PSF_HSM"],
                            obj["SIGMA_PSF_HSM"],
                        ],
                        jac,
                    )

                    new_e1_psf[ind] = g1_psf_tmp
                    new_e2_psf[ind] = g2_psf_tmp
                    new_sig_psf[ind] = sig_psf_tmp

                    g1_star_tmp, g2_star_tmp, sig_star_tmp = transform_shape(
                        [
                            obj["E1_STAR_HSM"],
                            obj["E2_STAR_HSM"],
                            obj["SIGMA_STAR_HSM"],
                        ],
                        jac,
                    )
                    new_e1_star[ind] = g1_star_tmp
                    new_e2_star[ind] = g2_star_tmp
                    new_sig_star[ind] = sig_star_tmp
                    k += 1

                # Build output catalogue using direct structured array assignment
                n_obj = len(psf_file)
                exp_cat = np.empty(n_obj, dtype=self._dt)
                exp_cat["X"] = psf_file["X"]
                exp_cat["Y"] = psf_file["Y"]
                exp_cat["RA"] = psf_file["RA"]
                exp_cat["DEC"] = psf_file["DEC"]
                exp_cat["E1_PSF_HSM"] = new_e1_psf
                exp_cat["E2_PSF_HSM"] = new_e2_psf
                exp_cat["SIGMA_PSF_HSM"] = new_sig_psf
                exp_cat["FLAG_PSF_HSM"] = psf_file["FLAG_PSF_HSM"]
                exp_cat["E1_STAR_HSM"] = new_e1_star
                exp_cat["E2_STAR_HSM"] = new_e2_star
                exp_cat["SIGMA_STAR_HSM"] = new_sig_star
                exp_cat["FLAG_STAR_HSM"] = psf_file["FLAG_STAR_HSM"]
                exp_cat["CCD_NB"] = ccd_id

                # Write output file immediately
                hdul = fits.HDUList()
                hdul.append(fits.PrimaryHDU())
                hdul.append(fits.BinTableHDU(exp_cat))
                try:
                    hdul.writeto(output_path, overwrite=True)
                except (FileNotFoundError, OSError):
                    pass  # Race condition with another thread
                hdul.close()

            else:
                # MCCD processing
                l2g = Loc2Glob()
                g2c = Glob2CCD(l2g)

                n_obj = len(psf_file)
                glob_x = psf_file["GLOB_POSITION_IMG_LIST"][:, 0]
                glob_y = psf_file["GLOB_POSITION_IMG_LIST"][:, 1]

                # Vectorized CCD lookup (replaces per-object loop)
                new_ccd_id = g2c.get_ccd_n_vectorized(glob_x, glob_y)

                # Pre-compute shifts for all CCDs (done once)
                x_shift = np.zeros(40)
                y_shift = np.zeros(40)
                for ccd_id in range(40):
                    x_shift[ccd_id], y_shift[ccd_id] = l2g.shift_coord(ccd_id)

                # Vectorized computation of local (x, y) from global coords
                # Handle invalid CCD IDs (-1) by using CCD 0 temporarily
                safe_ccd_id = np.where(new_ccd_id >= 0, new_ccd_id, 0)
                new_x = glob_x - x_shift[safe_ccd_id]
                new_y = glob_y - y_shift[safe_ccd_id]

                # Pre-allocate output arrays
                new_e1_psf = np.zeros(n_obj)
                new_e2_psf = np.zeros(n_obj)
                new_sig_psf = np.zeros(n_obj)
                new_e1_star = np.zeros(n_obj)
                new_e2_star = np.zeros(n_obj)
                new_sig_star = np.zeros(n_obj)
                new_flag_psf = np.zeros(n_obj)
                new_flag_star = np.zeros(n_obj)

                # Process by CCD (only load WCS for CCDs that have objects)
                unique_ccds = np.unique(new_ccd_id)
                setools_dir = f"{exp_run_dir}/{self._params['sub_dir_setools']}"
                for ccd_id in unique_ccds:
                    if ccd_id < 0:
                        # Mark objects with invalid CCD as flagged
                        mask = new_ccd_id == ccd_id
                        new_flag_psf[mask] = 16
                        new_flag_star[mask] = 16
                        continue

                    # Load WCS only for this CCD (lazy loading)
                    header_file_path = (
                        f"{setools_dir}/{self._params['file_pattern_setools']}"
                        f"-{exp_name}-{ccd_id}.fits"
                    )
                    try:
                        header_file = fits.getdata(header_file_path, 1)
                        header = fits.Header.fromstring(
                            "\n".join(header_file[0][0]), sep="\n"
                        )
                        this_wcs = galsim.AstropyWCS(header=header)
                    except Exception:
                        mask = new_ccd_id == ccd_id
                        new_flag_psf[mask] = 16
                        new_flag_star[mask] = 16
                        continue

                    mask = new_ccd_id == ccd_id
                    indices = np.where(mask)[0]
                    n_ccd_obj = len(indices)

                    # Extract data for this CCD's objects
                    ra_batch = psf_file["RA_LIST"][mask]
                    dec_batch = psf_file["DEC_LIST"][mask]
                    psf_mom_batch = psf_file["PSF_MOM_LIST"][mask]
                    star_mom_batch = psf_file["STAR_MOM_LIST"][mask]

                    # Pre-allocate Jacobian decomposition arrays
                    scales = np.zeros(n_ccd_obj)
                    shear_g1 = np.zeros(n_ccd_obj)
                    shear_g2 = np.zeros(n_ccd_obj)
                    thetas = np.zeros(n_ccd_obj)
                    flips = np.zeros(n_ccd_obj, dtype=bool)
                    valid = np.ones(n_ccd_obj, dtype=bool)

                    # Compute Jacobians for all objects in this CCD
                    for i in range(n_ccd_obj):
                        try:
                            jac = this_wcs.jacobian(
                                world_pos=galsim.CelestialCoord(
                                    ra=ra_batch[i] * galsim.degrees,
                                    dec=dec_batch[i] * galsim.degrees,
                                )
                            )
                            scale, shear, theta, flip = jac.getDecomposition()
                            scales[i] = scale
                            shear_g1[i] = shear.g1
                            shear_g2[i] = shear.g2
                            thetas[i] = theta.rad
                            flips[i] = flip
                        except Exception:
                            valid[i] = False

                    # Apply vectorized shape transform for valid objects
                    if valid.any():
                        # PSF shapes
                        e1_psf_v, e2_psf_v, sig_psf_v = transform_shape_vectorized(
                            psf_mom_batch[valid, 0],
                            psf_mom_batch[valid, 1],
                            psf_mom_batch[valid, 2],
                            scales[valid],
                            shear_g1[valid],
                            shear_g2[valid],
                            thetas[valid],
                            flips[valid],
                        )
                        # Star shapes
                        e1_star_v, e2_star_v, sig_star_v = transform_shape_vectorized(
                            star_mom_batch[valid, 0],
                            star_mom_batch[valid, 1],
                            star_mom_batch[valid, 2],
                            scales[valid],
                            shear_g1[valid],
                            shear_g2[valid],
                            thetas[valid],
                            flips[valid],
                        )

                        # Store results back using the original indices
                        valid_indices = indices[valid]
                        new_e1_psf[valid_indices] = e1_psf_v
                        new_e2_psf[valid_indices] = e2_psf_v
                        new_sig_psf[valid_indices] = sig_psf_v
                        new_flag_psf[valid_indices] = psf_mom_batch[valid, 3]
                        new_e1_star[valid_indices] = e1_star_v
                        new_e2_star[valid_indices] = e2_star_v
                        new_sig_star[valid_indices] = sig_star_v
                        new_flag_star[valid_indices] = star_mom_batch[valid, 3]

                    # Flag invalid objects
                    if not valid.all():
                        invalid_indices = indices[~valid]
                        new_flag_psf[invalid_indices] = 16
                        new_flag_star[invalid_indices] = 16

                # Build output catalogue using direct structured array assignment
                exp_cat = np.empty(n_obj, dtype=self._dt_mccd)
                exp_cat["X"] = new_x
                exp_cat["Y"] = new_y
                exp_cat["RA"] = psf_file["RA_LIST"]
                exp_cat["DEC"] = psf_file["DEC_LIST"]
                exp_cat["E1_PSF_HSM"] = new_e1_psf
                exp_cat["E2_PSF_HSM"] = new_e2_psf
                exp_cat["SIGMA_PSF_HSM"] = new_sig_psf
                exp_cat["FLAG_PSF_HSM"] = new_flag_psf
                exp_cat["E1_STAR_HSM"] = new_e1_star
                exp_cat["E2_STAR_HSM"] = new_e2_star
                exp_cat["SIGMA_STAR_HSM"] = new_sig_star
                exp_cat["FLAG_STAR_HSM"] = new_flag_star
                exp_cat["CCD_NB"] = new_ccd_id
                exp_cat["GLOB_X"] = glob_x
                exp_cat["GLOB_Y"] = glob_y

                # Write output file immediately
                hdul = fits.HDUList()
                hdul.append(fits.PrimaryHDU())
                hdul.append(fits.BinTableHDU(exp_cat))
                try:
                    hdul.writeto(output_path, overwrite=True)
                except (FileNotFoundError, OSError):
                    pass  # Race condition with another thread
                hdul.close()

            del psf_file


def run_convert(*args):

    # Create instance
    obj = Convert()

    obj.set_params_from_command_line(args)
    obj.update_params()

    obj.run()


def main(argv=None):
    """Main

    Main program
    """
    if argv is None:
        argv = sys.argv[1:]
    run_convert(*argv)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
