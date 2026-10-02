"""NGMIX.

This module contains a class for ngmix shape measurement.

:Authors: Lucie Baumont, Axel Guinot

"""

import os
import re
from collections import Counter
from typing import NamedTuple

import ngmix
import galsim
import numpy as np
from astropy.io import fits
from cs_util import size as cs_size
from modopt.math.stats import sigma_mad
from ngmix.observation import Observation, ObsList
from scipy.ndimage import binary_dilation
from scipy.spatial import cKDTree
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package.defect_interpolation import (
    fourfold,
    interpolable_defects,
    interpolate_defects,
)
from shapepipe.pipeline import file_io

# Neighbour treatments selectable with the BLEND_HANDLING option.
BLEND_HANDLINGS = ("noisefill", "uberseg")

# The epoch cuts (see :func:`prepare_postage_stamps`). Calibration outputs,
# not options: an epoch is dropped when more than EPOCH_MASKED_FRACTION_CUT
# of its stamp is in :func:`defect_mask`, or when a noise-filled
# (interpolated) defect pixel lies closer than EPOCH_CENTRAL_DEFECT_RADIUS
# (EPOCH_INTERPOLATED_DEFECT_RADIUS) pixels to the stamp centre (see
# :func:`central_defect_vetoes`).
# @sc [decision:shape_measurement.epoch_masked_fraction_cut]
EPOCH_MASKED_FRACTION_CUT = 1 / 3
# @sc [decision:shape_measurement.central_defect_veto]
EPOCH_CENTRAL_DEFECT_RADIUS = 10
# @sc [decision:shape_measurement.central_defect_veto]
EPOCH_INTERPOLATED_DEFECT_RADIUS = 7

# @sc [decision:shape_measurement.metacal_scheme]
METACAL_TYPES = ('noshear', '1p', '1m', '2p', '2m')

# Noise budget for the PSF observation's flat weight map (psf_wt =
# 1/PSF_NOISE**2). Mirrors the esheldon/aguinot pattern (sigma ~ 1e-5/1e-6);
# 1e-5 is the value Axel Guinot's #749 reproduction used. The fit is driven
# by the *relative* weighting of the likelihood vs the prior, so the exact
# value is non-critical once it is finite (validated on the digital twin: the
# recovered PSF shape/size are flat across 1e-4..1e-6). See
# make_ngmix_observation.
# @sc [decision:shape_measurement.psf_likelihood_noise]
PSF_NOISE = 1e-5


class MetacalResult(NamedTuple):
    """Return of :func:`do_ngmix_metacal`: the metacal fit plus both PSFs.

    A NamedTuple so the two PSF families are named, not positional — guarding
    against a reconv/orig transposition at call sites. Still a plain tuple, so
    ``resdict, psf_res, psf_orig_res = do_ngmix_metacal(...)`` keeps working.

    Attributes
    ----------
    resdict : dict
        MetacalBootstrapper result dict (one entry per metacal type).
    reconv : dict
        Averaged metacal *reconvolution*-kernel PSF (round, enlarged):
        :func:`average_multiepoch_psf`.
    orig : dict
        Averaged *original* image-PSF (psfex/mccd, its true shape and size):
        :func:`average_original_psf`.
    """

    resdict: dict
    reconv: dict
    orig: dict


def get_type_flags(fit):
    """Get Type Flags.

    Fit flags of one metacal type, reading absence of evidence of success
    as failure.

    @sc [decision:catalogue_assembly.failure_sentinels,label:convention] mcal-flags-zero-means-measured
    A flag of 0 means the fit ran, reported success and returned a finite
    shear; no default or fallback may produce 0. FLAGS_<SHEAR>, MCAL_FLAGS
    (OR) and MCAL_TYPES_FAIL (count) all derive from this function, and
    sp_validation selects galaxies on MCAL_FLAGS == 0 and
    MCAL_TYPES_FAIL == 0 as "measured". Every failure carries one of
    ngmix's own bits (``ngmix.flags``); ShapePipe adds none:

    - the fitter reported failure: its own ``flags``, unchanged;
    - the fitter reported success (``flags == 0``) without a finite shear
      ``g``, or the result has no ``flags``, or the type is absent:
      ``LM_FUNC_NOTFINITE`` (2**12). The shear the catalogue holds for
      such a type is NaN or a sentinel, never a finite measurement.

    An object ngmix never fit (no usable epoch, or its fit raised) has no
    metacal result at all; make_cat evaluates it as ``{}`` for every type,
    so it carries ``LM_FUNC_NOTFINITE`` in every flag column and fails all
    five types. ``NGMIX_N_EPOCH == 0`` tells it apart from a fitted object
    whose LM fit set the same bit.

    Parameters
    ----------
    fit : dict
        One metacal type's fit result; ``{}`` when the type is absent.

    Returns
    -------
    int
        The fit's own ``flags``, or ``ngmix.flags.LM_FUNC_NOTFINITE`` when
        the result is absent, lacks ``flags``, or claims success
        (``flags == 0``) without a finite shear ``g``.
    """
    notfinite = ngmix.flags.LM_FUNC_NOTFINITE
    flags = int(fit.get('flags', notfinite))
    g = np.asarray(fit.get('g', (np.nan, np.nan)), dtype=float)
    if flags == 0 and not np.all(np.isfinite(g)):
        return notfinite
    return flags


def get_mcal_flags(res):
    """Get Metacal Flags.

    Object-level metacal flags: the bitwise OR of :func:`get_type_flags`
    over :data:`METACAL_TYPES` (the NGMIX_MCAL_FLAGS column).

    Parameters
    ----------
    res : dict
        MetacalBootstrapper result dict with one entry per metacal type.

    Returns
    -------
    int
        OR of all per-type flags; 0 only if every type was measured.
    """
    return int(np.bitwise_or.reduce(
        [get_type_flags(res.get(name, {})) for name in METACAL_TYPES]
    ))


def get_mcal_types_fail(res):
    """Get Metacal Types Fail.

    Number of metacal types (0-5) with nonzero :func:`get_type_flags` (the
    NGMIX_MCAL_TYPES_FAIL column).

    Parameters
    ----------
    res : dict
        MetacalBootstrapper result dict with one entry per metacal type.

    Returns
    -------
    int
        Count of failed metacal types.
    """
    return sum(
        get_type_flags(res.get(name, {})) != 0 for name in METACAL_TYPES
    )


def log_run_health(w_log, count, n_fitted, n_flagged):
    """Log Run Health.

    Log an error when a run's metacal fits failed wholesale: either no
    object fitted at all, or every fitted object carries nonzero
    ``mcal_flags``.

    @sc [label:operations] run-health-logs-not-raises
    Wholesale metacal failure is logged at error level, never raised: one
    empty edge tile or broken input must not abort a multi-tile campaign
    job, and the error line is the signal to catch in review.

    Parameters
    ----------
    w_log : logging.Logger
        Logging instance
    count : int
        Number of objects considered for fitting
    n_fitted : int
        Number of objects that were fitted (present in the results list)
    n_flagged : int
        Number of fitted objects whose ``mcal_flags`` ended up nonzero

    """
    if count > 0 and n_fitted == 0:
        w_log.error(
            f'ngmix: all {count} objects failed the metacal fit'
            ' (0 fitted); writing an empty catalogue. Expected only for a'
            ' tile with no usable epochs; otherwise check the vignettes,'
            ' PSFs and ngmix installation.'
        )
    if n_fitted > 0 and n_flagged == n_fitted:
        w_log.error(
            f'ngmix: 100% of {n_fitted} fitted objects carry nonzero'
            ' mcal_flags -- the metacal fit failed wholesale; outputs'
            ' are unusable.'
        )


def empty_metacal_output():
    """Empty Metacal Output.

    The five-HDU-shaped output dict :meth:`Ngmix.compile_results` returns
    for zero fitted objects: one empty list per column, per metacal type.
    Factored out so a tile with no measurable objects gets the identical
    catalogue shape whether it is discovered by :meth:`Ngmix.process`
    (partly-empty tile, one skipped object at a time) or by one of
    ``ngmix_runner``'s all-empty-store guards (wholesale-empty tile, before
    any object is read).

    Returns
    -------
    dict
        ``{metacal_type: {column: []}}``, matching an empty
        :meth:`Ngmix.compile_results` call.

    Raises
    ------
    ValueError
        If the hardcoded HDU name list drifts out of sync with
        :data:`METACAL_TYPES`.
    """
    # Output HDU order. Same set as METACAL_TYPES, but kept in this
    # fixed order so output catalogues stay byte-reproducible; the check
    # below guards against the two lists silently diverging.
    names = ["1m", "1p", "2m", "2p", "noshear"]
    if set(names) != set(METACAL_TYPES):
        raise ValueError(
            "compile_results metacal type list is out of sync with"
            + " METACAL_TYPES"
        )
    names2 = [
        'id',
        'n_epoch_model',
        'mcal_types_fail',
        'neighbour_flag',
        'nfev_fit',
        # galaxy
        'g1',
        'g1_err',
        'g2',
        'g2_err',
        'T',
        'T_err',
        'flux',
        'flux_err',
        's2n',
        'mag',
        'mag_err',
        'flags',
        'mcal_flags',
        # original image PSF (psfex/mccd), fit by average_original_psf
        'g1_psf_orig',
        'g2_psf_orig',
        'g1_err_psf_orig',
        'g2_err_psf_orig',
        'T_psf_orig',
        'T_err_psf_orig',
        # metacal reconvolution kernel, fit by average_multiepoch_psf
        'g1_psf_reconv',
        'g2_psf_reconv',
        'g1_err_psf_reconv',
        'g2_err_psf_reconv',
        'T_psf_reconv',
        'T_err_psf_reconv',
    ]
    return {k: {kk: [] for kk in names2} for k in names}


def write_ngmix_fits(output_path, output_dict):
    """Write Ngmix Fits.

    Write a compiled ngmix results dict to a fresh output FITS file, one
    HDU per metacal type. The file must not already exist; an existing
    output is appended to by :meth:`Ngmix.save_results`, not this function.

    Parameters
    ----------
    output_path : str
        Path of the FITS file to create
    output_dict : dict
        Compiled results, as returned by :meth:`Ngmix.compile_results` or
        :func:`empty_metacal_output`

    Raises
    ------
    IndexError
        If ``output_dict`` does not have exactly five HDUs
    """
    n_hdu = len(output_dict.keys())
    if n_hdu != 5:
        raise IndexError(
            f"FITS output file data has {n_hdu} HDUs,"
            + " expected are 5"
        )
    f_out = file_io.FITSCatalogue(
        output_path, open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
    )
    for key in output_dict.keys():
        f_out.save_as_fits(output_dict[key], ext_name=key.upper())


def write_empty_tile_output(output_dir, file_number_string, w_log, count):
    """Write Empty Tile Output.

    Write ngmix's empty-tile product for a guard that fires before any
    stamp is read: the same five-HDU empty catalogue and run-health error
    line that :meth:`Ngmix.process` writes once every object in a tile has
    been skipped, without constructing an ``Ngmix`` instance (which would
    open the galaxy vignette store) or reading any vignette.

    Used by ``ngmix_runner``'s all-empty-store guards, which must return
    before the galaxy vignette store is opened at all -- reading its
    many-epoch, many-object arrays is what triggers a C-level malloc crash
    when the store is large (see the PSF-empty guard's own comment).

    Parameters
    ----------
    output_dir : str
        Output directory
    file_number_string : str
        File numbering scheme
    w_log : logging.Logger
        Logging instance
    count : int
        Number of objects considered for fitting (see :func:`log_run_health`);
        for these guards, the number of entries in the empty store.
    """
    log_run_health(w_log, count, n_fitted=0, n_flagged=0)
    output_path = f"{output_dir}/ngmix{file_number_string}.fits"
    if not os.path.exists(output_path):
        write_ngmix_fits(output_path, empty_metacal_output())


def check_wcs_centroid_offset(centroid_source, tile_cat, gal_vign_cat):
    """Check WCS Centroid Offset.

    Fail once, up front, when ``centroid_source="wcs"`` would have no
    coadd-centroid offset to place the galaxy Jacobian at.

    @sc [label:coupling] wcs-centroid-needs-offset
    ``centroid_source="wcs"`` reads the ``OFFSET`` the stamp extractor
    (:func:`shapepipe.modules.vignetmaker_package.vignetmaker.get_stamps`)
    writes into every vignette epoch entry; vignettes cut before that
    extractor carry none. Left unchecked, :func:`make_ngmix_observation`
    raises for every object in turn and :meth:`Ngmix.process`'s per-object
    exception handling turns the whole tile into a silently empty
    catalogue. OFFSET is a property of the extraction run, not of any one
    object, so the first object with epochs speaks for the whole vignette
    file: checking it is enough, and scanning every object would only cost
    more sqlitedict unpickling for the same answer.

    Parameters
    ----------
    centroid_source : {"wcs", "hsm"}
        The configured centroid source; a no-op unless it is ``"wcs"``.
    tile_cat : Tile_cat
        Tile catalogue, read for its object ID order.
    gal_vign_cat : Mapping
        Galaxy vignette store, keyed by ``str(obj_id)``.

    Raises
    ------
    ValueError
        If ``centroid_source == "wcs"`` and the first object with epochs
        has an epoch entry with no ``OFFSET``.
    """
    if centroid_source != "wcs":
        return
    for obj_id in tile_cat.obj_id:
        gal_obj = gal_vign_cat[str(obj_id)]
        if gal_obj == 'empty' or not gal_obj:
            continue
        first_epoch = next(iter(gal_obj.values()))
        if 'OFFSET' not in first_epoch:
            raise ValueError(
                "centroid_source='wcs' requires the coadd-centroid OFFSET"
                " the stamp extractor writes into every vignette epoch,"
                " but this tile's vignettes carry none: re-extract the"
                " stamps with the current vignetmaker, or set"
                " centroid_source='hsm'."
            )
        return


def get_prior(pixel_scale, rng, T_range=None, F_range=None):
    """Build ngmix joint prior for a 6-parameter galaxy model.

    Parameters
    ----------
    pixel_scale : float
        Pixel scale in arcsec (sets centroid prior width).
    rng : numpy.random.RandomState
        Random state for all priors.
    T_range : tuple of float, optional
        (min, max) for flat size prior; default (-1, 1e3).
    F_range : tuple of float, optional
        (min, max) for flat flux prior; default (-100, 1e9).

    Returns
    -------
    ngmix.joint_prior.PriorSimpleSep

    @sc [decision:shape_measurement.fit_priors]
    """
    if T_range is None:
        T_range = [-1.0, 1.0e3]
    if F_range is None:
        F_range = [-100.0, 1.0e9]

    cen_prior = ngmix.priors.CenPrior(
        cen1=0.0, cen2=0.0,
        sigma1=pixel_scale, sigma2=pixel_scale,
        rng=rng,
    )
    g_prior = ngmix.priors.GPriorBA(sigma=0.4, rng=rng)
    T_prior = ngmix.priors.FlatPrior(minval=T_range[0], maxval=T_range[1], rng=rng)
    F_prior = ngmix.priors.FlatPrior(minval=F_range[0], maxval=F_range[1], rng=rng)

    return ngmix.joint_prior.PriorSimpleSep(
        cen_prior=cen_prior,
        g_prior=g_prior,
        T_prior=T_prior,
        F_prior=F_prior,
    )


def chunk_rows(n_obj, row_min, row_max):
    """Catalogue rows of one ngmix chunk.

    A chunk is a closed range of 1-based row positions in the tile
    catalogue, independent of the ``NUMBER`` values those rows carry, so a
    partition of ``1..n_obj`` covers every object once however ``NUMBER``
    is ordered or spaced. A bound ``<= 0`` is unbounded on that side.

    Parameters
    ----------
    n_obj : int
        Number of rows in the tile catalogue
    row_min, row_max : int
        First and last row of the chunk (1-based, inclusive)

    Returns
    -------
    range
        0-based row indices of the chunk

    """
    start = row_min - 1 if row_min > 0 else 0
    stop = min(row_max, n_obj) if row_max > 0 else n_obj
    return range(start, max(start, stop))


def position_seed(ra, dec, ccd):
    """Deterministic RNG seed from an object's sky position (ngmix#796).

    Position seeding gives the same object the same RNG stream in each image
    branch, provided its sky position falls in the same seed box. It also makes
    the result independent of how the tile is split into
    ``ID_OBJ_MIN``/``ID_OBJ_MAX`` row chunks, which is why it is now the only
    mode.

    Box math (kept exactly as Fabian's issue #796)::

        box_x = floor(ra  * 3600 / 3) + (ccd + 1)
        box_y = floor(dec * 3600 / 3) + (ccd + 2)

    The 3-arcsec boxes make the seed robust to detection-order changes and to
    small centroid jitter within a box: the same physical object lands in the
    same box across branches and so gets the same seed. The ``+ccd`` offsets
    disambiguate exposure overlap.

    **Deviation from #796.** Fabian combines the boxes as ``box_x + box_y``,
    which collides along anti-diagonals — e.g. boxes ``(10, 20)`` and
    ``(11, 19)`` both give seed 30, so two distinct objects share a noise
    stream. We instead combine them with a Cantor pairing (a bijection on
    non-negative integers), after a zig-zag fold that maps signed box indices
    (dec can be negative) onto non-negative ones. The result is reduced mod
    ``2**32`` to land in the valid ``numpy.random.RandomState`` seed range
    ``[0, 2**32)``. The box math is untouched; only the collision-prone sum is
    replaced.

    Parameters
    ----------
    ra, dec : float
        Object sky position in degrees (first-epoch coordinates; all epochs of
        one object share them).
    ccd : int
        CCD number of the first epoch.

    Returns
    -------
    int
        Seed in ``[0, 2**32)`` for ``numpy.random.RandomState``.

    @sc [decision:shape_measurement.ngmix_seed_mode]
    """
    box_x = int(np.floor((ra * 3600) / 3) + (ccd + 1))
    box_y = int(np.floor((dec * 3600) / 3) + (ccd + 2))
    # Zig-zag fold signed -> non-negative (0,-1,1,-2,... -> 0,1,2,3,...) so the
    # Cantor pairing, which is a bijection on the non-negative integers, stays
    # injective for southern (dec<0) boxes.
    zx = 2 * box_x if box_x >= 0 else -2 * box_x - 1
    zy = 2 * box_y if box_y >= 0 else -2 * box_y - 1
    cantor = (zx + zy) * (zx + zy + 1) // 2 + zy
    return cantor % (2 ** 32)


class Tile_cat():
    """Tile_cat.

    catalog measured on a tile

    Parameters
    ----------
    cat_path : str
        Path to the tile SExtractor catalogue. Its optional ``SEG_VIGNET``
        column, one integer coadd segmentation stamp per object on the grid
        of its ``VIGNET``, becomes ``self.seg`` for the ``"uberseg"`` blend
        handling; without it ``self.seg`` is ``None``.

    """
    def __init__(
        self,
        cat_path,
    ):
        self.cat_path = cat_path
        if cat_path:
            self.get_data(cat_path)

    def get_data(self, cat_path):
        tile_cat = file_io.FITSCatalogue(
            cat_path,
            SEx_catalogue=True,
        )
        tile_cat.open()
        data = tile_cat.get_data()
        cols = data.dtype.names

        self.obj_id = np.copy(data['NUMBER'])
        self.ra = np.copy(data['XWIN_WORLD'])
        self.dec = np.copy(data['YWIN_WORLD'])

        # Optional columns — may be absent in external (non-SExtractor) catalogs.
        # The stamp columns are the bulk of the table and are views into it,
        # not copies, so it is held once (prepare_postage_stamps copies each
        # object's stamp before changing it).
        self.flux = np.copy(data['FLUX_AUTO']) if 'FLUX_AUTO' in cols else None
        self.vign = data['VIGNET'] if 'VIGNET' in cols else None

        # Coadd-frame segmentation stamp (integer labels, the catalogue's
        # NUMBER), one per object on the grid of its VIGNET, overlaid
        # unchanged on every epoch for uberseg neighbour masking
        # (shapepipe#776).
        self.seg = data['SEG_VIGNET'] if 'SEG_VIGNET' in cols else None

        tile_cat.close()

class Postage_stamp():
    """Galaxy Postage Stamp.

    Class to hold catalog of postage stamps for a single galaxy

    Parameters
    ----------
    bkg_sub: bool

    megacam_flip: bool
    We probably want to put weight and flag options here too

    """
    def __init__(
        self,
        bkg_sub=True,
        megacam_flip=True

    ):
        self.gals = []
        self.psfs = []
        self.weights = []
        self.flags = []
        # Neighbour masks, one per epoch: the pixels marked -1e30 in the tile
        # VIGNET on other detections' footprints (off-tile markers get zero
        # weight instead; see split_tile_markers), MegaCam-flipped to the
        # epoch. noisefill zero-weights and noise-fills them, and the epoch
        # cuts do not count them (see prepare_ngmix_weights).
        self.neighbours = []
        self.bkg_rms = []
        # Segmentation stamps, one per epoch, used only by the "uberseg" blend
        # handling; empty under the default "noisefill". All epochs carry the
        # SAME coadd-frame seg stamp (shapepipe#776: one coadd seg per object,
        # no per-epoch reprojection), each MegaCam-flipped to match its galaxy
        # stamp so the overlay stays registered.
        self.segs = []
        self.jacobs = []
        # Per-epoch sub-pixel [row, col] coadd-centroid offset propagated from
        # the stamp extractor; the "wcs" centroid source places the Jacobian
        # origin there (unused by "hsm").
        self.offsets = []
        # The object's sky position, per epoch; seeds the per-object RNG (see
        # :func:`position_seed`).
        self.ra = []
        self.dec = []
        # CCD number of the first epoch, used only to build the per-object
        # position seed (see :func:`position_seed`).
        self.ccd = None
        self.epoch_cuts = Counter(
            considered=0, masked_fraction=0, central_veto=0
        )
        self.bkg_sub = bkg_sub
        self.megacam_flip = megacam_flip

class Vignet():
    """Vignet.

    Class to hold catalog of postage stamps

    Parameters
    ----------
    gal_vignet_path
    bkg_vignet_path
    psf_vignet_path
    weight_vignet_path
    flag_vignet_path
    f_wcs_path
    bkg_rms_vignet_path
    """
    def __init__(
        self,
        gal_vignet_path,
        bkg_vignet_path,
        psf_vignet_path,
        weight_vignet_path,
        flag_vignet_path,
        f_wcs_path,
        bkg_rms_vignet_path=None,
    ):
        self.f_wcs_file = SqliteDict(f_wcs_path)
        self.gal_vign_cat = SqliteDict(gal_vignet_path)
        self.bkg_vign_cat = SqliteDict(bkg_vignet_path) if bkg_vignet_path is not None else None
        self.psf_vign_cat = SqliteDict(psf_vignet_path)
        self.weight_vign_cat = SqliteDict(weight_vignet_path)
        self.flag_vign_cat = SqliteDict(flag_vignet_path)
        self.bkg_rms_vign_cat = (
            SqliteDict(bkg_rms_vignet_path)
            if bkg_rms_vignet_path is not None
            else None
        )

    def close(self):
        self.f_wcs_file.close()
        self.gal_vign_cat.close()
        if self.bkg_vign_cat is not None:
            self.bkg_vign_cat.close()
        self.flag_vign_cat.close()
        self.weight_vign_cat.close()
        self.psf_vign_cat.close()
        if self.bkg_rms_vign_cat is not None:
            self.bkg_rms_vign_cat.close()

class Ngmix(object):
    """Ngmix.

    Class to handle NGMIX shapepe measurement.

    Parameters
    ----------
    input_file_list : list
        Input files
    output_dir : str
        Output directory
    file_number_string : str
        File numbering scheme
    zero_point : float
        Photometric zero point
    f_wcs_path : str
        Path to merged single-exposure single-HDU headers
    w_log : logging.Logger
        Logging instance
    save_batch : int, optional
        Save output catalogue in batches of this size; detaul is ``-1`` (no
        batch save)
    id_obj_min : int, optional
        First catalogue row to process (1-based, see :func:`chunk_rows`),
        not used if the value is set to ``-1``; the default is ``-1``
    id_obj_max : int, optional
        Last catalogue row to process (1-based, inclusive), not used if the
        value is set to ``-1``; the default is ``-1``
    centroid_source : {"wcs", "hsm"}, optional
        How to place the galaxy Jacobian origin for the centroid prior. The
        default ``"wcs"`` places it at the coadd centroid: the sub-pixel
        offset the stamp extractor computed when it cut the stamp,
        propagated on the vignette. ``"hsm"`` re-centers on the
        adaptive-moment centroid measured from the stamp pixels. See
        :func:`make_ngmix_observation`.
    blend_handling : {"noisefill", "uberseg"}, optional
        Neighbour treatment. ``"noisefill"`` (default) zero-weights and
        noise-fills the pixels marked -1e30 in the tile VIGNET on other
        detections' footprints; ``"uberseg"`` ignores those markers and
        zeroes the weight of neighbour-side pixels from the coadd
        segmentation stamps, the tile catalogue's ``SEG_VIGNET`` column (see
        :class:`Tile_cat`), which it requires. Defect pixels are filled
        under both (see :func:`prepare_ngmix_weights`).
    dilate_neighbour : int, optional
        Neighbour-mask dilation iterations for ``"uberseg"`` (see
        :func:`uberseg_mask`); the default is ``1``.

    Notes
    -----
    The RNG is always per object and seeded from that object's sky position;
    :func:`position_seed` says what that buys.

    Raises
    ------
    IndexError
        If the length of the input file list is incorrect
    ValueError
        If ``blend_handling`` is unknown.

    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        file_number_string,
        zero_point,
        f_wcs_path,
        w_log,
        save_batch=-1,
        id_obj_min=-1,
        id_obj_max=-1,
        bkg_sub=True,
        centroid_source="wcs",
        blend_handling="noisefill",
        dilate_neighbour=1,
        metacal_psf="fitgauss",
    ):

        # Base count = catalogue + vignets, excluding the f_wcs headers (passed
        # separately). One fewer when the background vignet is absent
        # (``bkg_sub=False``, image sims). An extra trailing slot may carry the
        # optional background-rms vignet (#779), so both counts are valid.
        n_base = 6 if bkg_sub else 5
        if len(input_file_list) not in {n_base, n_base + 1}:
            raise IndexError(
                f"Input file list has length {len(input_file_list)},"
                + f" required is {n_base} or {n_base + 1}"
            )

        if blend_handling not in BLEND_HANDLINGS:
            raise ValueError(
                f"Unknown BLEND_HANDLING '{blend_handling}'; expected one of"
                + f" {BLEND_HANDLINGS}"
            )

        self._tile_cat_path = input_file_list[0]
        if bkg_sub:
            bkg_path, psf_path, weight_path, flag_path = (
                input_file_list[2], input_file_list[3],
                input_file_list[4], input_file_list[5],
            )
        else:
            bkg_path, psf_path, weight_path, flag_path = (
                None, input_file_list[2],
                input_file_list[3], input_file_list[4],
            )
        bkg_rms_vignet_path = (
            input_file_list[n_base]
            if len(input_file_list) == n_base + 1
            else None
        )
        self._vignet_cat = Vignet(
            input_file_list[1],
            bkg_path,
            psf_path,
            weight_path,
            flag_path,
            f_wcs_path,
            bkg_rms_vignet_path,
        )

        self._output_dir = output_dir
        self._file_number_string = file_number_string

        self._zero_point = zero_point

        self._f_wcs_path = f_wcs_path

        self._save_batch = save_batch
        self._id_obj_min = id_obj_min
        self._id_obj_max = id_obj_max
        self._bkg_sub = bkg_sub
        self._centroid_source = centroid_source
        self._blend_handling = blend_handling
        self._dilate_neighbour = dilate_neighbour
        self._metacal_psf = metacal_psf

        self._w_log = w_log

        self._w_log.info(
            'Per-object RNG seeded from sky position (ngmix#796): results are'
            ' invariant to how the tile is split into object chunks'
        )

    @classmethod
    def MegaCamFlip(self, vign, ccd_nb):
        """Flip for MegaCam.

        MegaPipe has CCDs that are upside down. This function flips the
        postage stamps in these CCDs. TO DO: This will give incorrect results
        when used with THELI ccds.  Fix this.

        Parameters
        ----------
        vign : numpy.ndarray
            Array containing the postage stamp to flip
        ccd_nb : int
            ID of the CCD containing the postage stamp

        Returns
        -------
        numpy.ndarray
            The flipped postage stamp

        @sc [decision:shape_measurement.megacam_ccd_flip]
        """
        if ccd_nb < 18 or ccd_nb in [36, 37]:
            # swap x axis so origin is on top-right
            return np.rot90(vign, k=2)
        else:
            # swap y axis so origin is on bottom-left
            return vign

    def compile_results(self, results):
        """Compile Results.

        Prepare the results of NGMIX before saving. TO DO: add snr_r and T_r
        This needs to be updated
        Parameters
        ----------
        results : dict
            Results of NGMIX metacal

        Returns
        -------
        dict
            Compiled results ready to be written to a file.

            Two PSF column families — each carrying ellipticity *and* size,
            for *different* PSFs (shapepipe#749). See the function docstrings
            for what each family IS:

            * ``*_psf_orig`` (``g1``/``g2`` + ``*_err``, ``T``) — the original
              image PSF, fit by :func:`average_original_psf`.
            * ``*_psf_reconv`` — the metacal reconvolution kernel, fit by
              :func:`average_multiepoch_psf`.

        Raises
        ------
        KeyError
            If SNR key not found

        """
        # Column layout (HDU names and per-type columns) lives in
        # empty_metacal_output, shared with the runner's all-empty-store
        # guards so every zero-object catalogue has the identical shape.
        output_dict = empty_metacal_output()
        names = list(output_dict.keys())
        for idx in range(len(results)):
            # Object-level quality columns, derived from the same per-type
            # flags as the ``flags`` column below (see get_type_flags).
            mcal_flags = get_mcal_flags(results[idx])
            mcal_types_fail = get_mcal_types_fail(results[idx])
            for name in names:
                fit = results[idx].get(name, {})
                flags = get_type_flags(fit)

                # ngmix 2.x does not raise on fit failure: after ntry the
                # result keeps flags != 0 and carries none of the
                # measurement keys (g, g_cov, T, T_err, flux, flux_err,
                # s2n). NaN-fill those (and an absent type) so failed types
                # are recorded with their flags instead of crashing the tile
                # on a KeyError.
                flux = fit.get("flux", np.nan)
                flux_err = fit.get("flux_err", np.nan)
                g = np.asarray(fit.get("g", (np.nan, np.nan)))
                g_cov = np.asarray(
                    fit.get("g_cov", np.full((2, 2), np.nan))
                )
                T_gal = fit.get("T", np.nan)
                T_gal_err = fit.get("T_err", np.nan)

                mag = -2.5 * np.log10(flux) + self._zero_point
                mag_err = np.abs(-2.5 * flux_err / (flux * np.log(10)))

                output_dict[name]["id"].append(results[idx]["obj_id"])
                output_dict[name]["n_epoch_model"].append(
                    results[idx]["n_epoch_model"]
                )
                output_dict[name]["mcal_types_fail"].append(mcal_types_fail)
                # Per-object blend flag (see process()); replicated across all
                # shear types like id / n_epoch_model / mcal_types_fail.
                output_dict[name]["neighbour_flag"].append(
                    results[idx]["neighbour_flag"]
                )
                # ngmix 2.x reports the solver's function-evaluation count
                # (nfev, ~tens-hundreds; -1 on some failures), not the v1
                # 1-5 retry count, so the column is named accordingly. Fits
                # that fail before nfev is ever set (e.g. the NaN-filled
                # branch above) fall back to the same -1 sentinel, keeping
                # the column int64 (FITS 'K') across every batch: nfev_fit
                # is int64 with -1 meaning failed/absent.
                output_dict[name]["nfev_fit"].append(
                    fit.get("nfev", -1)
                )
                # The two PSF families are object-level (one value per
                # object, not per shear type) and self-named: every key
                # below is copied straight through from compile-loop input to
                # output, so the column name *is* the value's provenance.
                #   *_psf_orig   = original image PSF (average_original_psf)
                #   *_psf_reconv = reconvolution kernel (average_multiepoch_psf)
                for psf_key in (
                    'g1_psf_orig', 'g2_psf_orig',
                    'g1_err_psf_orig', 'g2_err_psf_orig',
                    'T_psf_orig', 'T_err_psf_orig',
                    'g1_psf_reconv', 'g2_psf_reconv',
                    'g1_err_psf_reconv', 'g2_err_psf_reconv',
                    'T_psf_reconv', 'T_err_psf_reconv',
                ):
                    output_dict[name][psf_key].append(results[idx][psf_key])

                output_dict[name]["g1"].append(g[0])
                output_dict[name]["g2"].append(g[1])
                output_dict[name]["g1_err"].append(np.sqrt(g_cov[0, 0]))
                output_dict[name]["g2_err"].append(np.sqrt(g_cov[1, 1]))
                output_dict[name]["T"].append(T_gal)
                output_dict[name]["T_err"].append(T_gal_err)
                output_dict[name]["flux"].append(flux)
                output_dict[name]["flux_err"].append(flux_err)
                output_dict[name]["mag"].append(mag)
                output_dict[name]["mag_err"].append(mag_err)

                if "s2n" in fit:
                    output_dict[name]["s2n"].append(fit["s2n"])
                elif "s2n_r" in fit:
                    output_dict[name]["s2n"].append(fit["s2n_r"])
                elif flags != 0:
                    output_dict[name]["s2n"].append(np.nan)
                else:
                    raise KeyError("No SNR key (s2n, s2n_r) found in results")

                output_dict[name]["flags"].append(flags)
                output_dict[name]["mcal_flags"].append(mcal_flags)

        return output_dict

    def get_output_path(self, directory):
        """Get Output Path.

        Return path of output ngmix catalogue file.

        Parameters
        ----------
        directoy: str
            directory name

        Returns
        -------
        str
            output path

        """
        return f"{directory}/ngmix{self._file_number_string}.fits"


    def save_results(self, output_dict):
        """Save Results.

        Save the results into a FITS file.

        Parameters
        ----------
        output_dict: dict
            Dictionary containing the results

        """
        n_hdu = len(output_dict.keys())
        if n_hdu != 5:
            raise IndexError(
                f"FITS output file data has {n_hdu} HDUs,"
                + " expected are 5"
            )

        output_name = self.get_output_path(self._output_dir)
        if not os.path.exists(output_name):
            write_ngmix_fits(output_name, output_dict)
            return

        with fits.open(output_name, mode='update') as hdul:

            # Iterate through HDUs (assuming they are Binary Table HDUs with data)
            for idx, hdu in enumerate(hdul):
                if isinstance(hdu, fits.BinTableHDU):  # Check for table data

                    # HDU extension name
                    ext_name = hdu.name.lower()
                    if ext_name not in output_dict:
                        raise ValueError(
                            f"HDU extension {ext_name} from existing FITS"
                            + " file not found in data"
                        )

                    # Existing data
                    existing_data = hdu.data
                    existing_dtype = hdu.data.dtype

                    # New data
                    new_data = output_dict[ext_name]

                    # Verify that all column names in existing_data exist
                    # in new_data_dict
                    if not all(
                        colname in new_data
                        for colname in existing_dtype.names
                    ):
                        raise ValueError(
                            "Mismatch between existing columns"
                            + f" ({existing_dtype.names}) and new data"
                            + f" columns ({list(new_data)})."
                        )

                    # New data to be appended
                    structured_data = np.zeros(
                        len(next(iter(new_data.values()))),
                        dtype=existing_data.dtype,
                    )
                    for colname in existing_data.dtype.names:
                        structured_data[colname] = new_data[colname]

                    # Combine existing and new data
                    updated_data = np.append(existing_data, structured_data)

                    # Update the data in the HDU
                    hdu.data = updated_data

            # Save changes to the FITS file
            hdul.flush()

    @classmethod
    def check_key(self, expccd_name_tmp, vign_cat, vignet_path):
        if expccd_name_tmp not in vign_cat:
            raise KeyError(
                f"Key '{expccd_name_tmp}' (exposure CCD ID from PSF postage stamp list)"
                + " not found in postage stamp database"
                + f" file '{vignet_path}'"
            )

    def _check_central_seg_label(self, seg, obj_id):
        """Cross-check the seg stamp's centre pixel against the object's NUMBER.

        The authoritative central label is ``obj_id`` (the SExtractor NUMBER);
        this is a diagnostic on the coadd-vs-epoch overlay registration, not a
        gate. Two outcomes:

        * The seg stamp does not contain ``obj_id`` anywhere — the object's own
          footprint fell outside the cutout (bad centroid / dropped detection),
          so uberseg has no central object to keep. Raise; the per-object
          ``try/except`` in :meth:`process` drops it, loud in the log.
        * The centre-pixel label is nonzero and ``!= obj_id`` — a few-pixel
          registration offset put a neighbour on the centre pixel. Harmless for
          the mask (which keys off ``obj_id``), but a signal worth surfacing:
          log a warning and proceed.

        Parameters
        ----------
        seg : numpy.ndarray
            Coadd segmentation stamp for this object (integer labels).
        obj_id : int
            The object's SExtractor NUMBER — the authoritative central label.
        """
        if not np.any(seg == obj_id):
            raise ValueError(
                f"seg stamp for object NUMBER={obj_id} contains that label"
                + " nowhere; the object's footprint fell outside the cutout"
                + " (bad centroid / registration). Cannot define the central"
                + " object for uberseg."
            )
        cy, cx = seg.shape[0] // 2, seg.shape[1] // 2
        centre_label = int(seg[cy, cx])
        if centre_label != 0 and centre_label != obj_id:
            self._w_log.info(
                f"uberseg: object NUMBER={obj_id} but seg centre pixel carries"
                + f" label {centre_label}; a coadd-vs-epoch registration offset"
                + " is placing a neighbour on the centre. Proceeding with the"
                + " NUMBER (authoritative)."
            )

    def log_mean_ellipticity(self):
        """Log mean ellipticity from NOSHEAR HDU to the run log.

        Reports <e1>, <e2> with standard errors for all objects and for
        objects passing the default metacal cuts (flags==0, mcal_flags==0,
        10 < SNR < 500, T/Tpsf > 0.5).
        """
        output_path = self.get_output_path(self._output_dir)
        try:
            with fits.open(output_path) as hdul:
                d = hdul['NOSHEAR'].data
                g1 = d['g1'].astype(float)
                g2 = d['g2'].astype(float)
                flags = d['flags']
                mcal_flags = d['mcal_flags']
                s2n = d['s2n'].astype(float)
                T = d['T'].astype(float)
                Tpsf = d['Tpsf'].astype(float)
        except Exception as e:
            self._w_log.warning(f"Could not compute mean ellipticity: {e}")
            return

        n_total = len(g1)
        if n_total == 0:
            self._w_log.info("Mean ellipticity: no objects in output catalogue")
            return

        def _log_stats(g1_sel, g2_sel, label):
            n = len(g1_sel)
            if n == 0:
                self._w_log.info(f"  {label}: 0 objects")
                return
            mean_g1 = g1_sel.mean()
            mean_g2 = g2_sel.mean()
            err_g1 = g1_sel.std() / np.sqrt(n)
            err_g2 = g2_sel.std() / np.sqrt(n)
            self._w_log.info(
                f"  {label} (N={n}):"
                f"  <e1> = {mean_g1:+.4e} +/- {err_g1:.4e},"
                f"  <e2> = {mean_g2:+.4e} +/- {err_g2:.4e}"
            )

        self._w_log.info(f"Mean ellipticity (NOSHEAR, N_total={n_total}):")
        _log_stats(g1, g2, "no cuts")

        with np.errstate(invalid='ignore'):
            mask = (
                (flags == 0)
                & (mcal_flags == 0)
                & (s2n >= 10.0)
                & (s2n <= 500.0)
                & (T / Tpsf >= 0.5)
            )
        _log_stats(g1[mask], g2[mask], "SNR in [10, 500], T/Tpsf > 0.5")

    def process(self):
        """Process.

        Funcion to processs NGMIX.
        organizes object cutouts from detection catalog in image, 
        weight, and flag files
        per object: 
            gathers wcs and psf info from exposures
            background subtracts (make this an option)
            scales by relative zeropoints
            runs metacal convolutions and ngmix fitting
        Returns
        -------
        dict
            Dictionary containing the NGMIX metacal results

        Raises
        ------
        ValueError
            If ``centroid_source == "wcs"`` and the vignette catalogue
            carries no coadd-centroid OFFSET (see
            :func:`check_wcs_centroid_offset`).

        @sc [decision:shape_measurement.fit_initialisation,decision:shape_measurement.ngmix_seed_mode]
        """
        tile_cat = Tile_cat(self._tile_cat_path)
        vignet_cat = self._vignet_cat

        # Fail before the per-object loop, whose try/except would otherwise
        # drop every object one by one (shapepipe#776).
        if self._blend_handling == "uberseg" and tile_cat.seg is None:
            raise ValueError(
                "BLEND_HANDLING = uberseg needs the tile catalogue's"
                + f" SEG_VIGNET column, which {self._tile_cat_path} lacks;"
                + " write it at tile detection (SEG_VIGNET = True)."
            )

        check_wcs_centroid_offset(
            self._centroid_source, tile_cat, vignet_cat.gal_vign_cat
        )

        final_res = []

        count = 0
        n_empty_cat = 0
        n_no_epoch = 0
        n_ngmix_fail = 0
        n_fitted = 0
        n_flagged = 0
        id_first = -1
        id_last = -1
        count_batch = 0
        epoch_cuts = Counter(considered=0, masked_fraction=0, central_veto=0)
        n_emptied = 0
        saved_batch_cumul = 0

        rows = chunk_rows(
            len(tile_cat.obj_id), self._id_obj_min, self._id_obj_max
        )
        for i_tile in rows:
            obj_id = tile_cat.obj_id[i_tile]
            if id_first == -1:
                id_first = obj_id
            id_last = obj_id
            count += 1

            # Skip objects with no multi-epoch PSF or vignet data.
            # Read each store once here and pass the dicts down: every
            # sqlitedict access unpickles the object's whole all-epoch dict.
            psf_obj = vignet_cat.psf_vign_cat[str(obj_id)]
            # Avoid allocating galaxy stamp arrays when there is no PSF coverage.
            if psf_obj == 'empty' or not psf_obj:
                n_empty_cat += 1
                continue
            gal_obj = vignet_cat.gal_vign_cat[str(obj_id)]
            if gal_obj == 'empty' or not gal_obj:
                n_empty_cat += 1
                continue

            stamp = prepare_postage_stamps(
                vignet_cat,
                obj_id,
                i_tile,
                tile_cat,
                self._bkg_sub,
                psf_obj,
                gal_obj,
                blend_handling=self._blend_handling,
            )
            epoch_cuts.update(stamp.epoch_cuts)

            if len(stamp.gals) == 0:
                n_no_epoch += 1
                n_emptied += stamp.epoch_cuts["considered"] > 0
                continue

            # Per-object RNG, seeded from (ra, dec, ccd) — see
            # :func:`position_seed`. The prior is rebuilt from that same RNG
            # because the guesser draws its initial guess via prior.sample()
            # (ngmix guessers.py), which consumes the RNG the prior was
            # CONSTRUCTED with: a per-object rng alone would leave the guess
            # drawing from a shared stream and break the invariance. The
            # centroid prior is one pixel wide, in this object's own epochs.
            obj_rng = np.random.RandomState(
                position_seed(stamp.ra[0], stamp.dec[0], stamp.ccd)
            )
            obj_prior = get_prior(stamp_pixel_scale(stamp.jacobs), obj_rng)

            try:
                flux_guess = (
                    tile_cat.flux[i_tile]
                    if tile_cat.flux is not None
                    else 1.0
                )
                # The central object's segmentation label IS its SExtractor
                # NUMBER (obj_id): seg labels are the NUMBERs of the same SE run
                # that produced the tile catalogue, so obj_id is authoritative
                # and is used for both the uberseg mask and the neighbour flag.
                # The centre pixel is only a diagnostic: a few-pixel coadd-vs-
                # epoch overlay offset can put a NEIGHBOUR's label there, which
                # would invert the mask — so we never trust it to define the
                # central object. Cross-check it and warn on a mismatch (a
                # registration signal), but proceed with obj_id (shapepipe#776).
                if self._blend_handling == "uberseg":
                    self._check_central_seg_label(tile_cat.seg[i_tile], obj_id)
                res, psf_res, psf_orig_res = do_ngmix_metacal(
                    stamp,
                    obj_prior,
                    flux_guess,
                    obj_rng,
                    centroid_source=self._centroid_source,
                    blend_handling=self._blend_handling,
                    object_number=obj_id,
                    dilate_neighbour=self._dilate_neighbour,
                    metacal_psf=self._metacal_psf,
                )
            except Exception as ee:
                self._w_log.info(
                    f'ngmix failed for object ID={obj_id}.\nMessage: {ee}'
                )
                n_ngmix_fail += 1
                continue

            res['obj_id'] = obj_id
            # Neighbour flag: does the coadd seg stamp hold any non-central,
            # non-zero label? Systematics-test hook (shapepipe#776); computed
            # from the raw coadd seg, 0 when uberseg is off / seg absent.
            res['neighbour_flag'] = int(
                seg_has_neighbour(tile_cat.seg[i_tile], obj_id)
                if getattr(tile_cat, "seg", None) is not None
                else 0
            )
            # epochs that survived the PSF fit and entered the model,
            # not the number of epochs submitted (v1 contract)
            res['n_epoch_model'] = psf_res['n_epoch']
            # The mcal flag columns are derived from the per-type results in
            # compile_results; here they only feed the run-health count.
            if get_mcal_flags(res) != 0:
                n_flagged += 1
            # Two distinct PSF families (shapepipe#749), each carrying its own
            # ellipticity AND size: the metacal reconvolution kernel (psf_res)
            # and the original image PSF (psf_orig_res). Tag both into res from
            # ONE template so the families stay symmetric — same quantities,
            # parallel ``*_psf_{family}`` names — and cannot drift apart by a
            # hand-edit to one and not the other.
            for family, psf in (("reconv", psf_res), ("orig", psf_orig_res)):
                res[f"g1_psf_{family}"] = psf["g_psf"][0]
                res[f"g2_psf_{family}"] = psf["g_psf"][1]
                res[f"g1_err_psf_{family}"] = psf["g_psf_err"][0]
                res[f"g2_err_psf_{family}"] = psf["g_psf_err"][1]
                res[f"T_psf_{family}"] = psf["T_psf"]
                res[f"T_err_psf_{family}"] = psf["T_psf_err"]
            final_res.append(res)
            n_fitted += 1
            count_batch += 1

            if self._save_batch > 0 and count_batch == self._save_batch:
                res_dict = self.compile_results(final_res)
                self.save_results(res_dict)
                saved_batch_cumul += count_batch
                self._w_log.info(
                    f"Batch-saved {count_batch} ({len(res_dict)} valid) objects,"
                    + f" cumul={saved_batch_cumul}"
                )
                final_res = []
                count_batch = 0

            if count % 1000 == 0:
                self._w_log.info(
                    f"Progress: {count} iterated, {n_empty_cat} empty catalog,"
                    + f" {n_no_epoch} no valid epoch, {n_ngmix_fail} fit failed,"
                    + f" {n_fitted} fitted"
                )

        self._w_log.info(
            f"ngmix loop finished: {count} iterated"
            + f" (id {id_first}/{id_last}),"
            + f" {n_empty_cat} empty catalog,"
            + f" {n_no_epoch} no valid epoch,"
            + f" {n_ngmix_fail} fit failed,"
            + f" {n_fitted} fitted"
        )

        self._w_log.info(
            "epoch cuts:"
            + f" considered={epoch_cuts['considered']}"
            + f" masked_fraction={epoch_cuts['masked_fraction']}"
            + f" central_veto={epoch_cuts['central_veto']}"
            + f" objects_emptied={n_emptied}"
        )
        log_run_health(self._w_log, count, n_fitted, n_flagged)

        vignet_cat.close()

        # Put all results together
        res_dict = self.compile_results(final_res)

        # Save results
        self.save_results(res_dict)

        # Log mean ellipticity statistics
        self.log_mean_ellipticity()

def prepare_postage_stamps(
    vignet,
    obj_id,
    i_tile,
    tile_cat,
    bkg_sub=True,
    psf_obj=None,
    gal_obj=None,
    blend_handling="noisefill",
):
    """Gather one object's epoch stamps, dropping epochs its defects spoil.

    @sc [decision:shape_measurement.central_defect_veto,decision:shape_measurement.epoch_masked_fraction_cut] epoch-cut-on-defect-mask
    An epoch is dropped when more than ``EPOCH_MASKED_FRACTION_CUT`` of its
    stamp lies in :func:`defect_mask`, the set :func:`prepare_ngmix_weights`
    zero-weights and fills, or when :func:`central_defect_vetoes` finds a
    defect near the object. The quarter-turn weight copies of interpolated
    pixels keep their light and are not counted.

    @sc [decision:shape_measurement.central_defect_veto,decision:shape_measurement.epoch_masked_fraction_cut,decision:shape_measurement.blend_handling] neighbour-markers-are-not-defects
    The tile VIGNET's -1e30 markers on other detections' footprints
    (:func:`split_tile_markers`) are the epoch's neighbour mask
    (``stamp.neighbours``), not defects: every epoch shares the tile VIGNET,
    so counting them would drop every epoch of a blended object.

    @sc [decision:shape_measurement.central_defect_veto,decision:shape_measurement.epoch_masked_fraction_cut,decision:shape_measurement.defect_fill] off-tile-pixels-are-defects
    Beyond the tile's edge the epoch holds the object's own light, cut off by
    the tile, so those pixels get zero exposure weight and join the defect
    set under every ``blend_handling``. On a 51-px stamp an object within
    about 8.5 px of the tile edge fails the masked-fraction cut.

    Parameters
    ----------
    vignet : Vignet
        Per-object vignet stores.
    obj_id : int
        Object ID (SExtractor ``NUMBER``).
    i_tile : int
        Row of the object in ``tile_cat``.
    tile_cat : Tile_cat
        Tile catalogue.
    bkg_sub : bool, optional
        Subtract the background vignet; the default is ``True``.
    psf_obj, gal_obj : dict, optional
        The object's PSF and galaxy vignet dicts, if already read.
    blend_handling : {"noisefill", "uberseg"}, optional
        The neighbour treatment :func:`prepare_ngmix_weights` will apply,
        which decides which defects it can interpolate and so their veto
        radius; the default is ``"noisefill"``.

    Returns
    -------
    Postage_stamp
        The surviving epochs' stamps.
    """
    # define per-object lists of individual exposures to go into ngmix
    stamp = Postage_stamp(bkg_sub=bkg_sub)
    # Read each store's per-object dict ONCE: every sqlitedict access
    # unpickles the object's whole all-epoch dict, so keeping these out of
    # the epoch loop below saves O(n_epoch) full unpickles per store.
    if psf_obj is None:
        psf_obj = vignet.psf_vign_cat[str(obj_id)]
    if gal_obj is None:
        gal_obj = vignet.gal_vign_cat[str(obj_id)]
    bkg_obj = (
        vignet.bkg_vign_cat[str(obj_id)]
        if stamp.bkg_sub and vignet.bkg_vign_cat is not None
        else None
    )
    flag_obj = vignet.flag_vign_cat[str(obj_id)]
    weight_obj = vignet.weight_vign_cat[str(obj_id)]
    bkg_rms_obj = (
        vignet.bkg_rms_vign_cat[str(obj_id)]
        if vignet.bkg_rms_vign_cat is not None
        else None
    )
    wcs_cache = {}
    #identify exposure and ccd number from psf catalog
    psf_expccd_names = list(psf_obj.keys())
    for expccd_name in psf_expccd_names:
        exp_name, ccd_n = re.split('-', expccd_name)

        gal_vign = gal_obj[expccd_name]['VIGNET']

        if np.all(gal_vign == 0):
            continue
        
        if stamp.bkg_sub:
            bkg_vign = bkg_obj[expccd_name]['VIGNET']
            gal_vign_sub_bkg = background_subtract(
                gal_vign,
                bkg_vign
            )
        else:
            gal_vign_sub_bkg = gal_vign

        # Skip epochs where sigma_mad=0 (CCD edge: mostly-zero stamp).
        # prepare_ngmix_weights divides by sig_noise; zero sigma causes NaN/inf
        # weights that corrupt GalSim's C-level FFT allocations.
        if not sigma_mad(gal_vign_sub_bkg) > 0:
            continue

        tile_vign = (
            np.copy(tile_cat.vign[i_tile])
            if tile_cat.vign is not None
            else None
        )
        if stamp.megacam_flip and tile_vign is not None:
            tile_vign = Ngmix.MegaCamFlip(tile_vign, int(ccd_n))

        # Coadd-frame segmentation stamp for this object, MegaCam-flipped to
        # match the (flipped) galaxy stamp so the overlay stays registered.
        # The SAME coadd seg rides every epoch (shapepipe#776: no reprojection).
        tile_seg = (
            np.copy(tile_cat.seg[i_tile])
            if getattr(tile_cat, "seg", None) is not None
            else None
        )
        if stamp.megacam_flip and tile_seg is not None:
            tile_seg = Ngmix.MegaCamFlip(tile_seg, int(ccd_n))

        flag_vign = flag_obj[expccd_name]['VIGNET']
        # Off-tile pixels are defects (off-tile-pixels-are-defects); the
        # other -1e30 markers are neighbours (neighbour-markers-are-not-defects).
        neighbour, off_tile = split_tile_markers(tile_vign, np.shape(gal_vign))
        weight_vign = np.where(off_tile, 0, weight_obj[expccd_name]['VIGNET'])
        bkg_rms_vign = (
            bkg_rms_obj[expccd_name]['VIGNET']
            if bkg_rms_obj is not None
            else None
        )
        defect = defect_mask(weight_vign, flag_vign, bkg_rms_vign)
        stamp.epoch_cuts["considered"] += 1
        if defect.mean() > EPOCH_MASKED_FRACTION_CUT:
            stamp.epoch_cuts["masked_fraction"] += 1
            continue
        removed = neighbour if blend_handling == "noisefill" else None
        if central_defect_vetoes(
            defect, interpolable_defects(defect, removed)
        ):
            stamp.epoch_cuts["central_veto"] += 1
            continue

        # One unpickle per exposure (all CCDs), reused across this object's
        # epochs; the cache is per-call, so bounded by the object's exposures
        if exp_name not in wcs_cache:
            wcs_cache[exp_name] = vignet.f_wcs_file[exp_name]
        ccd_wcs = wcs_cache[exp_name][int(ccd_n)]

        epoch_wcs = ccd_wcs['WCS']
        jacob = get_galsim_jacobian(
            epoch_wcs,
            tile_cat.ra[i_tile],
            tile_cat.dec[i_tile]
        )

        header = fits.Header.fromstring(ccd_wcs['header'])

        # rescale by relative zero-points
        (
            gal_vign_scaled,
            weight_vign_scaled,
            bkg_rms_vign_scaled,
        ) = rescale_epoch_fluxes(
            gal_vign_sub_bkg,
            weight_vign,
            header,
            bkg_rms_vign,
        )

        # gather postage stamps in all of the epochs
        stamp.gals.append(gal_vign_scaled)
        stamp.psfs.append(psf_obj[expccd_name]['VIGNET'])
        stamp.weights.append(weight_vign_scaled)
        stamp.flags.append(flag_vign)
        stamp.neighbours.append(neighbour)
        stamp.bkg_rms.append(bkg_rms_vign_scaled)
        if tile_seg is not None:
            stamp.segs.append(tile_seg)
        stamp.jacobs.append(jacob)
        # Coadd-centroid offset the stamp extractor used, propagated on the
        # galaxy vignette. Consumed only by the "wcs" centroid source (see
        # make_ngmix_observation), which raises if it is missing; the "hsm"
        # path ignores it, so read it leniently rather than coupling hsm to a
        # field it never uses.
        stamp.offsets.append(gal_obj[expccd_name].get('OFFSET'))
        stamp.ra.append(tile_cat.ra[i_tile])
        stamp.dec.append(tile_cat.dec[i_tile])
        # CCD of the first surviving epoch — Fabian's coord_list[0] convention
        # for the position seed. All epochs of one object share the ra/dec
        # above, so first-epoch CCD pins one deterministic seed stream.
        if stamp.ccd is None:
            stamp.ccd = int(ccd_n)

    return stamp

def split_tile_markers(tile_vign, shape):
    """Split the tile VIGNET's -1e30 markers into neighbour and off-tile.

    @sc [decision:shape_measurement.blend_handling,decision:shape_measurement.defect_fill] off-tile-is-marked-border-rows-and-columns
    The tile VIGNET holds -1e30 on the footprints of other detections and on
    stamp pixels beyond the tile's edge (SExtractor writes it). A stamp
    clipped by the tile's rectangle loses whole rows and whole columns from
    its border, so the off-tile pixels are the union of the runs of
    entirely -1e30 rows and columns that start at a stamp border. The remaining markers are
    neighbour pixels: a footprint touching the stamp border, and a footprint
    that completes an interior row or column beside an off-tile band, stay
    neighbours.

    Parameters
    ----------
    tile_vign : numpy.ndarray or None
        Tile VIGNET stamp, oriented like the epoch; ``None`` marks nothing.
    shape : tuple of int
        Stamp shape, used when ``tile_vign`` is ``None``.

    Returns
    -------
    numpy.ndarray of bool
        Neighbour pixels.
    numpy.ndarray of bool
        Off-tile pixels.
    """
    if tile_vign is None:
        return np.zeros(shape, dtype=bool), np.zeros(shape, dtype=bool)
    marker = tile_vign == -1e30

    def border_runs(full):
        # Lines in the unbroken run of marked lines from either border.
        lead = np.logical_and.accumulate(full)
        trail = np.logical_and.accumulate(full[::-1])[::-1]
        return lead | trail

    off_tile = (
        border_runs(marker.all(axis=1))[:, None]
        | border_runs(marker.all(axis=0))[None, :]
    )
    return marker & ~off_tile, off_tile


def background_subtract(gal,bkg):
    """background subtraction.
        
    Parameters
    ----------
    gal : numpy.ndarray
        galaxy image
    bkg : numpy.ndarray
        background
        
    Returns
    -------
    numpy.ndarray
        background subtracted galaxy
    @sc [decision:shape_measurement.galaxy_pixel_weights]
    """

    # background subtraction
    gal_vign_sub_bkg = gal - bkg

    return gal_vign_sub_bkg

def rescale_epoch_fluxes(gal, weight, header, bkg_rms=None):
    """rescale epochs by relative zeropoints to be on the same flux scale
        
    Parameters
    ----------
    gal : numpy.ndarray
        background subtracted galaxy image
    weight : numpy.ndarray
        weight image
    header : 
        image header
    bkg_rms : numpy.ndarray, optional
        Background RMS image
        
    Returns
    -------
    numpy.ndarray
        rescaled galaxy image
    numpy.ndarray
        rescaled weight image
    numpy.ndarray or None
        rescaled background RMS image
    @sc [decision:shape_measurement.epoch_flux_rescaling]
    """
    Fscale = header['FSCALE']

    gal_scaled = gal * Fscale
    weight_scaled = weight * 1 / Fscale ** 2
    bkg_rms_scaled = bkg_rms * Fscale if bkg_rms is not None else None

    return gal_scaled, weight_scaled, bkg_rms_scaled

def get_galsim_jacobian(wcs, ra, dec):
    """Get local wcs.
    This produces a galsim jacobian at a point.  We call it local_wcs because we convert to a ngmix object to create the jacobian later.
    TO DO: can we do this within ngmix?

    Parameters
    ----------
    wcs : astropy.wcs.WCS
        WCS object for which we want the Jacobian
    ra : float
        RA position of the center of the vignet (in degrees)
    dec : float
        Dec position of the center of the vignet (in degress)

    Returns
    -------
    galsim.wcs.BaseWCS.jacobian
        Jacobian of the WCS at the required position

    """
    g_wcs = galsim.fitswcs.AstropyWCS(wcs=wcs)
    world_pos = galsim.CelestialCoord(
        ra=ra * galsim.angle.degrees,
        dec=dec * galsim.angle.degrees,
    )
    galsim_jacob = g_wcs.jacobian(world_pos=world_pos)

    return galsim_jacob


def stamp_pixel_scale(jacobs):
    """Pixel scale (arcsec) of one object's epochs, for its centroid prior.

    Each epoch's linear scale is ``sqrt(|det J|)``: the side of a square
    pixel of the same sky area, blind to the CCD's rotation and flip. The
    object's scale is the mean over its epochs, because the centroid prior is
    one Gaussian in the sky frame those epochs share.

    Parameters
    ----------
    jacobs : list of galsim.JacobianWCS
        The object's per-epoch Jacobians, from :func:`get_galsim_jacobian`.

    Returns
    -------
    float
        Pixel scale in arcsec.

    @sc [decision:shape_measurement.fit_priors]
    """
    return float(np.mean([np.sqrt(abs(jac.pixelArea())) for jac in jacobs]))


def get_noise(gal, weight, guess, pixel_scale, thresh=1.2):
    """Get Noise.
    TO DO: modify guess, pixel scale
    Compute the sigma of the noise from an object postage stamp.
    Use a guess on the object size, ellipticity and flux to create a window
    function.

    Parameters
    ----------
    gal : numpy.ndarray
        Galaxy image
    weight : numpy.ndarray
        Weight image
    guess : list
        Gaussian parameters fot the window function
        ``[x0, y0, g1, g2, T, flux]``
    pixel_scale : float
        Pixel scale of the galaxy image
    thresh : float, optional
        Threshold to cut the window function,
        cut = ``thresh`` * sigma_noise;  the default is ``1.2``

    Returns
    -------
    float
        Sigma of the noise on the galaxy image

    """
    img_shape = gal.shape
    m_weight = weight != 0

    sig_tmp = sigma_mad(gal[m_weight])

    gauss_win = galsim.Gaussian(sigma=cs_size.T_to_sigma(guess[4]), flux=guess[5])
    gauss_win = gauss_win.shear(g1=guess[2], g2=guess[3])
    gauss_win = gauss_win.drawImage(
        nx=img_shape[0], ny=img_shape[1], scale=pixel_scale
    ).array

    m_weight = weight[gauss_win < thresh * sig_tmp] != 0

    sig_noise = sigma_mad(gal[gauss_win < thresh * sig_tmp][m_weight])

    return sig_noise


def central_seg_label(seg):
    """Centre-pixel label of a seg stamp — a *diagnostic*, not the central id.

    The authoritative central label is the object's SExtractor ``NUMBER``
    (``obj_id``), which is plumbed straight through: segmentation labels ARE
    the NUMBERs of the same SE run, so ``obj_id`` names the central footprint
    directly (shapepipe#776). This helper reads the centre pixel only as a
    registration cross-check — a few-pixel coadd-vs-epoch overlay offset can
    put a NEIGHBOUR's label on the centre pixel, which is exactly why the
    centre pixel must NOT define the central object (trusting it there would
    invert the uberseg mask). See :meth:`Ngmix._check_central_seg_label`.

    Parameters
    ----------
    seg : numpy.ndarray
        Segmentation stamp (integer SExtractor labels; 0 for sky).

    Returns
    -------
    int
        Label at the centre pixel — the central object's segmentation id.

    Raises
    ------
    ValueError
        If the centre pixel is sky (label 0): the cutout is not centred on any
        detected footprint (bad centroid / dropped detection), so uberseg
        cannot define a central object. Fail fast and loud; the object is then
        dropped by the per-object ``try/except`` in :meth:`Ngmix.process`.
    """
    cy, cx = seg.shape[0] // 2, seg.shape[1] // 2
    label = int(seg[cy, cx])
    if label == 0:
        raise ValueError(
            "central_seg_label: centre pixel is sky (label 0); the"
            + " segmentation stamp is not centred on a detected object."
        )
    return label


def seg_has_neighbour(seg, object_number):
    """Whether a seg stamp holds any non-central, non-zero label.

    The per-object neighbour flag (shapepipe#776, a systematics-test hook):
    ``True`` when the coadd segmentation stamp contains a footprint other than
    the central object's, i.e. the object is blended.

    Parameters
    ----------
    seg : numpy.ndarray
        Segmentation stamp (integer SExtractor labels; 0 for sky).
    object_number : int
        Central object's segmentation label.

    Returns
    -------
    bool
        ``True`` if any non-zero label differs from ``object_number``.
    """
    labels = seg[seg != 0]
    return bool(labels.size and np.any(labels != object_number))


def uberseg_mask(seg, object_number, dilate_neighbour=0):
    """Neighbour-side pixels of a stamp: the UberSeg mask.

    UberSeg is Erin Sheldon's MEDS/ngmix neighbour mask (``esheldon/meds``,
    https://github.com/esheldon/meds — ``MEDS.get_uberseg`` /
    ``meds._uberseg.uberseg_tree``, in MEDS itself "adapted from Niall MacCrann
    and Joe Zuntz", and used in the DES shear pipeline). Each stamp pixel is
    assigned to the object whose segmentation footprint it lies nearest to — a
    nearest-segment Voronoi partition — and the pixels assigned to a neighbour
    are masked.

    This function **reimplements** the partition rather than depending on
    ``meds``: it is a five-line ``scipy.spatial.cKDTree`` nearest-neighbour
    query, and pulling in the whole MEDS-file library — plus ``esutil`` and a
    C-extension build (``fitsio`` is already in the ShapePipe stack; ``esutil``
    is not) — for one standard geometric operation is disproportionate. The
    reimplementation was validated bit-for-bit against Sheldon's own reference
    ``get_uberseg`` over a battery of synthetic and random multi-object seg
    maps: the two masks agree on every pixel except exact Voronoi-boundary ties
    (a measure-zero tie-break convention — ``argmin`` first-min vs cKDTree
    order — scientifically inert). See the ``uberseg`` fiber's ``validation/``
    harness for the equivalence proof and figures. The cKDTree stands in for
    the C k-d tree; results are identical up to that boundary convention.

    Because the partition is by distance to the nearest footprint, the pixels
    surviving around a compact central object form a single connected,
    roughly circular core; the "circularisation" is emergent geometry, not a
    separate aperture. :func:`prepare_ngmix_weights` zeroes the weight of the
    masked pixels and keeps their image values, so metacal shears the
    neighbour's light along with the target's.

    Parameters
    ----------
    seg : numpy.ndarray
        Segmentation stamp: 0 for sky, the SExtractor object number for each
        detected object's footprint.
    object_number : int
        Segmentation label of the central object — its SExtractor ``NUMBER``
        (``obj_id``), authoritative because seg labels are the NUMBERs of the
        same SE run. Not the centre-pixel label, which a few-pixel coadd-vs-
        epoch offset can steal for a neighbour and so invert this mask.
    dilate_neighbour : int, optional
        Enlarge the mask by this many binary-dilation iterations of the
        neighbour footprints (4-connected, ~one pixel per iteration) on top of
        the Voronoi partition, to absorb the few-pixel coadd-vs-epoch
        registration offset the coadd-seg overlay accepts (shapepipe#776,
        decision on seg source). ``0`` (the default) is the pure Sheldon
        UberSeg mask. The dilation only ever adds pixels, so over-masking a
        boundary pixel costs a little central-object S/N but never leaks
        neighbour flux into the fit.

    Returns
    -------
    numpy.ndarray of bool
        ``True`` on pixels nearer a neighbour's footprint than the central
        object's; all ``False`` when the stamp holds no neighbour.
    @sc [decision:shape_measurement.blend_handling]
    """
    seg = np.asarray(seg)
    masked = np.zeros(seg.shape, dtype=bool)

    obj_pix = np.argwhere(seg != 0)
    labels = seg[seg != 0]
    # No neighbour footprint on the stamp: nothing to mask (cf. MEDS' early
    # ``len(np.unique(seg)) == 2`` return).
    if np.all(labels == object_number):
        return masked

    # Nearest segmentation pixel for every stamp pixel; mask wherever that
    # nearest footprint belongs to a neighbour rather than to the central
    # object.
    grid = np.indices(seg.shape).reshape(2, -1).T
    _, nearest = cKDTree(obj_pix).query(grid)
    masked |= labels[nearest].reshape(seg.shape) != object_number

    if dilate_neighbour > 0:
        masked |= binary_dilation(
            (seg != 0) & (seg != object_number),
            iterations=dilate_neighbour,
        )

    return masked


def defect_mask(weight, flag, bkg_rms=None):
    """Defect pixels of one epoch stamp.

    @sc [decision:shape_measurement.defect_fill,label:physics] defect-set-unsymmetrized
    A defect is a pixel with zero exposure weight (off-tile pixels
    included), a nonzero exposure flag, or, when a background RMS map is
    given, a non-finite or non-positive RMS. This one set is zero-weighted
    and filled by :func:`prepare_ngmix_weights` and counted by the epoch
    cuts. It is not ORed with its rotations: a four-fold fill quadruples the
    filled area near the object and with it m (astra option
    ``defect_fill.symmetrized_4fold_noise``).

    Parameters
    ----------
    weight : numpy.ndarray
        Exposure weight stamp.
    flag : numpy.ndarray
        Exposure flag stamp.
    bkg_rms : numpy.ndarray, optional
        Background RMS stamp.

    Returns
    -------
    numpy.ndarray of bool
        ``True`` on defect pixels.
    """
    defect = (weight == 0) | (flag != 0)
    if bkg_rms is not None:
        defect |= ~(np.isfinite(bkg_rms) & (bkg_rms > 0))
    return defect


def central_defect_vetoes(defect, interpolated):
    """Whether a defect near the stamp centre drops the epoch.

    @sc [decision:shape_measurement.central_defect_veto,decision:shape_measurement.defect_fill] veto-radius-follows-the-fill
    The epoch is dropped when an interpolated defect pixel lies closer than
    ``EPOCH_INTERPOLATED_DEFECT_RADIUS`` to the stamp centre, or a
    noise-filled one closer than ``EPOCH_CENTRAL_DEFECT_RADIUS``: the
    smallest radii at which the defects kept recover shear within
    |m| < 1% and |c| < 5e-4. ``interpolated`` is the set
    :func:`prepare_ngmix_weights` interpolates, so each pixel is vetoed at
    the radius of the fill it gets. The veto reads only masks, so it selects
    on nothing that responds to shear. Calibration: astra decision
    ``shape_measurement.central_defect_veto``; guarded by
    ``tests/science/test_defect_recovery.py``.

    Parameters
    ----------
    defect : numpy.ndarray of bool
        Defect mask of one epoch stamp (:func:`defect_mask`).
    interpolated : numpy.ndarray of bool
        The defect pixels that are interpolated (:func:`interpolable_defects`
        with the epoch's removed neighbour pixels); the rest are noise-filled.

    Returns
    -------
    bool
        ``True`` if the epoch should be dropped.
    """
    rows, cols = np.indices(defect.shape)
    distance = np.hypot(
        rows - (defect.shape[0] - 1) / 2, cols - (defect.shape[1] - 1) / 2
    )
    radius = np.where(
        interpolated,
        EPOCH_INTERPOLATED_DEFECT_RADIUS,
        EPOCH_CENTRAL_DEFECT_RADIUS,
    )
    return bool(np.any(defect & (distance < radius)))


def prepare_ngmix_weights(
    gal, weight, flag, rng, bkg_rms=None,
    blend_handling="noisefill", seg=None, object_number=None,
    dilate_neighbour=0, neighbour=None,
):
    """Build one epoch's image, weight map and noise image for ngmix.

    Every stamp pixel falls in one of three classes. A clean pixel keeps its
    light and its weight. A pixel whose light is replaced by noise keeps
    neither: a defect that is not interpolated, or under ``"noisefill"`` a
    marked neighbour pixel. A pixel whose light stays at zero weight is a
    hole in the likelihood only: an interpolated defect and its quarter-turn
    copies, or under ``"uberseg"`` a neighbour-side pixel.

    @sc [decision:shape_measurement.defect_fill,decision:shape_measurement.blend_handling] defects-filled-whatever-the-blend-handling
    Every pixel of :func:`defect_mask` gets weight 0 and is filled the same
    way under every ``blend_handling``. The short runs of
    :func:`interpolable_defects` take a Clough-Tocher interpolant of the
    pixels whose light the image keeps (:func:`interpolate_defects`), in the
    image and the noise image alike; the other defects take noise at the
    background RMS. Metacal shears the whole image whatever the weights, so
    no raw defect value may reach it.

    @sc [decision:shape_measurement.blend_handling] noisefill-fills-markers
    Under ``"noisefill"`` the ``neighbour`` pixels (the tile VIGNET's -1e30
    neighbour markers) get weight 0 and noise. Their light is gone, so they
    neither support the interpolant nor receive it.

    @sc [decision:shape_measurement.blend_handling] uberseg-ignores-markers
    Under ``"uberseg"`` the markers are ignored: the pixels of
    :func:`uberseg_mask` get weight 0 and keep their light, which is real
    sky that metacal shears with the target and which supports the
    interpolant.

    @sc [decision:shape_measurement.weight_symmetrization,label:physics] symmetrized-weight-holes
    A one-sided zero-weight hole pulls the Gaussian fit toward its side, so
    the weight is also zeroed on the three quarter-turn copies, about the
    stamp centre, of every interpolated pixel, and their light stays.
    Noise-filled pixels and the uberseg neighbour side are not symmetrized.
    Measurements: astra decision ``shape_measurement.weight_symmetrization``;
    guarded by ``tests/science/test_defect_recovery.py``.

    Parameters
    ----------
    gal : numpy.ndarray
        Background-subtracted galaxy stamp.
    weight : numpy.ndarray
        Exposure weight stamp; zero marks a defect.
    flag : numpy.ndarray
        Exposure flag stamp; nonzero marks a defect.
    rng : numpy.random.RandomState
        Random state for the noise realisations (seeded per object; see
        :func:`position_seed`).
    bkg_rms : numpy.ndarray, optional
        Per-pixel background RMS map. If supplied, clean pixels use
        ``1 / bkg_rms**2`` as the ngmix inverse variance, and non-finite or
        non-positive values mark defects. Otherwise every clean pixel gets
        ``1 / sigma_mad(gal)**2``.
    blend_handling : {"noisefill", "uberseg"}, optional
        Neighbour treatment. ``"noisefill"`` (default) zero-weights and
        noise-fills the ``neighbour`` pixels. ``"uberseg"`` ignores
        ``neighbour``, zeroes the weight of the pixels of
        :func:`uberseg_mask` and keeps their raw image values.
    seg : numpy.ndarray, optional
        Segmentation map on the stamp grid (object NUMBERs). Required for
        ``blend_handling="uberseg"``; ignored otherwise.
    object_number : int, optional
        Central object's segmentation label. Required for
        ``blend_handling="uberseg"``; ignored otherwise.
    dilate_neighbour : int, optional
        Neighbour-mask dilation iterations, passed to :func:`uberseg_mask`
        under ``blend_handling="uberseg"``; ignored otherwise.
    neighbour : numpy.ndarray of bool, optional
        The epoch's neighbour mask (``Postage_stamp.neighbours``); read only
        under ``blend_handling="noisefill"``. ``None`` marks no pixel.

    Returns
    -------
    numpy.ndarray
        Galaxy image with defect pixels, and under noisefill marked
        neighbour pixels, filled.
    numpy.ndarray
        Inverse-variance weight map for ngmix.
    numpy.ndarray
        Noise image: an independent realisation over the whole stamp, for
        metacal's ``fixnoise``, interpolated where the galaxy image is.

    Raises
    ------
    ValueError
        If ``blend_handling`` is unknown, or ``"uberseg"`` lacks ``seg`` or
        ``object_number``.
    RuntimeError
        If the interpolant does not reach a pixel :func:`interpolable_defects`
        selected (degenerate support, or non-finite image values in it).
    @sc [decision:masking.pixel_mask_source,decision:shape_measurement.blend_handling,decision:shape_measurement.defect_fill,decision:shape_measurement.galaxy_pixel_weights,decision:shape_measurement.weight_symmetrization]
    """
    if blend_handling not in BLEND_HANDLINGS:
        raise ValueError(
            f"Unknown blend_handling '{blend_handling}'; expected one of"
            + f" {BLEND_HANDLINGS}"
        )
    if blend_handling == "uberseg" and (seg is None or object_number is None):
        raise ValueError(
            "blend_handling='uberseg' requires a segmentation map and the"
            + " central object_number; none reached prepare_ngmix_weights."
            + " The tile catalogue's SEG_VIGNET column carries the map (see"
            + " CosmoStat/shapepipe#776)."
        )

    defect = defect_mask(weight, flag, bkg_rms)
    no_pixel = np.zeros_like(defect)
    # Neighbour pixels whose light noisefill replaces (noisefill-fills-markers)
    # or whose weight uberseg zeroes (uberseg-ignores-markers).
    removed_neighbour, neighbour_side = no_pixel, no_pixel
    if blend_handling == "noisefill" and neighbour is not None:
        removed_neighbour = np.asarray(neighbour, dtype=bool)
    elif blend_handling == "uberseg":
        neighbour_side = uberseg_mask(seg, object_number, dilate_neighbour)
    # Pixels whose light the image keeps raw.
    clean = ~(defect | removed_neighbour)
    weighted = clean & ~neighbour_side

    if bkg_rms is None:
        sig_noise = sigma_mad(gal)
        # Guard the degenerate constant stamp (sigma_mad == 0): 0 * inf
        # would otherwise put NaN in a fully-masked weight map.
        weight_map = (
            weighted.astype(float) / sig_noise ** 2
            if sig_noise > 0
            else np.zeros_like(gal, dtype=float)
        )
    else:
        weight_map = np.zeros_like(gal, dtype=float)
        weight_map[weighted] = 1.0 / bkg_rms[weighted] ** 2
        # Per-pixel noise sigma for the realisations below: metacal's
        # fixnoise bookkeeping (1/w + 1/w_noise) assumes the noise image
        # is a faithful realisation of the per-pixel variance the weights
        # claim; a scalar sigma there mis-reports errors and erodes the
        # inverse-variance advantage whenever the RMS map actually varies.
        # Pixels without a valid RMS take the median over clean pixels.
        valid_rms = np.isfinite(bkg_rms) & (bkg_rms > 0)
        sig_noise = (
            np.where(valid_rms, bkg_rms, np.median(bkg_rms[clean]))
            if clean.any()
            else sigma_mad(gal)
        )

    # Guard: sig_noise=0 means galaxy is at the CCD edge (mostly-zero stamp).
    # Division by zero would make weight_map NaN/inf, crashing GalSim C code.
    # np.all keeps the guard valid for the per-pixel bkg_rms path (any zero
    # pixel would produce the same NaN weight the scalar case guards against).
    if not np.all(sig_noise > 0):
        return np.zeros_like(gal), np.zeros_like(weight_map), np.zeros_like(gal)

    noise_img = rng.standard_normal(gal.shape) * sig_noise
    noise_img_gal = rng.standard_normal(gal.shape) * sig_noise
    gal_filled = np.where(clean, gal, noise_img_gal).astype(gal.dtype)
    interpolated = interpolable_defects(defect, removed_neighbour)
    if interpolated.any():
        # One operator for the image and the noise image, so fixnoise
        # mirrors the science image's interpolated noise.
        filled = interpolate_defects(
            [gal, noise_img], ~clean, interpolated
        )[:, interpolated]
        if not np.all(np.isfinite(filled)):
            raise RuntimeError(
                "The defect interpolant is not finite on pixels"
                + " interpolable_defects selected."
            )
        gal_filled[interpolated] = filled[0]
        noise_img[interpolated] = filled[1]
        weight_map[fourfold(interpolated)] = 0.0

    return gal_filled, weight_map, noise_img


def make_ngmix_observation(
    gal, weight, flag, psf, wcs, rng,
    bkg_rms=None, centroid_source="wcs", offset=None,
    blend_handling="noisefill", seg=None, object_number=None,
    dilate_neighbour=0, neighbour=None,
):
    """Build an ngmix Observation for a single galaxy epoch.

    The galaxy Jacobian origin sets where the centroid prior is centered, so
    it must sit on the object. Two ways to place it, selected by
    ``centroid_source``:

    * ``"wcs"`` (default) — the **coadd centroid**: place the origin at the
      object's catalogue sky position as a sub-pixel ``offset`` from the
      stamp center, with no shape measurement. The offset is not recomputed
      here — it is the value the stamp extractor already used to round the
      extraction pixel
      (:func:`shapepipe.modules.vignetmaker_package.vignetmaker.get_stamps`),
      propagated on the vignette. One projection, one rounding: the
      extraction and the centroid prior cannot disagree near a rounding tie.
      Stable for both galaxies and stars, and correct when the object sits
      off the stamp center.
    * ``"hsm"`` — re-center on the HSM adaptive-moment centroid measured from
      the stamp pixels, following the light rather than the astrometry.
      Noisier, notably for stars; the option for stamps that carry no
      propagated offset.

    Parameters
    ----------
    gal : numpy.ndarray
    weight : numpy.ndarray
    flag : numpy.ndarray
    psf : numpy.ndarray
    wcs : galsim.BaseWCS
        Local WCS Jacobian at the object position.
    rng : numpy.random.RandomState
        Random state for the noise realisations (seeded per object; see
        :func:`position_seed`).
    bkg_rms : numpy.ndarray, optional
        Per-pixel background RMS map.
    centroid_source : {"wcs", "hsm"}, optional
        How to place the galaxy Jacobian origin; the default is ``"wcs"``.
    offset : array_like, optional
        Sub-pixel ``[row, col]`` coadd-centroid offset propagated from the
        stamp extractor. Required for ``centroid_source="wcs"`` (ignored for
        ``"hsm"``).
    blend_handling : {"noisefill", "uberseg"}, optional
        Neighbour treatment passed through to :func:`prepare_ngmix_weights`;
        the default ``"noisefill"`` zero-weights and noise-fills the
        ``neighbour`` pixels.
    seg : numpy.ndarray, optional
        Segmentation map on the stamp grid. Required for
        ``blend_handling="uberseg"`` (ignored otherwise).
    object_number : int, optional
        Central object's segmentation label. Required for
        ``blend_handling="uberseg"`` (ignored otherwise).
    dilate_neighbour : int, optional
        Neighbour-mask dilation iterations passed through to
        :func:`prepare_ngmix_weights` under ``blend_handling="uberseg"``.
    neighbour : numpy.ndarray of bool, optional
        Neighbour mask passed through to :func:`prepare_ngmix_weights`.

    Returns
    -------
    ngmix.observation.Observation
    @sc [decision:shape_measurement.centroid_source,decision:shape_measurement.psf_likelihood_noise]
    """
    psf_jacob = ngmix.Jacobian(
        row=(psf.shape[0] - 1) / 2,
        col=(psf.shape[1] - 1) / 2,
        wcs=wcs,
    )
    # A model fit always needs a weight: without one ngmix defaults to unit
    # weights on a unit-flux PSF stamp, so the galaxy GPriorBA(0.4) prior
    # handed to the diagnostic PSF fitter swamps the (tiny) likelihood and the
    # recovered PSF shape collapses toward the prior (#749). PSFEx/MCCD models
    # already carry a little noise, so a finite flat weight matching the
    # esheldon/aguinot noise budget (psf_noise ~ 1e-5) suffices to restore the
    # likelihood — no explicit stamp-noise injection needed (#774). This is the
    # PSF analogue of the galaxy weight built just below; force float so an
    # integer-typed stamp cannot truncate the weight (matching that build).
    # Metacal's own PSF fit is prior-free (AdmomFitter) and so is insensitive
    # to this flat weight's scale, leaving the calibration untouched.
    psf_wt = np.full_like(psf, 1.0 / PSF_NOISE ** 2, dtype=float)
    psf_obs = Observation(psf, weight=psf_wt, jacobian=psf_jacob)

    gal_masked, weight_map, noise_img = prepare_ngmix_weights(
        gal, weight, flag, rng, bkg_rms=bkg_rms,
        blend_handling=blend_handling, seg=seg, object_number=object_number,
        dilate_neighbour=dilate_neighbour, neighbour=neighbour,
    )

    if centroid_source == "hsm":
        # Re-center the Jacobian on the HSM adaptive-moment centroid (pixel
        # offset from the stamp center), measured on the filled image so no
        # raw defect value pulls it; fall back to the stamp center if HSM
        # fails.
        try:
            _hsm = galsim.hsm.FindAdaptiveMom(
                galsim.Image(gal_masked, scale=1.0), strict=False
            )
            if _hsm.error_message != "":
                raise galsim.hsm.GalSimHSMError(_hsm.error_message)
            _cen = _hsm.moments_centroid - galsim.Image(gal, scale=1.0).center
            cen_row, cen_col = _cen.y, _cen.x
        except Exception:
            cen_row, cen_col = 0.0, 0.0
    elif centroid_source == "wcs":
        # Coadd centroid: use the sub-pixel offset the stamp extractor already
        # computed and rounded against (propagated on the vignette), rather
        # than re-projecting the sky position through the WCS and re-rounding:
        # one projection, one rounding, so a milli-pixel WCS disagreement
        # cannot flip a rounding tie and put the prior a whole pixel off.
        if offset is None:
            raise ValueError(
                "centroid_source='wcs' requires the coadd-centroid offset "
                + "propagated from the stamp extractor (the vignette's "
                + "OFFSET), but none was given: re-extract the stamps with "
                + "the current vignetmaker, or opt in to "
                + "centroid_source='hsm'"
            )
        cen_row, cen_col = float(offset[0]), float(offset[1])
    else:
        raise ValueError(
            f"Unknown centroid_source '{centroid_source}'; expected"
            + " 'hsm' or 'wcs'"
        )

    gal_jacob = ngmix.Jacobian(
        row=(gal.shape[0] - 1) / 2 + cen_row,
        col=(gal.shape[1] - 1) / 2 + cen_col,
        wcs=wcs,
    )

    return Observation(
        gal_masked,
        weight=weight_map,
        jacobian=gal_jacob,
        psf=psf_obs,
        noise=noise_img,
    )

def _average_psf_fits(results_and_weights):
    """Weight-average a set of per-epoch ngmix PSF-fit results.

    Shared core for both PSF families this module exports: the metacal
    reconvolution kernel (:func:`average_multiepoch_psf`) and the original
    image PSF (:func:`average_original_psf`). Epochs whose PSF fit failed
    (``flags != 0``, carrying only flags/pars and no T/g) are dropped.

    Parameters
    ----------
    results_and_weights : iterable of (dict, float)
        Per-epoch ``(result, weight)`` pairs, where ``result`` is an ngmix
        Fitter result with keys ``flags``, ``g``, ``g_err``, ``T``,
        ``T_err`` and ``weight`` is the epoch's averaging weight.

    Returns
    -------
    dict
        Keys ``g_psf``, ``g_psf_err``, ``T_psf``, ``T_psf_err`` (weighted
        averages over the surviving epochs) and ``n_epoch`` (their count).
    @sc [decision:shape_measurement.psf_epoch_averaging]
    """
    n_epoch_used = 0
    wsum = 0
    g_psf_sum = np.array([0., 0.])
    g_psf_err_sum = np.array([0., 0.])
    T_psf_sum = 0
    T_psf_err_sum = 0
    for result, weight in results_and_weights:
        if result['flags'] != 0:
            continue
        n_epoch_used += 1
        wsum += weight
        g_psf_sum += result['g'] * weight
        g_psf_err_sum += result['g_err'] * weight
        T_psf_sum += result['T'] * weight
        T_psf_err_sum += result['T_err'] * weight

    if wsum == 0:
        raise ZeroDivisionError('Sum of weights = 0, division by zero')

    return {
        'g_psf': g_psf_sum / wsum,
        'g_psf_err': g_psf_err_sum / wsum,
        'T_psf': T_psf_sum / wsum,
        'T_psf_err': T_psf_err_sum / wsum,
        'n_epoch': n_epoch_used,
    }


def average_multiepoch_psf(obsdict):
    """Average the metacal *reconvolution* PSF over epochs.

    The PSF carried by each metacal observation (``obs.psf``) is the
    Gaussian reconvolution kernel that metacal fit and convolved back in —
    round by construction and slightly enlarged relative to the original
    PSF. This is the kernel defining the sheared galaxy images, exported to
    the reconvolution-kernel columns
    (``NGMIX_G1/G2_PSF_RECONV``, ``NGMIX_T_PSF_RECONV``). The independent fit
    of the *original* image PSF is :func:`average_original_psf`.

    Parameters
    ----------
    obsdict : dict
        Observation dict returned by MetacalBootstrapper.go().

    Returns
    -------
    dict
        Keys: 'g_psf', 'g_psf_err', 'T_psf', 'T_psf_err' (weighted
        averages over the epochs whose PSF fit succeeded) and 'n_epoch'
        (the number of those surviving epochs).

    @sc [decision:shape_measurement.psf_epoch_averaging]
    """
    # ignore_failed_psf=True drops failed-PSF epochs from the galaxy fit but
    # keeps them in obsdict; _average_psf_fits skips them on flags != 0.
    return _average_psf_fits(
        (obs.psf.meta['result'], obs.weight.sum())
        for obs in obsdict['noshear']
    )


def average_original_psf(gal_obs_list, psf_runner):
    """Fit and average the *original* image PSF over epochs.

    The original PSF is the psfex/mccd model stamp handed to ngmix
    (``gal_obs.psf``), fit here with the same ``psf_runner`` (hence the same
    fit prior and guesser) used inside metacal, but on the PSF *before*
    metacal's reconvolution. Exported to the original-PSF columns
    (``NGMIX_G1/G2_PSF_ORIG``, ``NGMIX_T_PSF_ORIG``). Distinct from the
    reconvolution-kernel fit (:func:`average_multiepoch_psf`): the original
    PSF retains its true ellipticity and size, whereas the reconvolution
    kernel is round and enlarged by construction. This is the PSF whose true
    shape and size enter object-wise PSF-leakage diagnostics.

    Weighted by the raw galaxy inverse variance (``gal_obs.weight.sum()`` per
    epoch); :func:`average_multiepoch_psf` uses the same scheme but the
    fixnoise-combined metacal-image weight instead.

    The fit runs on a *copy* of each PSF observation so ``gal_obs.psf`` —
    the object metacal later deep-copies and consumes via
    ``boot.go(gal_obs_list)`` — is never mutated. ``PSFRunner.go`` sets both
    ``.meta['result']`` and, on success, the ``.gmix`` attribute of the
    observation it fits; were that ``gal_obs.psf`` itself, the stray gmix
    would survive metacal's deep copy and be reused as the
    ``MetacalFitGaussPSF`` fallback when its own admom+ML PSF fits both fail,
    silently rescuing objects the base branch dropped (``BootPSFFailure``) and
    changing the galaxy/shear result set. Fitting a copy closes this
    PSF-aliasing channel; combined with :func:`do_ngmix_metacal` seeding this
    pre-fit from a *snapshot* of the metacal RNG (so the pre-fit does not
    advance it), the add-column refactor stays bit-identical on the galaxy
    results.

    Parameters
    ----------
    gal_obs_list : ngmix.observation.ObsList
        Per-epoch galaxy observations; each ``gal_obs.psf`` is the original
        (pre-metacal) PSF observation to fit, with no further ``.psf`` of
        its own so the runner fits the stamp itself. Left pristine.
    psf_runner : ngmix.runners.PSFRunner
        The module's PSF runner (built by :func:`make_runners` from the
        shared ``prior``).

    Returns
    -------
    dict
        Same keys as :func:`average_multiepoch_psf`.
    @sc [decision:shape_measurement.psf_epoch_averaging]
    """
    def fit(gal_obs):
        # Fit a COPY so gal_obs.psf stays pristine for metacal — see docstring.
        # Failed fits keep flags != 0 and are dropped by _average_psf_fits.
        psf_obs = gal_obs.psf.copy()
        psf_runner.go(psf_obs)
        return psf_obs.meta['result'], gal_obs.weight.sum()

    return _average_psf_fits(fit(gal_obs) for gal_obs in gal_obs_list)


def make_runners(prior, flux_guess, rng):
    """Build the module's galaxy and PSF runners.

    Single source of truth for the fitter configuration (Gaussian galaxy
    and PSF models, guessers, retry counts), shared by the metacal
    bootstrap below and by validation tests that fit module-built
    observations directly.

    Parameters
    ----------
    prior : ngmix.joint_prior.PriorSimpleSep
        Priors for the fitting parameters.
    flux_guess : float
        Initial flux guess.
    rng : numpy.random.RandomState
        Random state for the guessers.

    Returns
    -------
    tuple
        (runner, psf_runner) : ngmix.runners.Runner, ngmix.runners.PSFRunner
    @sc [decision:shape_measurement.fit_initialisation,decision:shape_measurement.fit_priors,decision:shape_measurement.galaxy_model]
    """
    fitter = ngmix.fitting.Fitter(model='gauss', prior=prior)
    guesser = ngmix.guessers.TPSFFluxAndPriorGuesser(rng=rng, T=0.25, prior=prior)

    psf_fitter = ngmix.fitting.Fitter(model='gauss', prior=prior)
    psf_guesser = ngmix.guessers.TFluxGuesser(rng=rng, T=0.25, prior=prior, flux=flux_guess)

    return (
        ngmix.runners.Runner(fitter=fitter, guesser=guesser, ntry=5),
        ngmix.runners.PSFRunner(fitter=psf_fitter, guesser=psf_guesser, ntry=2),
    )


def do_ngmix_metacal(
    stamp, prior, flux_guess, rng, centroid_source="wcs",
    blend_handling="noisefill", object_number=None, dilate_neighbour=0,
    metacal_psf="fitgauss",
):
    """Do Ngmix Metacal.

    Performs metacalibration on a single multi-epoch object and returns the
    joint shape measurement with NGMIX.

    Parameters
    ----------
    stamp : Postage_stamp
        Postage stamps for all epochs of one galaxy.
    prior : ngmix.joint_prior.PriorSimpleSep
        Priors for the fitting parameters.
    flux_guess : float
        Initial flux guess.
    rng : numpy.random.RandomState
        Random state for guesses and priors.
    centroid_source : {"wcs", "hsm"}, optional
        How to place the galaxy Jacobian origin; passed through to
        :func:`make_ngmix_observation`. The default is ``"wcs"`` (the
        coadd-centroid offset in ``stamp.offsets``, propagated from the stamp
        extractor); ``"hsm"`` uses the adaptive-moment centroid from the
        stamp pixels — see that function.
    blend_handling : {"noisefill", "uberseg"}, optional
        Neighbour treatment passed through to
        :func:`make_ngmix_observation`; the default ``"noisefill"``
        zero-weights and noise-fills the pixels of ``stamp.neighbours``.
        ``"uberseg"`` consumes ``stamp.segs`` and ``object_number``.
    object_number : int, optional
        Central object's segmentation label — its SExtractor ``NUMBER``
        (``obj_id``), authoritative because seg labels are the NUMBERs of the
        same SE run (shapepipe#776). Required for ``blend_handling="uberseg"``
        (ignored otherwise).
    dilate_neighbour : int, optional
        Neighbour-mask dilation iterations passed through to
        :func:`make_ngmix_observation` under ``blend_handling="uberseg"``.
    metacal_psf : {"fitgauss", "gauss", "dilate", "azgauss"}, optional
        The metacal reconvolution-kernel scheme passed to ngmix as
        ``metacal_pars['psf']``. The default ``"fitgauss"`` fits a Gaussian to
        the PSF and rounds it (``MetacalFitGaussPSF``); ``"gauss"`` reconvolves
        with a fixed round Gaussian sized from the PSF (``MetacalGaussPSF``);
        ``"dilate"`` dilates the original PSF; ``"azgauss"`` (ngmix >= 2.4.1) is
        a noise-robust variant of ``"gauss"``. This kernel is a SEPARATE object
        from the PSF-model fit (the ``psf_runner`` above) — it only sets the
        round PSF that metacal reconvolves with after shearing, so it moves the
        metacal *response* (and therefore the recovered shear) but never reaches
        the deconvolution, which is by the PSF image.

    Returns
    -------
    MetacalResult
        Named 3-tuple ``(resdict, reconv, orig)``: the MetacalBootstrapper
        result dict, the averaged metacal *reconvolution*-kernel PSF dict
        (:func:`average_multiepoch_psf`), and the averaged *original* image-PSF
        dict (:func:`average_original_psf`). The two PSF dicts share keys but
        describe different PSFs; the named fields guard against transposing
        them. Unpacks positionally as ``resdict, psf_res, psf_orig_res``.
    @sc [decision:shape_measurement.defect_fill,decision:shape_measurement.metacal_scheme,decision:shape_measurement.weight_symmetrization]
    """
    n_epoch = len(stamp.gals)
    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    gal_obs_list = ObsList()
    for n_e in range(n_epoch):
        bkg_rms = stamp.bkg_rms[n_e] if len(stamp.bkg_rms) > n_e else None
        gal_obs = make_ngmix_observation(
            stamp.gals[n_e],
            stamp.weights[n_e],
            stamp.flags[n_e],
            stamp.psfs[n_e],
            stamp.jacobs[n_e],
            rng,
            bkg_rms=bkg_rms,
            centroid_source=centroid_source,
            offset=stamp.offsets[n_e] if n_e < len(stamp.offsets) else None,
            blend_handling=blend_handling,
            seg=stamp.segs[n_e] if n_e < len(stamp.segs) else None,
            object_number=object_number,
            dilate_neighbour=dilate_neighbour,
            neighbour=(
                stamp.neighbours[n_e] if n_e < len(stamp.neighbours) else None
            ),
        )
        gal_obs_list.append(gal_obs)

    runner, psf_runner = make_runners(prior, flux_guess, rng)

    # Fit the ORIGINAL (psfex/mccd) PSF before metacal reconvolves it. Use a
    # psf_runner seeded from a snapshot of rng's state so this prefit does NOT
    # advance the rng that boot.go consumes below — keeping the galaxy/shear
    # results bit-identical to the no-prefit branch. average_original_psf fits a
    # COPY of each gal_obs.psf, so gal_obs_list also reaches boot.go pristine.
    prefit_rng = np.random.RandomState()
    prefit_rng.set_state(rng.get_state())
    _, prefit_psf_runner = make_runners(prior, flux_guess, prefit_rng)
    psf_orig_res = average_original_psf(gal_obs_list, prefit_psf_runner)

    metacal_pars = {
        'types': ['noshear', '1p', '1m', '2p', '2m'],
        'step': 0.01,
        'psf': metacal_psf,
        'fixnoise': True,
        'use_noise_image': True,
    }

    boot = ngmix.metacal.MetacalBootstrapper(
        runner=runner,
        psf_runner=psf_runner,
        ignore_failed_psf=True,
        rng=rng,
        **metacal_pars,
    )
    resdict, obsdict = boot.go(gal_obs_list)
    psf_res = average_multiepoch_psf(obsdict)
    del obsdict

    # Each FitModel in resdict retains the full metacal observations
    # (``.obs``, a plain attribute set by FitModel._set_obs) plus the pixel
    # arrays derived from them (``._pixels_list``) — several MB per object.
    # Results accumulate until the SAVE_BATCH flush, and compile_results
    # reads only scalar dict keys, so release the arrays now to keep
    # resident memory flat.
    for _fm in resdict.values():
        if hasattr(_fm, 'obs'):
            _fm.obs = None
        if hasattr(_fm, '_pixels_list'):
            _fm._pixels_list = None

    return MetacalResult(resdict=resdict, reconv=psf_res, orig=psf_orig_res)
