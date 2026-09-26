"""DEFECT INTERPOLATION.

Clough-Tocher interpolation of short defect runs before metacal, selected
with ``DEFECT_FILL = interpolate`` (see
:func:`shapepipe.modules.ngmix_package.ngmix.prepare_ngmix_weights`).

"""

import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.ndimage import binary_dilation, label
from scipy.spatial import QhullError

# Longest row or column run of defect pixels that is interpolated.
MAX_INTERPOLATED_RUN = 3

# Clean pixels within this Chebyshev distance (pixels) of the interpolated
# pixels support the interpolant.
SUPPORT_RADIUS = 4

_ROW_RUNS = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]])


def _short_row_runs(defect, max_run):
    """Pixels in row runs of at most ``max_run`` defects that stop short of
    both stamp borders."""
    labels, n_runs = label(defect, structure=_ROW_RUNS)
    if n_runs == 0:
        return np.zeros_like(defect)
    short = np.bincount(labels.ravel(), minlength=n_runs + 1) <= max_run
    short[0] = False
    short[labels[:, 0]] = False
    short[labels[:, -1]] = False
    return short[labels]


def interpolable_defects(defect, max_run=MAX_INTERPOLATED_RUN):
    """Defect pixels that ``DEFECT_FILL = interpolate`` interpolates.

    @sc [decision:shape_measurement.defect_fill] interpolable-defects
    A defect pixel is interpolated when its row or its column run of defect
    pixels is at most ``max_run`` (3) long and has clean pixels at both
    ends. That covers columns, 3-px bleeds and isolated pixels, the widths
    whose shear recovery is calibrated. Wider holes, and runs that reach the
    stamp border (edge bands, corners), have clean light on one side only;
    they are noise-filled and vetoed at the noise-fill radius (see
    :func:`~shapepipe.modules.ngmix_package.ngmix.central_defect_vetoes`).
    The rule reads only the mask and commutes with quarter turns of the
    stamp.

    Parameters
    ----------
    defect : numpy.ndarray of bool
        Defect mask of one epoch stamp.
    max_run : int, optional
        Longest interpolated run; the default is ``MAX_INTERPOLATED_RUN``.

    Returns
    -------
    numpy.ndarray of bool
        ``True`` on the defect pixels to interpolate.
    """
    defect = np.asarray(defect, dtype=bool)
    return (
        _short_row_runs(defect, max_run)
        | _short_row_runs(defect.T, max_run).T
    )


def fourfold(mask):
    """Union of a square stamp mask with its quarter turns.

    Parameters
    ----------
    mask : numpy.ndarray of bool
        Square mask.

    Returns
    -------
    numpy.ndarray of bool
        ``mask`` ORed with its rotations by 90, 180 and 270 degrees about
        the stamp centre.

    Raises
    ------
    ValueError
        If ``mask`` is not square.
    """
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 2 or mask.shape[0] != mask.shape[1]:
        raise ValueError(
            f"A quarter-turn orbit needs a square stamp, not {mask.shape}"
        )
    return mask | np.rot90(mask) | np.rot90(mask, 2) | np.rot90(mask, 3)


def _interpolate_once(planes, defect, target):
    """One Clough-Tocher interpolant of every plane at ``target``; NaN
    elsewhere and where the support cannot reach."""
    out = np.full(planes.shape, np.nan)
    support = binary_dilation(
        target, structure=np.ones((3, 3), dtype=bool),
        iterations=SUPPORT_RADIUS,
    ) & ~defect
    points = np.argwhere(support).astype(float)
    if len(points) < 3:
        return out
    query = np.argwhere(target)
    try:
        interpolant = CloughTocher2DInterpolator(
            points, planes[:, support].T, fill_value=np.nan,
        )
    except QhullError:
        return out
    out[:, query[:, 0], query[:, 1]] = interpolant(query.astype(float)).T
    return out


def interpolate_defects(planes, defect, target):
    """Replace the ``target`` pixels of every plane by a Clough-Tocher
    interpolant of the clean pixels around them.

    @sc [decision:shape_measurement.defect_fill] shared-rotation-averaged-interpolant
    The support is the clean pixels within ``SUPPORT_RADIUS`` (4 px) of the
    target; no defect pixel enters it, so defect values are never read. For
    each quarter turn of the stamp, one Delaunay triangulation of the support
    serves every plane, so the science image and the metacal noise image
    see the same linear operator and fixnoise mirrors the science image's
    interpolated noise. A regular grid's triangulation has degenerate
    diagonals, so one orientation has a preferred direction; averaging the
    four quarter-turned operators makes the fill commute with rotations of
    the stamp. That removes up to 6e-5 of c2 for a column or 3-px bleed
    6 px from the object.

    Parameters
    ----------
    planes : array_like
        Stamp planes, shape ``(n, ny, nx)``.
    defect : numpy.ndarray of bool
        Every defect pixel, shape ``(ny, nx)``; none of them supports the
        interpolant.
    target : numpy.ndarray of bool
        Defect pixels to interpolate (:func:`interpolable_defects`).

    Returns
    -------
    numpy.ndarray
        A copy of ``planes`` with ``target`` pixels interpolated; NaN at a
        target pixel whose support is degenerate in some orientation.
    """
    planes = np.asarray(planes, dtype=float)
    defect = np.asarray(defect, dtype=bool)
    target = np.asarray(target, dtype=bool)
    out = planes.copy()
    if not target.any():
        return out
    turns = [
        np.rot90(
            _interpolate_once(
                np.rot90(planes, k, axes=(1, 2)), np.rot90(defect, k),
                np.rot90(target, k),
            ),
            -k, axes=(1, 2),
        )[:, target]
        for k in range(4)
    ]
    out[:, target] = np.mean(turns, axis=0)
    return out
