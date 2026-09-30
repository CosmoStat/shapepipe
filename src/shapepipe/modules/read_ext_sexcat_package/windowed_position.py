"""WINDOWED POSITION.

SExtractor's windowed centroid (``XWIN_IMAGE``, ``YWIN_IMAGE``,
``FLAGS_WIN``), measured on a tile image at given starting positions, and the
background map SExtractor subtracts before measuring it.

The centroid is a numpy port of ``compute_winpos`` in SExtractor 2.25.0
(``src/winpos.c``), the version the tile SExtractor run uses; the
background is an approximation of its ``makeback`` (``src/back.c``). Both
work on 0-based pixel indices internally and take and return 1-based (FITS)
positions, as SExtractor does.

:Author: Cail Daley

"""

import numpy as np
from scipy.interpolate import CubicSpline

# SExtractor 2.25.0, src/winpos.h and src/growth.h.
WINPOS_NITERMAX = 16  # maximum number of steps
WINPOS_NSIG = 4  # measurement radius, in window sigmas
WINPOS_OVERSAMP = 11  # subpixels per pixel side on the aperture boundary
WINPOS_STEPMIN = 1e-4  # stop once the mean offset is below this, in pixels
WINPOS_FAC = 2.0  # step = WINPOS_FAC * mean offset (exact for a Gaussian)
GROWTH_MINHLRAD = 0.5  # minimum half-light radius, in pixels

# FLAGS_WIN bits. 1, 2, 4 and 8 are SExtractor's own, each of which makes it
# revert to the isophotal barycentre; 16 and 32 are ShapePipe's.
FLAGS_WIN_SINGULAR = 1  # windowed second moments are singular
FLAGS_WIN_NEGATIVE_MOMENT = 2  # a windowed second moment is negative
FLAGS_WIN_NEGATIVE_FLUX = 4  # windowed flux is not positive
FLAGS_WIN_WANDERED = 8  # centroid left the isophotal 1-sigma ellipse
FLAGS_WIN_UNMEASURED = 16  # no measurement: barycentre not finite or off image
FLAGS_WIN_NOT_CONVERGED = 32  # still moving after 16 steps; position kept
FLAGS_WIN_FALLBACK = (
    FLAGS_WIN_SINGULAR
    | FLAGS_WIN_NEGATIVE_MOMENT
    | FLAGS_WIN_NEGATIVE_FLUX
    | FLAGS_WIN_WANDERED
    | FLAGS_WIN_UNMEASURED
)

# Pixels per object and chunk, bounding the working arrays (~100 MB).
_CHUNK_PIXELS = 2_000_000


def _natural_spline(nodes, at, axis):
    """Natural cubic spline through equally spaced nodes 0, 1, ..., n - 1.

    Outside the nodes the end cubic continues, as in SExtractor's
    ``subbackline``; a single node is constant.
    """
    if nodes.shape[axis] == 1:
        return np.repeat(nodes, len(at), axis=axis)
    grid = np.arange(nodes.shape[axis])
    return CubicSpline(grid, nodes, axis=axis, bc_type="natural")(at)


def _mesh_mode(pix):
    """SExtractor's background estimate for the pixels of one mesh.

    ``backstat`` clips at the mean +- 2 sigma, and ``backguess`` iterates a
    3 sigma clip about the median, then takes the mode estimate
    2.5 median - 1.5 mean unless the distribution is too skewed, when it
    keeps the median. SExtractor does both on a quantised histogram; this
    does them on the pixel values.
    """
    mean, sigma = pix.mean(), pix.std()
    pix = pix[np.abs(pix - mean) <= 2 * sigma]
    lo, hi = -np.inf, np.inf
    for _ in range(100):
        sel = pix[(pix >= lo) & (pix <= hi)]
        med, mean, sigma = np.median(sel), sel.mean(), sel.std()
        new_lo, new_hi = med - 3 * sigma, med + 3 * sigma
        if new_lo == lo and new_hi == hi:
            break
        lo, hi = new_lo, new_hi
    if sigma > 0 and abs(mean - med) / sigma < 0.3:
        return 2.5 * med - 1.5 * mean
    return med


def sextractor_background(image, back_size=512, filter_size=9, good=None):
    """Background map as SExtractor's ``BACK_TYPE AUTO`` estimates it.

    The image is cut into ``back_size`` meshes (the last row and column
    smaller); each mesh with at least half its pixels good gets the clipped
    mode of :func:`_mesh_mode`, a bad mesh the mean of its nearest good
    meshes; a ``filter_size`` median filter, narrowed symmetrically at the
    edges, smooths the mesh map (``BACK_FILTTHRESH`` 0); and a bicubic
    natural spline through the mesh centres interpolates it to every pixel.

    Parameters
    ----------
    image : numpy.ndarray
        2-D image, shape ``(ny, nx)``
    back_size : int, optional
        ``BACK_SIZE``, the mesh side in pixels
    filter_size : int, optional
        ``BACK_FILTERSIZE``, the median-filter side in meshes
    good : numpy.ndarray, optional
        Boolean map of the pixels to use (SExtractor: weight above
        threshold); default all finite pixels

    Returns
    -------
    numpy.ndarray
        Background map, float32, shape ``(ny, nx)``

    @sc [decision:detection.background_model]
    """
    ny, nx = image.shape
    nbx, nby = (nx - 1) // back_size + 1, (ny - 1) // back_size + 1
    if good is None:
        good = np.isfinite(image)
    back = np.full((nby, nbx), np.nan)
    for j in range(nby):
        for i in range(nbx):
            cut = np.s_[j * back_size:(j + 1) * back_size,
                        i * back_size:(i + 1) * back_size]
            pix = image[cut][good[cut]]
            if pix.size >= 0.5 * image[cut].size:
                back[j, i] = _mesh_mode(pix.astype(np.float64))

    ok = np.argwhere(~np.isnan(back))
    for j, i in np.argwhere(np.isnan(back)):
        if len(ok) == 0:
            back[j, i] = 0.0
            continue
        d2 = ((ok - (j, i)) ** 2).sum(1)
        back[j, i] = back[tuple(ok[d2 == d2.min()].T)].mean()

    half = filter_size // 2
    filtered = np.empty_like(back)
    for j in range(nby):
        hy = min(half, j, nby - 1 - j)
        for i in range(nbx):
            hx = min(half, i, nbx - 1 - i)
            filtered[j, i] = np.median(
                back[j - hy:j + hy + 1, i - hx:i + hx + 1]
            )

    fy = np.arange(ny) / back_size - 0.5
    fx = np.arange(nx) / back_size - 0.5
    rows = _natural_spline(filtered, fy, axis=0)
    out = np.empty((ny, nx), np.float32)
    step = max(1, 4_000_000 // nx)
    for y0 in range(0, ny, step):
        out[y0:y0 + step] = _natural_spline(rows[y0:y0 + step], fx, axis=1)
    return out


def _aperture_fraction(dx, dy, raper2):
    """Fraction of each pixel's 11 x 11 subpixel grid inside the aperture.

    ``winpos.c`` counts the subpixel centres ``(dx + i/11, dy + j/11)``,
    ``i, j`` in ``-5..5``, with squared radius below ``raper2``; this counts
    them row by row in closed form.
    """
    half = WINPOS_OVERSAMP // 2
    sub = np.arange(-half, half + 1) / WINPOS_OVERSAMP
    rem = raper2[:, None] - (dy[:, None] + sub) ** 2
    s = np.sqrt(np.clip(rem, 0, None))
    lo = np.maximum(
        np.floor(WINPOS_OVERSAMP * (-s - dx[:, None])) + 1, -half
    )
    hi = np.minimum(np.ceil(WINPOS_OVERSAMP * (s - dx[:, None])) - 1, half)
    count = np.where(rem > 0, np.clip(hi - lo + 1, 0, None), 0)
    return count.sum(1) / WINPOS_OVERSAMP**2


def _blanked(label, own):
    """Whether pixels of ``label`` are blanked while measuring ``own``.

    SExtractor (``MASK_TYPE`` BLANK or CORRECT) sets every detected pixel
    to -BIG as the scan passes it, and pastes an object's own pixels back
    only while it measures that object: every footprint but the object's
    own is blanked, including those no catalogue object claims.
    """
    return (label != 0) & (label != own)


def _iterate(image, mx, my, sig, seg=None, own=None,
             n_iter_max=WINPOS_NITERMAX):
    """Run ``compute_winpos``'s loop for objects sharing one stamp size.

    ``mx``, ``my`` are 0-based starting positions and are not modified.
    With ``seg``, pixels :func:`_blanked` for the object take the value
    of their mirror image through the current centre, or 0 where the
    mirror is off the image or blanked too: SExtractor's ``MASK_TYPE
    CORRECT``. Returns the final
    0-based positions, the windowed flux ``tv``, the windowed second
    moments of the last step, and whether each object stopped on the step
    criterion.
    """
    ny, nx = image.shape
    n = len(mx)
    mx, my = mx.astype(np.float64), my.astype(np.float64)
    invtwosig2 = 1.0 / (2.0 * sig * sig)
    raper = WINPOS_NSIG * sig
    raper2 = raper * raper
    rintlim2 = np.clip(raper - 0.75, 0, None) ** 2
    rextlim2 = (raper + 0.75) ** 2
    half = int(np.ceil(raper.max())) + 2
    offs = np.arange(-half, half + 1)

    tv = np.zeros(n)
    mom = np.zeros((3, n))
    converged = np.zeros(n, bool)
    active = np.arange(n)
    for _ in range(n_iter_max):
        if active.size == 0:
            break
        amx, amy = mx[active], my[active]
        xs = np.rint(amx).astype(np.int64)[:, None] + offs
        ys = np.rint(amy).astype(np.int64)[:, None] + offs
        # The C box: [(int)(m - raper + 0.499999), (int)(m + raper
        # + 1.499999)) clipped to the image; (int) truncates, and the
        # clip makes truncation and floor agree.
        r = raper[active][:, None]
        in_x = (xs >= np.clip(np.trunc(amx[:, None] - r + 0.499999), 0, None))
        in_x &= xs < np.clip(np.trunc(amx[:, None] + r + 1.499999), None, nx)
        in_y = (ys >= np.clip(np.trunc(amy[:, None] - r + 0.499999), 0, None))
        in_y &= ys < np.clip(np.trunc(amy[:, None] + r + 1.499999), None, ny)

        dx = (xs - amx[:, None])[:, None, :]
        dy = (ys - amy[:, None])[:, :, None]
        r2 = dx * dx + dy * dy
        area = (r2 < rextlim2[active, None, None]).astype(np.float64)
        ring = (r2 > rintlim2[active, None, None]) & (area > 0)
        if ring.any():
            k, iy, ix = np.nonzero(ring)
            area[k, iy, ix] = _aperture_fraction(
                dx[k, 0, ix], dy[k, iy, 0], raper2[active][k]
            )
        area *= in_y[:, :, None] & in_x[:, None, :]
        cy = np.clip(ys, 0, ny - 1)[:, :, None]
        cx = np.clip(xs, 0, nx - 1)[:, None, :]
        pix = image[cy, cx]
        if seg is not None:
            lab = seg[cy, cx]
            blank = _blanked(lab, own[active, None, None]) & (area > 0)
            if blank.any():
                k, iy, ix = np.nonzero(blank)
                # (int)(2 m + 0.49999 - x), truncating as C does.
                x2 = np.trunc(2 * amx[k] + 0.49999 - xs[k, ix])
                y2 = np.trunc(2 * amy[k] + 0.49999 - ys[k, iy])
                x2, y2 = x2.astype(np.int64), y2.astype(np.int64)
                inside = (x2 >= 0) & (x2 < nx) & (y2 >= 0) & (y2 < ny)
                x2c, y2c = np.clip(x2, 0, nx - 1), np.clip(y2, 0, ny - 1)
                mlab = seg[y2c, x2c]
                usable = inside & ~_blanked(mlab, own[active][k])
                pix = pix.astype(np.float64)
                pix[k, iy, ix] = np.where(usable, image[y2c, x2c], 0.0)
        locpix = area * np.exp(-r2 * invtwosig2[active, None, None]) * pix

        atv = locpix.sum((1, 2))
        pos = atv > 0
        safe = np.where(pos, atv, 1.0)
        dxpos = (locpix * dx).sum((1, 2)) / safe
        dypos = (locpix * dy).sum((1, 2)) / safe
        tv[active] = atv
        mom[0, active] = (locpix * dx * dx).sum((1, 2)) / safe - dxpos**2
        mom[1, active] = (locpix * dy * dy).sum((1, 2)) / safe - dypos**2
        mom[2, active] = (locpix * dx * dy).sum((1, 2)) / safe - dxpos * dypos
        mx[active] += np.where(pos, dxpos * WINPOS_FAC, 0.0)
        my[active] += np.where(pos, dypos * WINPOS_FAC, 0.0)
        done = pos & (dxpos**2 + dypos**2 < WINPOS_STEPMIN**2)
        converged[active[done]] = True
        active = active[pos & ~done]
    return mx, my, tv, mom, converged


def windowed_positions(image, x_image, y_image, flux_radius, seg=None,
                       number=None, cxx=None, cyy=None, cxy=None):
    """SExtractor's windowed centroid, started from given positions.

    For each object, a Gaussian window of sigma
    ``2 * hl / 2.35`` (``hl`` = ``|flux_radius|``, at least 0.5 pixel)
    weights the image within ``WINPOS_NSIG`` sigmas, the aperture edge
    subsampled 11 x 11; the centre moves by ``WINPOS_FAC`` times the
    weighted mean offset, up to ``WINPOS_NITERMAX`` times, until that
    offset is below ``WINPOS_STEPMIN``. Given the segmentation map,
    neighbours' footprints take the mirror image of the object's side, as
    SExtractor's ``MASK_TYPE CORRECT`` does. Where SExtractor would revert
    to the isophotal position (``FLAGS_WIN`` 1, 2, 4, 8) this returns the
    input position, and so it does where there is nothing to measure (16).
    A centroid still moving after the last step keeps its position, as
    SExtractor's does, and is flagged 32.

    Parameters
    ----------
    image : numpy.ndarray
        Background-subtracted image, shape ``(ny, nx)``
    x_image, y_image : array_like
        1-based starting positions: SExtractor's ``X_IMAGE``, ``Y_IMAGE``
    flux_radius : array_like
        Half-light radius, SExtractor's ``FLUX_RADIUS`` for
        ``PHOT_FLUXFRAC`` 0.5, in pixels
    seg : numpy.ndarray, optional
        Segmentation map on the image grid, 0 for sky, labelled with
        ``number`` (:func:`read_ext_sexcat.relabel_seg`)
    number : array_like, optional
        Each object's label in ``seg``, required with it
    cxx, cyy, cxy : array_like, optional
        Isophotal ellipse (``CXX_IMAGE`` ...) for the flag-8 test; without
        them bit 8 is never set

    Returns
    -------
    numpy.ndarray
        ``XWIN_IMAGE``, 1-based
    numpy.ndarray
        ``YWIN_IMAGE``, 1-based
    numpy.ndarray
        ``FLAGS_WIN``, int16

    @sc [decision:preparation.object_position_columns]
    """
    ny, nx = image.shape
    x0 = np.asarray(x_image, dtype=np.float64)
    y0 = np.asarray(y_image, dtype=np.float64)
    fr = np.asarray(flux_radius, dtype=np.float64)
    n = len(x0)
    flags = np.zeros(n, np.int16)
    xwin, ywin = x0.copy(), y0.copy()

    ok = np.isfinite(x0) & np.isfinite(y0) & np.isfinite(fr)
    ok &= (x0 >= 0.5) & (x0 < nx + 0.5) & (y0 >= 0.5) & (y0 < ny + 0.5)
    flags[~ok] |= FLAGS_WIN_UNMEASURED
    sig = np.maximum(np.abs(np.where(ok, fr, 1.0)), GROWTH_MINHLRAD) * 2 / 2.35
    half = np.ceil(WINPOS_NSIG * sig).astype(np.int64) + 2

    own = None if seg is None else np.asarray(number)
    tv = np.zeros(n)
    mom = np.zeros((3, n))
    for h in np.unique(half[ok]):
        members = np.flatnonzero(ok & (half == h))
        size = (2 * h + 1) ** 2
        for c0 in range(0, len(members), max(1, _CHUNK_PIXELS // size)):
            idx = members[c0:c0 + max(1, _CHUNK_PIXELS // size)]
            mx, my, tv[idx], mom[:, idx], conv = _iterate(
                image, x0[idx] - 1.0, y0[idx] - 1.0, sig[idx], seg,
                None if own is None else own[idx],
            )
            xwin[idx], ywin[idx] = mx + 1.0, my + 1.0
            stopped_on_flux = tv[idx] <= 0
            flags[idx[~conv & ~stopped_on_flux]] |= FLAGS_WIN_NOT_CONVERGED

    mx2, my2, mxy = mom
    measured = ok
    flags[measured & (mx2 * my2 - mxy * mxy < 0)] |= FLAGS_WIN_SINGULAR
    flags[measured & ((mx2 < 0) | (my2 < 0))] |= FLAGS_WIN_NEGATIVE_MOMENT
    flags[measured & (tv <= 0)] |= FLAGS_WIN_NEGATIVE_FLUX
    if cxx is not None:
        dx, dy = xwin - x0, ywin - y0
        outside = (
            np.asarray(cxx) * dx * dx
            + np.asarray(cyy) * dy * dy
            + np.asarray(cxy) * dx * dy
        ) > 1.0
        flags[measured & (dx * dx > 1) & (dy * dy > 1) & outside] |= (
            FLAGS_WIN_WANDERED
        )
    revert = (flags & FLAGS_WIN_FALLBACK) != 0
    xwin[revert], ywin[revert] = x0[revert], y0[revert]
    return xwin, ywin, flags


def isophotal_ellipse(a_world, b_world, theta_j2000, wcs, x_image, y_image):
    """``CXX_IMAGE``, ``CYY_IMAGE``, ``CXY_IMAGE`` from the world ellipse.

    Maps the isophotal ellipse (``A_WORLD``, ``B_WORLD`` in degrees,
    ``THETA_J2000`` east of north) to pixels through the local Jacobian of
    ``wcs`` at each 1-based position, and returns its inverse-covariance
    coefficients in SExtractor's convention
    (``CXX x^2 + CYY y^2 + CXY x y = 1`` on the 1-sigma ellipse).
    """
    x = np.asarray(x_image, dtype=np.float64)
    y = np.asarray(y_image, dtype=np.float64)
    eps = 0.5
    ra, dec = wcs.all_pix2world(
        np.concatenate([x, x + eps, x]), np.concatenate([y, y, y + eps]), 1
    )
    n = len(x)
    cosd = np.cos(np.deg2rad(dec[:n]))
    # d(east, north)/d(x, y), degrees per pixel; east = +RA.
    dra = (ra[n:] - np.tile(ra[:n], 2) + 180.0) % 360.0 - 180.0
    ddec = dec[n:] - np.tile(dec[:n], 2)
    jac = np.empty((n, 2, 2))
    jac[:, 0, 0], jac[:, 0, 1] = dra[:n] * cosd / eps, dra[n:] * cosd / eps
    jac[:, 1, 0], jac[:, 1, 1] = ddec[:n] / eps, ddec[n:] / eps

    theta = np.deg2rad(np.asarray(theta_j2000, dtype=np.float64))
    u = np.stack([np.sin(theta), np.cos(theta)], -1)
    v = np.stack([np.cos(theta), -np.sin(theta)], -1)
    a2 = np.asarray(a_world, dtype=np.float64) ** 2
    b2 = np.asarray(b_world, dtype=np.float64) ** 2
    cov_w = (a2[:, None, None] * u[:, :, None] * u[:, None, :]
             + b2[:, None, None] * v[:, :, None] * v[:, None, :])
    jinv = np.linalg.inv(jac)
    cov = jinv @ cov_w @ np.swapaxes(jinv, 1, 2)
    icov = np.linalg.inv(cov)
    return icov[:, 0, 0], icov[:, 1, 1], 2 * icov[:, 0, 1]
