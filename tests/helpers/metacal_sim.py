"""Shared metacal-on-a-sim plumbing for the fast science guardrails."""

import numpy as np

from shapepipe.modules.ngmix_package.ngmix import (
    Postage_stamp,
    do_ngmix_metacal,
    get_prior,
)
from shapepipe.testing.simulate import make_data

METACAL_STEP = 0.01   # == metacal_pars['step'] (ngmix.py:1393); R denominator
PIXEL_SCALE = 0.1857


def build_stamp(sim, transform=None):
    """Fill a ``Postage_stamp`` from a ``make_data`` sim, optionally transformed.

    ``transform`` (if given) is applied to every gal/psf/weight/flag image —
    used by the geometric-invariance arms (e.g. ``np.rot90(a, 2)`` for the
    MegaCam 180 deg flip).
    """
    gals, psfs, _, weights, flags, jacobs = sim
    if transform is not None:
        gals, psfs, weights, flags = (
            [transform(a) for a in gals], [transform(a) for a in psfs],
            [transform(a) for a in weights], [transform(a) for a in flags])
    s = Postage_stamp(bkg_sub=False, megacam_flip=False)
    s.gals, s.psfs, s.weights, s.flags, s.jacobs = gals, psfs, weights, flags, jacobs
    return s


def recover(seed, shear, psf_shear=(0.0, 0.0), gal_hlr=0.3, psf_fwhm=0.55,
            transform=None, n_epochs=2, img_size=51, noise=1e-4):
    """Run metacal on one controlled sim; return recovered g and response R.

    Returns a dict with ``g`` (len-2), ``R`` (2x2), and the flattened scalars
    ``g1``/``g2``/``R11``/``R12``/``R21``/``R22`` the guardrails assert on.
    """
    rng = np.random.RandomState(seed)
    prior = get_prior(PIXEL_SCALE, rng)
    sim = make_data(rng=np.random.RandomState(seed + 100), shear=shear,
                    psf_shear=psf_shear, noise=noise, n_epochs=n_epochs,
                    img_size=img_size, gal_hlr=gal_hlr, psf_fwhm=psf_fwhm)
    res, _, _ = do_ngmix_metacal(build_stamp(sim, transform), prior, 1.0, rng)
    step = METACAL_STEP
    R = np.array([
        [(res["1p"]["g"][0] - res["1m"]["g"][0]) / (2 * step),
         (res["2p"]["g"][0] - res["2m"]["g"][0]) / (2 * step)],
        [(res["1p"]["g"][1] - res["1m"]["g"][1]) / (2 * step),
         (res["2p"]["g"][1] - res["2m"]["g"][1]) / (2 * step)]])
    g = np.array(res["noshear"]["g"])
    g = g[::-1]  # FAULT-INJECTION: g1<->g2 swap
    return dict(g=g, R=R, g1=g[0], g2=g[1],
                R11=R[0, 0], R12=R[0, 1], R21=R[1, 0], R22=R[1, 1])
