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
    return dict(g=g, R=R, g1=g[0], g2=g[1],
                R11=R[0, 0], R12=R[0, 1], R21=R[1, 0], R22=R[1, 1])


# --------------------------------------------------------------------------- #
# Paired few-seed m estimator — the statistically-powered ladder rung.
# --------------------------------------------------------------------------- #
BOOTSTRAP_B = 1000       # bootstrap-over-seeds replicates for the m1 error.
DEGENERACY_K = 5         # |Rbar| < K*se(Rbar) => response consistent with zero.


def _paired_arm(seed, shear, gal_hlr, psf_fwhm, n_epochs, img_size, noise):
    """One ±γ arm of a paired seed: recovered g1 and response R11.

    The pairing rests on a fresh ``RandomState(seed)`` for the fit and
    ``RandomState(seed + 1000)`` into ``make_data``, so the two arms of a seed
    (``shear = (+γ,0)`` and ``(−γ,0)``) see a BIT-IDENTICAL noise realisation
    and sub-pixel shift — galsim ``.shear()`` consumes no rng — and their
    difference cancels the shared shape noise. Same seed offset and construction
    as ``run_breakdown_grid.py``'s ``one_stamp`` so the two agree by design.
    Returns ``(g1, R11, ok)`` with ``ok`` false on a non-finite / flagged fit.
    """
    rng = np.random.RandomState(seed)
    prior = get_prior(PIXEL_SCALE, rng)
    sim = make_data(rng=np.random.RandomState(seed + 1000), shear=shear,
                    noise=noise, n_epochs=n_epochs, img_size=img_size,
                    gal_hlr=gal_hlr, psf_fwhm=psf_fwhm)
    res, _, _ = do_ngmix_metacal(build_stamp(sim), prior, 1.0, rng)
    ns = res["noshear"]
    finite = bool(np.all(np.isfinite(ns["g"])))
    R11 = ((res["1p"]["g"][0] - res["1m"]["g"][0]) / (2 * METACAL_STEP)
           if finite else np.nan)
    ok = finite and int(ns["flags"]) == 0
    return (ns["g"][0] if finite else np.nan), R11, ok


def paired_m1(gamma, gal_hlr, psf_fwhm, n_seeds, seed0, n_epochs=2,
              img_size=51, noise=1e-4, response_scale=1.0):
    """Paired few-seed multiplicative bias ``m1`` at one ladder rung.

    Ports ``run_breakdown_grid.py``'s v2 estimator onto the ladder: drives a
    ``g1p``/``g1m`` arm per seed with a bit-identical noise realisation (see
    :func:`_paired_arm`), then

        ``m1 = (⟨g1[g1p]⟩ − ⟨g1[g1m]⟩) / (2 γ Rbar11) − 1``,

    with ``Rbar11`` the mean R11 over both arms. The ±γ difference cancels the
    per-seed shape noise that made the old single-seed estimator's σ_m ≈ the
    tolerance itself. Error is a bootstrap over the seed *pairs* (one seed
    resample per replicate, both arms resampled together so the pairing
    survives — resampling arms independently would silently inflate σ_m to the
    unpaired level). ``response_scale`` multiplies Rbar (fault-injection hook:
    a wrong metacal response shifts every rung's m by a near-constant offset).

    Returns a dict with ``m``, ``m_err``, ``R11`` (= scaled Rbar), ``R11_err``,
    ``g1_recovered`` (mean over both arms, back at +γ sign), ``n_pairs``, and
    ``degenerate`` (True when ``|Rbar| < DEGENERACY_K·se(Rbar)`` — m/m_err then
    set null rather than dividing by a response consistent with zero).
    """
    seeds = [seed0 + s for s in range(n_seeds)]
    p = np.array([_paired_arm(s, (gamma, 0.0), gal_hlr, psf_fwhm, n_epochs,
                              img_size, noise) for s in seeds])
    m_arm = np.array([_paired_arm(s, (-gamma, 0.0), gal_hlr, psf_fwhm, n_epochs,
                                  img_size, noise) for s in seeds])
    mask = (p[:, 2] > 0) & (m_arm[:, 2] > 0)          # pair-drop policy
    gp, gm = p[mask, 0], m_arm[mask, 0]
    Rp, Rm = p[mask, 1], m_arm[mask, 1]
    npair = int(mask.sum())

    Rbar = float(0.5 * (Rp.mean() + Rm.mean())) * response_scale
    boot = np.random.default_rng(seed0)               # deterministic per rung
    idx = boot.integers(0, npair, size=(BOOTSTRAP_B, npair))
    bR = 0.5 * (Rp[idx].mean(axis=1) + Rm[idx].mean(axis=1)) * response_scale
    seR = float(np.std(bR))
    degenerate = abs(Rbar) < DEGENERACY_K * seR
    bm = (gp[idx].mean(axis=1) - gm[idx].mean(axis=1)) / (2 * gamma * bR) - 1.0
    m = None if degenerate else float(
        (gp.mean() - gm.mean()) / (2 * gamma * Rbar) - 1.0)
    m_err = None if degenerate else float(np.std(bm))
    return dict(m=m, m_err=m_err, R11=Rbar, R11_err=seR,
                g1_recovered=float(0.5 * (gp.mean() - gm.mean())),
                n_pairs=npair, degenerate=bool(degenerate))
