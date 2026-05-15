#!/usr/bin/env python
"""centroid_bias_v2.py

Metacal multiplicative-bias validation for ShapePipe ngmix_v2.0 interface.

Uses the updated ngmix v2.4+ API:
  - ngmix.fitting.GalsimFitter  (replaces GalsimSimple)
  - priors require rng argument
  - gal centroid via inline galsim HSM  (replaces get_guess)
  - get_all_metacal without cheatnoise/symmetrize_psf

Simulation data is imported from shapepipe.testing.simulate (stable,
shared with centroid_bias.py).

Usage
-----
    python centroid_bias_v2.py [--ntrial N] [--seed S] [--noise SIGMA]
"""

import argparse

import galsim
import numpy as np
import ngmix
from ngmix import Observation, ObsList
from ngmix.fitting.galsim_results import GalsimFitModel
from ngmix.gexceptions import GMixRangeError
from numpy.random import uniform as urand

from shapepipe.modules.ngmix_package.ngmix import get_noise
from shapepipe.testing.simulate import make_data
from modopt.math.stats import sigma_mad


# ---------------------------------------------------------------------------
# Jacobian-shift fix for new ngmix GalsimFitModel
#
# GalsimFitModel.make_model applies shift = pars[0:2] only.  When the ngmix
# Jacobian center is recentered to the galaxy centroid (row0, col0 ≠ canonical),
# get_galsim_wcs() drops the center information, so the k-space
# InterpolatedImage still places the galaxy at its pixel offset from the
# canonical center.  The model shift must include this Jacobian center offset
# (in sky units) so that pars[0:2] = 0 at the solution and the centroid prior
# does not bias the shape.  This replicates Axel Guinot's _set_jacobian_shift
# fix from the aguinot/ngmix fork.
# ---------------------------------------------------------------------------

class FixedGalsimFitModel(GalsimFitModel):
    """GalsimFitModel with Jacobian-shift fix applied in make_model.

    The parent's _fill_models calls make_model(band_pars) without band/ep,
    so we override _fill_models to pass them explicitly.
    """

    def _set_jacobian_shift(self, obs_in):
        """Compute sky offset of Jacobian center from canonical image center."""

        def get_shift(obs):
            nrow, ncol = obs.image.shape
            canonical = (np.array((ncol, nrow)) - 1.0) / 2.0
            jrow, jcol = obs.jacobian.get_cen()
            offset = np.array((jcol, jrow)) - canonical
            offset *= obs.jacobian.get_scale()
            return offset

        if isinstance(obs_in, Observation):
            self._jac_shift = [[get_shift(obs_in)]]
        elif isinstance(obs_in, ObsList):
            self._jac_shift = [[get_shift(o) for o in obs_in]]
        else:  # MultiBandObsList
            self._jac_shift = [
                [get_shift(o) for o in obs_list] for obs_list in obs_in
            ]

    def __init__(self, obs, model, guess, prior=None):
        # Must set _jac_shift before super().__init__() because the parent
        # calls make_model() (via _check_guess) during initialisation.
        self._set_jacobian_shift(obs)
        super().__init__(obs, model, guess, prior=prior)

    def make_model(self, pars, band=0, ep=0):
        model = self.make_round_model(pars)
        shift = pars[0:2] + self._jac_shift[band][ep]
        try:
            model = model.shear(g1=pars[2], g2=pars[3])
        except ValueError as err:
            raise GMixRangeError(str(err))
        model = model.shift(shift)
        return model

    def _fill_models(self, pars):
        """Override to pass band and ep to make_model."""
        try:
            for band, kobs_list in enumerate(self.mb_kobs):
                band_pars = self.get_band_pars(pars, band)
                for ep, kobs in enumerate(kobs_list):
                    gal = self.make_model(band_pars, band=band, ep=ep)
                    meta = kobs.meta
                    kmodel = meta["kmodel"]
                    gal._drawKImage(kmodel)
                    if kobs.has_psf():
                        kmodel *= kobs.psf.kimage
        except RuntimeError as err:
            raise GMixRangeError(str(err))


class FixedGalsimFitter(ngmix.fitting.GalsimFitter):
    """GalsimFitter that uses FixedGalsimFitModel (Jacobian-shift fix)."""

    def _make_fit_model(self, obs, guess):
        return FixedGalsimFitModel(
            obs=obs, model=self.model, guess=guess, prior=self.prior,
        )


# ---------------------------------------------------------------------------
# Priors  (ngmix v2.4+ API: all priors require rng)
# ---------------------------------------------------------------------------

def get_prior(pixel_scale, rng=None):
    """Build ngmix joint prior for 6-parameter galaxy model."""
    if rng is None:
        rng = np.random.default_rng()
    g_prior = ngmix.priors.GPriorBA(0.4, rng=rng)
    cen_prior = ngmix.priors.CenPrior(0.0, 0.0, pixel_scale, pixel_scale, rng=rng)
    T_prior = ngmix.priors.FlatPrior(-10.0, 1.0e6, rng=rng)
    F_prior = ngmix.priors.FlatPrior(-1.0e4, 1.0e9, rng=rng)
    return ngmix.joint_prior.PriorSimpleSep(cen_prior, g_prior, T_prior, F_prior)


# ---------------------------------------------------------------------------
# Fitting  (ngmix v2.4+ API: GalsimFitter, obs passed to go())
# ---------------------------------------------------------------------------

def make_galsimfit(obs, model, guess0, prior=None, ntry=5):
    """Fit an observation with GalsimFitter, retrying with perturbed guesses.

    Parameters
    ----------
    obs : ngmix.Observation or ObsList
    model : str
    guess0 : numpy.ndarray  shape (6,)
    prior : ngmix prior, optional
    ntry : int, optional

    Returns
    -------
    dict-like
        Fit result with flags, g, T, T_err, pars, pars_err, etc.

    Raises
    ------
    ngmix.gexceptions.BootGalFailure
    """
    limit = 0.1
    guess = np.copy(guess0)
    fres = {"flags": 1}

    for it in range(ntry):
        guess[0:5] += urand(low=-limit, high=limit)
        guess[5:] *= 1 + urand(low=-limit, high=limit)
        try:
            fitter = FixedGalsimFitter(model, prior=prior)
            fres = fitter.go(obs, guess)
            if fres["flags"] == 0 and "T" not in fres:
                fres["T"] = fres["pars"][4]
                fres["T_err"] = np.sqrt(fres["pars_cov"][4, 4])
        except Exception:
            continue
        if fres["flags"] == 0:
            break

    if fres["flags"] != 0:
        raise ngmix.gexceptions.BootGalFailure(
            "Failed to fit galaxy image with GalsimFitter"
        )
    fres["ntry"] = it + 1
    return fres


# ---------------------------------------------------------------------------
# Metacal pipeline  (ngmix v2.4+ API)
# ---------------------------------------------------------------------------

def do_ngmix_metacal(
    gals, psfs, psfs_sigma, weights, flags, jacob_list,
    prior, pixel_scale, sig_noise=None,
):
    """Run metacalibration on a multi-epoch object.

    Parameters
    ----------
    gals, psfs, psfs_sigma, weights, flags : lists
        Per-epoch image arrays and metadata.
    jacob_list : list of galsim Jacobians
    prior : ngmix joint prior
    pixel_scale : float  [arcsec]
    sig_noise : float or None

    Returns
    -------
    dict
        Metacal result dict keyed by shear type ('noshear', '1p', '1m').
    """
    n_epoch = len(gals)
    if n_epoch == 0:
        raise ValueError("0 epochs to process")

    gal_obs_list = ObsList()
    T_guess_psf = []
    psf_res_gT = {
        "g_PSFo": np.array([0.0, 0.0]),
        "g_err_PSFo": np.array([0.0, 0.0]),
        "T_PSFo": 0.0,
        "T_err_PSFo": 0.0,
    }
    gal_guess = []
    gal_guess_flag = True
    wsum = 0

    for n_e in range(n_epoch):
        psf_jacob = ngmix.Jacobian(
            row=(psfs[0].shape[0] - 1) / 2,
            col=(psfs[0].shape[1] - 1) / 2,
            wcs=jacob_list[n_e],
        )
        psf_obs = Observation(psfs[n_e], jacobian=psf_jacob)
        psf_T = psfs_sigma[n_e] * 1.17741 * pixel_scale

        weight_map = np.copy(weights[n_e])
        weight_map[flags[n_e] != 0] = 0.0
        weight_map[weight_map != 0] = 1

        psf_guess = np.array([0.0, 0.0, 0.0, 0.0, psf_T, 1.0])
        try:
            psf_res = make_galsimfit(psf_obs, "gauss", psf_guess)
        except Exception:
            continue

        # Gal centroid + size guess via adaptive moments (inline HSM)
        try:
            _gim = galsim.Image(gals[n_e], scale=pixel_scale)
            _hsm = galsim.hsm.FindAdaptiveMom(_gim, strict=False)
            if _hsm.error_message != "":
                raise galsim.hsm.GalSimHSMError(_hsm.error_message)
            _cen = _hsm.moments_centroid - _gim.center
            _size = 2 * (_hsm.moments_sigma * pixel_scale) ** 2
            gal_guess_tmp = np.array([_cen.x, _cen.y, 0.0, 0.0, _size, _hsm.moments_amp])
        except Exception:
            gal_guess_flag = False
            gal_guess_tmp = np.array([0.0, 0.0, 0.0, 0.0, 1, 100])

        gal_jacob = ngmix.Jacobian(
            row=(gals[n_e].shape[0] - 1) / 2 + gal_guess_tmp[1],
            col=(gals[n_e].shape[1] - 1) / 2 + gal_guess_tmp[0],
            wcs=jacob_list[n_e],
        )

        if sig_noise is None:
            sig_noise = (
                get_noise(gals[n_e], weight_map, gal_guess_tmp, pixel_scale)
                if gal_guess_flag
                else sigma_mad(gals[n_e])
            )

        noise_img = np.random.randn(*gals[n_e].shape) * sig_noise
        noise_img_gal = np.random.randn(*gals[n_e].shape) * sig_noise

        gal_masked = np.copy(gals[n_e])
        if (weight_map == 0).any():
            gal_masked[weight_map == 0] = noise_img_gal[weight_map == 0]
        weight_map *= 1 / sig_noise ** 2

        w_tmp = np.sum(weight_map)
        psf_res_gT["g_PSFo"] += psf_res["g"] * w_tmp
        psf_res_gT["g_err_PSFo"] += (
            np.array([psf_res["pars_err"][2], psf_res["pars_err"][3]]) * w_tmp
        )
        psf_res_gT["T_PSFo"] += psf_res["T"] * w_tmp
        psf_res_gT["T_err_PSFo"] += psf_res["T_err"] * w_tmp
        wsum += w_tmp

        gal_obs = Observation(
            gal_masked, weight=weight_map,
            jacobian=gal_jacob, psf=psf_obs, noise=noise_img,
        )
        if gal_guess_flag:
            final_gal_guess = np.copy(gal_guess_tmp)
            final_gal_guess[:2] = 0
            gal_guess.append(final_gal_guess)

        gal_obs_list.append(gal_obs)
        T_guess_psf.append(psf_T)
        gal_guess_flag = True

    if wsum == 0:
        raise ZeroDivisionError("Sum of weights = 0")

    for key in psf_res_gT:
        psf_res_gT[key] /= wsum

    fail_get_guess = len(gal_guess) == 0
    gal_pars = [0.0, 0.0, 0.0, 0.0, 1, 100] if fail_get_guess else np.mean(gal_guess, 0)

    metacal_pars = {
        "types": ["noshear", "1p", "1m"],
        "step": 0.01,
        "psf": "gauss",
        "fixnoise": True,
        "use_noise_image": True,
    }
    Tguess = np.mean(T_guess_psf)

    obs_dict_mcal = ngmix.metacal.get_all_metacal(gal_obs_list, **metacal_pars)
    res = {"mcal_flags": 0}

    for key in sorted(obs_dict_mcal):
        fres = make_galsimfit(obs_dict_mcal[key], "gauss", gal_pars, prior=prior)
        res["mcal_flags"] |= fres["flags"]
        tres = dict(fres)

        wsum_psf = 0
        Tpsf_sum = 0
        gpsf_sum = np.zeros(2)
        for obs in obs_dict_mcal[key]:
            psf_obs_fit = getattr(obs, "psf_nopix", obs.psf)
            try:
                psf_res = make_galsimfit(
                    psf_obs_fit, "gauss",
                    np.array([0.0, 0.0, 0.0, 0.0, Tguess, 1.0]),
                )
            except Exception:
                continue
            tw = obs.weight.sum()
            wsum_psf += tw
            gpsf_sum += np.array(psf_res["g"]) * tw
            Tpsf_sum += psf_res["T"] * tw

        tres["gpsf"] = gpsf_sum / wsum_psf if wsum_psf > 0 else np.zeros(2)
        tres["Tpsf"] = Tpsf_sum / wsum_psf if wsum_psf > 0 else 0.0
        res[key] = tres

    res.update(psf_res_gT)
    res["moments_fail"] = fail_get_guess
    return res


# ---------------------------------------------------------------------------
# Analysis utilities  (identical to centroid_bias.py — ngmix-version agnostic)
# ---------------------------------------------------------------------------

def progress(total, miniters=1):
    """Minimal progress printer."""
    sl = str(len(str(total)))
    fmt = "%" + sl + "d/%" + sl + "d %3d%%"
    last_print_n = 0
    last_len = 0
    for i in range(total):
        yield i
        num = i + 1
        if i == 0 or num == total or num - last_print_n >= miniters:
            meter = fmt % (num, total, 100 * num // total)
            print("\r" + meter + " " * max(last_len - len(meter), 0),
                  end="", flush=True)
            last_len = len(meter)
            last_print_n = num
    print(flush=True)


def make_struct(res, shear_type):
    """Pack a metacal result into a structured array row."""
    dt = [
        ("flags", "i4"), ("shear_type", "U7"),
        ("s2n", "f8"), ("g", "f8", 2), ("T", "f8"), ("Tpsf", "f8"),
    ]
    data = np.zeros(1, dtype=dt)
    data["shear_type"] = shear_type
    data["flags"] = res["flags"]
    if res["flags"] == 0:
        data["s2n"] = res["flux"] / res["flux_err"]
        data["g"] = res["g"]
        data["T"] = res["T"]
        data["Tpsf"] = res["Tpsf"]
    else:
        data["s2n"] = data["g"] = data["T"] = data["Tpsf"] = np.nan
    return data


def select(data, shear_type):
    """Return indices where flags==0 and shear_type matches."""
    return np.where(
        (data["flags"] == 0) & (data["shear_type"] == shear_type)
    )[0]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(ntrial=50, seed=42, sig_noise=1e-10, n_epochs=3, pixel_scale=0.1857):
    prior = get_prior(pixel_scale, rng=np.random.default_rng(seed))
    shear_types = ["noshear", "1p", "1m"]

    rng = np.random.RandomState(seed)
    seeds = rng.randint(0, 2 ** 30, size=ntrial)

    dlist_p, dlist_m = [], []

    for i in progress(ntrial, miniters=10):
        d_ = []
        for shear_true in [(0.02, 0.00), (-0.02, 0.00)]:
            gals, psfs, psfs_sigmas, weights, flags, jacob_lists = make_data(
                rng=np.random.RandomState(seeds[i]),
                noise=sig_noise,
                shear=shear_true,
                n_epochs=n_epochs,
                share_shift=False,
            )
            try:
                resdict = do_ngmix_metacal(
                    gals, psfs, psfs_sigmas, weights, flags, jacob_lists,
                    prior=prior, pixel_scale=pixel_scale, sig_noise=sig_noise,
                )
                stt = [make_struct(resdict[st], st) for st in shear_types]
            except Exception:
                continue
            d_.append(np.hstack(stt))

        if len(d_) != 2:
            continue
        dlist_p.extend(d_[0])
        dlist_m.extend(d_[1])

    print()

    data_p = np.hstack(dlist_p)
    data_m = np.hstack(dlist_m)

    w_p = select(data_p, "noshear")
    R11_p = np.atleast_2d(
        (data_p["g"][select(data_p, "1p"), 0] - data_p["g"][select(data_p, "1m"), 0]) / 0.02
    ).T
    w_m = select(data_m, "noshear")
    R11_m = np.atleast_2d(
        (data_m["g"][select(data_m, "1p"), 0] - data_m["g"][select(data_m, "1m"), 0]) / 0.02
    ).T

    shear_ = (data_p["g"][w_p] - data_m["g"][w_m]) / (R11_p + R11_m)
    shear = np.mean(shear_, axis=0)
    shear_err = np.std(shear_, axis=0) / np.sqrt(len(shear_))

    m = shear[0] / 0.02 - 1
    merr = shear_err[0] / 0.02
    s2n = data_p["s2n"][w_p].mean()

    print("S/N: %g" % s2n)
    print("R11: %g %g" % (np.mean(R11_p), np.mean(R11_m)))
    print("m[1e-3, 3sigmas]: %g +/- %g (99.7%% conf)" % (m / 1e-3, merr * 3 / 1e-3))
    print("c[1e-5, 3sigmas]: %g +/- %g (99.7%% conf)" % (shear[1] / 1e-5, shear_err[1] * 3 / 1e-5))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Metacal centroid bias validation (v2 interface)")
    parser.add_argument("--ntrial", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--noise", type=float, default=1e-10)
    parser.add_argument("--n-epochs", type=int, default=3)
    parser.add_argument("--pixel-scale", type=float, default=0.1857)
    args = parser.parse_args()

    main(
        ntrial=args.ntrial,
        seed=args.seed,
        sig_noise=args.noise,
        n_epochs=args.n_epochs,
        pixel_scale=args.pixel_scale,
    )
