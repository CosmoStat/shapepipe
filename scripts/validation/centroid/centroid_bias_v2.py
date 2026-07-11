#!/usr/bin/env python
"""centroid_bias_v2.py

Metacal multiplicative-bias validation using the ShapePipe ngmix v2.0 interface.

Simulation from shapepipe.testing.simulate; fitting via
shapepipe.modules.ngmix_package.ngmix (Postage_stamp, do_ngmix_metacal, get_prior).

Usage
-----
    python centroid_bias_v2.py [--ntrial N] [--seed S] [--noise SIGMA]
"""

import argparse

import numpy as np

from shapepipe.modules.ngmix_package.ngmix import (
    Postage_stamp,
    do_ngmix_metacal,
    get_prior,
)
from shapepipe.testing.simulate import make_data


# ---------------------------------------------------------------------------
# Analysis utilities
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
        ("s2n", "f8"), ("g", "f8", 2), ("T", "f8"),
    ]
    data = np.zeros(1, dtype=dt)
    data["shear_type"] = shear_type
    data["flags"] = res["flags"]
    if res["flags"] == 0:
        data["s2n"] = res.get("s2n", np.nan)
        data["g"] = res["g"]
        data["T"] = res.get("T", np.nan)
    else:
        data["s2n"] = data["g"] = data["T"] = np.nan
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

    dg = 0.02

    rng = np.random.RandomState(seed)
    prior = get_prior(pixel_scale, rng)
    shear_types = ["noshear", "1p", "1m"]

    seeds = rng.randint(0, 2 ** 30, size=ntrial)

    dlist_p, dlist_m = [], []

    for i in progress(ntrial, miniters=10):
        d_ = []
        for shear_true in [(dg, 0.00), (-dg, 0.00)]:
            gals, psfs, psfs_sigmas, weights, flags, jacob_lists = make_data(
                rng=np.random.RandomState(seeds[i]),
                noise=sig_noise,
                shear=shear_true,
                n_epochs=n_epochs,
                share_shift=False,
            )

            stamp = Postage_stamp(bkg_sub=False, megacam_flip=False)
            stamp.gals = gals
            stamp.psfs = psfs
            stamp.weights = weights
            stamp.flags = flags
            stamp.jacobs = jacob_lists

            try:
                resdict, psf_res = do_ngmix_metacal(stamp, prior, 1.0, rng)
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
        (data_p["g"][select(data_p, "1p"), 0] - data_p["g"][select(data_p, "1m"), 0]) / dg
    ).T
    w_m = select(data_m, "noshear")
    R11_m = np.atleast_2d(
        (data_m["g"][select(data_m, "1p"), 0] - data_m["g"][select(data_m, "1m"), 0]) / dg
    ).T

    shear_ = (data_p["g"][w_p] - data_m["g"][w_m]) / (R11_p + R11_m)
    shear = np.mean(shear_, axis=0)
    shear_err = np.std(shear_, axis=0) / np.sqrt(len(shear_))

    m = shear[0] / dg - 1
    merr = shear_err[0] / dg
    s2n = data_p["s2n"][w_p].mean()

    print("S/N: %g" % s2n)
    print("R11: %g %g" % (np.mean(R11_p), np.mean(R11_m)))
    print("m[1e-3, 3sigmas]: %g +/- %g (99.7%% conf)" % (m / 1e-3, merr * 3 / 1e-3))
    print("c[1e-5, 3sigmas]: %g +/- %g (99.7%% conf)" % (shear[1] / 1e-5, shear_err[1] * 3 / 1e-5))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Metacal centroid bias validation (v2 interface)"
    )
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
