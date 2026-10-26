"""Multi-epoch defect recovery and outer-distance support for the survey model.

Uses the original full-matrix estimator, with exactly one defective epoch.
No quarter-turn orientation averaging of the defect itself is performed.
"""
import argparse
import json
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np

from bias_vs_distance import GALAXIES, INTERPOLATED, KINDS, geometry


def summarize(samples):
    e, r, a, rm = [np.array([s[i] for s in samples]) for i in range(4)]
    sv = float(np.linalg.svd(rm.mean(axis=0), compute_uv=False).min())
    if sv <= .1:
        raise ValueError(f'minimum response singular value {sv} <= 0.1')
    c = np.linalg.solve(r.mean(axis=0), e.mean(axis=0))
    matrix = np.linalg.solve(rm.mean(axis=0), a.mean(axis=0)) - np.eye(2)
    draw = np.random.RandomState(91).randint(len(e), size=(1000, len(e)))
    cb = np.linalg.solve(r[draw].mean(axis=1), e[draw].mean(axis=1)[..., None])[..., 0]
    mb = np.linalg.solve(rm[draw].mean(axis=1), a[draw].mean(axis=1)) - np.eye(2)
    return dict(c=c.tolist(), m=np.diag(matrix).tolist(), matrix=matrix.tolist(),
                c_err=cb.std(axis=0).tolist(), m_err=np.diag(mb.std(axis=0)).tolist(),
                min_response=sv,
                seed_groups=[dict(e=s[0].tolist(), R=s[1].tolist(), A=s[2].tolist(),
                                  Rm=s[3].tolist()) for s in samples])


def response(bad, hlr=.3, psf=.7, psf_shear=(0., 0.), n_epochs=1,
             fill='interpolate', seeds=range(6)):
    from shapepipe.modules.ngmix_package import ngmix as ngm
    from shapepipe.testing.simulate import make_data
    from tests.helpers.metacal_sim import build_stamp
    samples = []
    for seed in seeds:
        def arm(axis, sign, rotation=0):
            shear = [0., 0.]
            if axis >= 0:
                shear[axis] = sign * .02
            ps = tuple(v * (-1 if rotation else 1) for v in psf_shear)
            data = list(make_data(rng=np.random.RandomState(seed + 100), shear=shear,
                                 psf_shear=ps, noise=1e-4, n_epochs=n_epochs,
                                 img_size=51, gal_hlr=hlr, psf_fwhm=psf,
                                 return_centers=True))
            centres = data.pop()
            offsets = [np.array([c.y - 26, c.x - 26]) for c in centres]
            if rotation:
                data[0] = [np.rot90(a).copy() for a in data[0]]
                data[1] = [np.rot90(a).copy() for a in data[1]]
                offsets = [np.array([-d[1], d[0]]) for d in offsets]
            data[0][0] = np.where(bad, 1e3, data[0][0])
            data[4] = [bad.astype(np.int32)] + [np.zeros_like(bad, dtype=np.int32)
                                               for _ in range(n_epochs - 1)]
            stamp = build_stamp(data)
            stamp.offsets = offsets
            stamp.bkg_rms = [np.full((51, 51), 1e-4) for _ in range(n_epochs)]
            rng = np.random.RandomState(seed)
            res, _, _ = ngm.do_ngmix_metacal(stamp, ngm.get_prior(.1857, rng), 1., rng,
                                            centroid_source='wcs', defect_fill=fill)
            assert all(res[t]['flags'] == 0 for t in ngm.METACAL_TYPES), {
                t: res[t]['flags'] for t in ngm.METACAL_TYPES}
            e = np.asarray(res['noshear']['g'])
            R = np.column_stack([(np.asarray(res[p]['g']) - res[m]['g']) / .02
                                 for p, m in (('1p', '1m'), ('2p', '2m'))])
            return e, R
        null = [arm(-1, 0, k) for k in (0, 1)]
        e0 = np.mean([r[0] for r in null], axis=0)
        r0 = np.mean([r[1] for r in null], axis=0)
        derivatives, responses = [], []
        for axis in (0, 1):
            p, m = arm(axis, 1), arm(axis, -1)
            derivatives.append((p[0] - m[0]) / .04)
            responses.extend([p[1], m[1]])
        samples.append((e0, r0, np.column_stack(derivatives), np.mean(responses, axis=0)))
    return summarize(samples)


def one(job):
    kind, fill, distance, hlr, psf, psf_shear, ne = job
    row = dict(kind=kind, fill=fill, d=distance, hlr=hlr, psf=psf,
               psf_shear=psf_shear, n_epochs=ne)
    start = time.monotonic()
    bad = np.zeros((51, 51), dtype=bool) if kind == 'clean' else geometry(kind, distance)
    try:
        row.update(response(bad, hlr, psf, psf_shear, ne, fill))
    except Exception as exc:
        import traceback
        row.update(error=repr(exc), traceback=traceback.format_exc())
    row['seconds'] = time.monotonic() - start
    return row


def jobs():
    out = []
    for hlr, psf in GALAXIES:
        for e in ((0., 0.), (.05, .02)):
            for kind in INTERPOLATED:
                out += [(kind, 'interpolate', d, hlr, psf, e, 1)
                        for d in (8, 10, 13, 16, 20)]
            for kind in KINDS:
                out += [(kind, 'noise', d, hlr, psf, e, 1) for d in (14, 16, 20)]
    for ne in (1, 2, 3, 5, 7):
        for hlr, psf in GALAXIES:
            out.append(('clean', 'interpolate', 0, hlr, psf, (0., 0.), ne))
        for d in (4, 5, 6):
            out.append(('bleed', 'interpolate', d, .3, .7, (0., 0.), ne))
        for hlr, psf in ((.3, .7), (.7, .9)):
            for d in (6, 8, 10, 12):
                out.append(('bleed', 'noise', d, hlr, psf, (0., 0.), ne))
        for kind, d in (('cluster4', 0), ('cluster4', 2), ('bleed', 0), ('bleed', 2)):
            out.append((kind, 'noise', d, .3, .7, (0., 0.), ne))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', default='/scratch/cdaley/veto-survey/fixtures.json')
    p.add_argument('--nproc', type=int, default=14)
    p.add_argument('--jobs', help='Optional explicit JSON job list for targeted follow-up')
    args = p.parse_args()
    rows, work = [], json.load(open(args.jobs)) if args.jobs else jobs()
    with open(args.out + '.jsonl', 'w', buffering=1) as journal:
        with get_context('spawn').Pool(args.nproc) as pool:
            for r in pool.imap_unordered(one, work, chunksize=1):
                rows.append(r)
                journal.write(json.dumps(r) + '\n')
                print(f'{len(rows)}/{len(work)} {r["kind"]} {r["fill"]} '
                      f'd={r["d"]} N={r["n_epochs"]} {r.get("error", "ok")}', flush=True)
    Path(args.out).write_text(json.dumps(rows))


if __name__ == '__main__':
    main()
