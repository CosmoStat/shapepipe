"""Read-only 188.307 object/epoch census for counterfactual veto rules.

Run ONLY in an allocation and the campaign SIF. All outputs go to --out.
SQLite connections use immutable read-only URIs, not SqliteDict writers.
"""
import argparse
import json
import pickle
import sqlite3
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation, find_objects, label

from shapepipe.modules.ngmix_package import ngmix as ngm

RADIAL = np.hypot(*(np.indices((51, 51)) - 25))
KINDS = ('pixel', 'column', 'bleed', 'finite_bleed', 'cluster4', 'edge', 'blob')
EP_DTYPE = [('oi', 'i4'), ('exposure', 'U20'), ('eligible', '?'), ('fraction', 'f4'),
            ('d_interp', 'f4'), ('d_noise', 'f4'), ('hit0', '?'), ('hit2', '?'),
            ('interp_hit0', '?'), ('interp_hit2', '?'), ('iso_bits', 'i8'),
            ('known', '?'), ('ivar', 'f8'), ('psf_T', 'f4'), ('psf_g1', 'f4'),
            ('psf_g2', 'f4')]
CP_DTYPE = [('ei', 'i4'), ('fill', 'U12'), ('kind', 'U16'), ('distance', 'f4'),
            ('height', 'i2'), ('width', 'i2'), ('npix', 'i4'), ('hit0', '?'),
            ('hit2', '?')]


def epoch_rules(r):
    """Keep masks, sharing the 10% fraction cut outside this function.

    Unknown/clipped ISO footprints use the fixed veto in D/E, never a
    fictitious clean footprint. C's object-level cut is applied later.
    """
    fixed = r['d_interp'] >= 7 and r['d_noise'] >= 13
    known = r['known']
    return dict(A=fixed, B=True,
                D0=(not r['hit0']) if known else fixed,
                D2=(not r['hit2']) if known else fixed,
                E0=(r['d_noise'] >= 13 and not r['interp_hit0']) if known else fixed,
                E2=(r['d_noise'] >= 13 and not r['interp_hit2']) if known else fixed)


def classify_epoch(flag, weight, rms, tile_vign, seg, number):
    neighbour, off_tile = ngm.split_tile_markers(tile_vign, flag.shape)
    weight = np.where(off_tile, 0, weight)
    defect = ngm.defect_mask(weight, flag, rms)
    interp = ngm.interpolated_defects(defect, neighbour, 'uberseg', 'interpolate')
    noise = defect & ~interp
    foot = seg == number
    known = bool(foot.any() and not (foot[0].any() or foot[-1].any()
                                    or foot[:, 0].any() or foot[:, -1].any()))
    yy, xx = np.indices((5, 5)) - 2
    foot2 = binary_dilation(foot, structure=xx * xx + yy * yy <= 4)
    r = dict(fraction=float(defect.mean()),
             d_interp=float(RADIAL[interp].min()) if interp.any() else np.inf,
             d_noise=float(RADIAL[noise].min()) if noise.any() else np.inf,
             hit0=bool((defect & foot).any()), hit2=bool((defect & foot2).any()),
             interp_hit0=bool((interp & foot).any()),
             interp_hit2=bool((interp & foot2).any()), known=known,
             iso_bits=int(np.bitwise_or.reduce(flag[foot].astype(np.int64), initial=0)))
    parts = []
    for fill, mask in (('interpolate', interp), ('noise', noise)):
        labs, _ = label(mask, np.ones((3, 3)))
        for k, sl in enumerate(find_objects(labs), 1):
            part = labs[sl] == k
            h, w = part.shape
            n = int(part.sum())
            if (off_tile[sl] & part).any():
                kind = 'edge'
            elif n <= 6:
                kind = 'pixel'
            elif w == 1 and h >= 8:
                kind = 'column'
            elif w <= 3 and h >= 20:
                kind = 'bleed'
            elif w <= 3 and h >= 4:
                kind = 'finite_bleed'
            elif w >= 4 and h >= 8:
                kind = 'cluster4'
            else:
                kind = 'blob'
            parts.append(dict(fill=fill, kind=kind,
                              distance=float(RADIAL[sl][part].min()),
                              height=h, width=w, npix=n,
                              hit0=bool((foot[sl] & part).any()),
                              hit2=bool((foot2[sl] & part).any())))
    return r, parts


TILE = None
PATHS = None
HEADERS = None
DB = None


def initialize():
    global DB
    DB = {key: sqlite3.connect(f'file:{PATHS[key]}?mode=ro&immutable=1', uri=True)
          for key in ('image', 'weight', 'flag', 'rms', 'psf', 'background')}


def get(key, number):
    row = DB[key].execute('SELECT value FROM unnamed WHERE key=?', (str(number),)).fetchone()
    return pickle.loads(bytes(row[0])) if row else None


def chunk(indices):
    epochs, components, skips = [], [], dict(no_psf=0, missing=0, zeros=0, mad=0)
    for oi in indices:
        row = TILE[oi]
        num = int(row['NUMBER'])
        psfs = get('psf', num)
        if not isinstance(psfs, dict):
            skips['no_psf'] += 1
            continue
        stores = {k: get(k, num) for k in ('image', 'weight', 'flag', 'rms', 'background')}
        if any(not isinstance(s, dict) for s in stores.values()):
            skips['missing'] += 1
            continue
        for name, psf in psfs.items():
            if name not in stores['image']:
                continue
            gal = stores['image'][name]['VIGNET']
            if np.all(gal == 0):
                skips['zeros'] += 1
                continue
            sub = ngm.background_subtract(gal, stores['background'][name]['VIGNET'])
            if not ngm.sigma_mad(sub) > 0:
                skips['mad'] += 1
                continue
            ccd = int(name.split('-')[1])
            tv = ngm.Ngmix.MegaCamFlip(row['VIGNET'], ccd)
            seg = ngm.Ngmix.MegaCamFlip(row['SEG_VIGNET'], ccd)
            flag = stores['flag'][name]['VIGNET']
            rms = stores['rms'][name]['VIGNET']
            r, parts = classify_epoch(flag, stores['weight'][name]['VIGNET'], rms, tv, seg, num)
            scale = HEADERS[name]
            good = np.isfinite(rms) & (rms > 0) & (flag == 0)
            ivar = float(np.median(1 / (rms[good] * scale) ** 2)) if good.any() else 0.
            shape = psf['SHAPES']
            r.update(oi=int(oi), exposure=str(name), eligible=r['fraction'] <= .10,
                     ivar=ivar, psf_T=float(shape.get('HSM_T_PSF', np.nan)),
                     psf_g1=float(shape.get('HSM_G1_PSF', np.nan)),
                     psf_g2=float(shape.get('HSM_G2_PSF', np.nan)))
            ei = len(epochs)
            epochs.append(tuple(r[k] for k, _ in EP_DTYPE))
            components.extend(tuple(dict(ei=ei, **p)[k] for k, _ in CP_DTYPE) for p in parts)
    return np.array(epochs, dtype=EP_DTYPE), np.array(components, dtype=CP_DTYPE), skips


def main():
    global TILE, PATHS, HEADERS
    import fitsio
    from astropy.io import fits
    from multiprocessing import get_context
    from sqlitedict import SqliteDict
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='/scratch/cdaley/veto-survey')
    parser.add_argument('--nproc', type=int, default=12)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    metadata = json.load(open('/scratch/cdaley/defectflag-check/188.307.json'))
    PATHS = metadata['inputs']
    cols = ['NUMBER', 'FLUX_RADIUS', 'SNR_WIN', 'MAG_AUTO', 'FLAGS', 'N_EPOCH',
            'X_IMAGE', 'Y_IMAGE', 'VIGNET', 'SEG_VIGNET']
    TILE = fitsio.read(PATHS['sexcat'], ext='LDAC_OBJECTS', columns=cols)
    HEADERS = {}
    with SqliteDict(PATHS['headers'], flag='r') as d:
        for exposure in d:
            if exposure == 'TILE_ID':
                continue
            for ccd, item in enumerate(d[exposure]):
                try:
                    HEADERS[f'{exposure}-{ccd}'] = float(fits.Header.fromstring(item['header'])['FSCALE'])
                except (KeyError, TypeError):
                    continue
    chunks = [np.arange(i, min(i + 256, len(TILE))) for i in range(0, len(TILE), 256)]
    eparts, cparts, skipped = [], [], {}
    offset = 0
    with get_context('fork').Pool(args.nproc, initializer=initialize) as pool:
        for i, (e, c, skip) in enumerate(pool.imap(chunk, chunks)):
            c['ei'] += offset
            offset += len(e)
            eparts.append(e)
            cparts.append(c)
            for k, v in skip.items():
                skipped[k] = skipped.get(k, 0) + v
            if i % 10 == 0:
                print(f'objects {min((i+1)*256,len(TILE))}/{len(TILE)} epochs {offset}', flush=True)
    epochs, components = np.concatenate(eparts), np.concatenate(cparts)
    scalar_cols = [k for k in cols if k not in ('VIGNET', 'SEG_VIGNET')]
    objects = np.empty(len(TILE), dtype=[(k, TILE.dtype[k]) for k in scalar_cols])
    for k in scalar_cols:
        objects[k] = TILE[k]
    # Match the existing independently generated fixed-veto diagnostics.
    prior = np.load('/scratch/cdaley/defectflag-check/188.307.npz')
    assert np.array_equal(objects['NUMBER'], prior['number'])
    fixed = epochs['eligible'] & (epochs['d_interp'] >= 7) & (epochs['d_noise'] >= 13)
    nfixed = np.bincount(epochs['oi'][fixed], minlength=len(objects))
    considered = np.bincount(epochs['oi'], minlength=len(objects))
    check = dict(objects=len(objects), considered=len(epochs),
                 masked_fraction=int((~epochs['eligible']).sum()),
                 central_veto=int((epochs['eligible'] & ~fixed).sum()),
                 fixed_kept=int(fixed.sum()),
                 fixed_cache_mismatches=int((nfixed != prior['n_kept']).sum()),
                 candidate_cache_mismatches=int((considered != prior['n_candidate']).sum()),
                 skipped=skipped, inputs=PATHS)
    print(json.dumps(check, indent=2), flush=True)
    assert check['fixed_cache_mismatches'] == 0
    np.savez_compressed(out / 'prevalence.npz', objects=objects, epochs=epochs, components=components)
    (out / 'prevalence_checks.json').write_text(json.dumps(check, indent=2))


if __name__ == '__main__':
    main()
