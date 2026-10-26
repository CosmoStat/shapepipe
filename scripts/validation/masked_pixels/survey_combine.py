"""Project measured prevalence through fixture responses; no catalogue refits.

Output is a one-tile, equal-object-weight transport model, NOT a measured
survey calibration. Statistical errors and model envelopes are separate.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from bias_vs_distance import GALAXIES, INTERPOLATED
from survey_prevalence import epoch_rules

BASE = Path(__file__).resolve().parent
ARMS = ('A', 'B', 'C', 'D0', 'D2', 'E0', 'E2')


def bootstrap_response(row, draws):
    groups = row['seed_groups']
    e, r, a, rm = [np.asarray([g[k] for g in groups]) for k in ('e', 'R', 'A', 'Rm')]
    c = np.linalg.solve(r[draws].mean(axis=1), e[draws].mean(axis=1)[..., None])[..., 0]
    m = np.linalg.solve(rm[draws].mean(axis=1), a[draws].mean(axis=1)) - np.eye(2)
    return np.column_stack((m[:, 0, 0], m[:, 1, 1], c))


def vector(row):
    return np.asarray(row['m'] + row['c'])


def fixture_geometry():
    """Map observed FLUX_RADIUS to observed radii of the three fixtures.

    This avoids treating convolved SExtractor FLUX_RADIUS as intrinsic hlr.
    It is still a proxy: isolated exponential profiles, no deblending/noise.
    """
    import galsim
    radii, psf_t = [], []
    fine_scale = .1857 / 4
    y, x = np.indices((601, 601)) - 300
    radius = np.hypot(y, x).ravel() * fine_scale
    order = np.argsort(radius)
    for hlr, fwhm in GALAXIES:
        psf = galsim.Moffat(beta=2.5, fwhm=fwhm)
        gal = galsim.Convolve(psf, galsim.Exponential(half_light_radius=hlr))
        image = gal.drawImage(nx=601, ny=601, scale=fine_scale).array.ravel()
        cumulative = np.cumsum(image[order])
        radii.append(float(np.interp(cumulative[-1] / 2, cumulative, radius[order])))
        mom = psf.drawImage(nx=51, ny=51, scale=.1857).FindAdaptiveMom()
        psf_t.append(float(2 * (mom.moments_sigma * .1857) ** 2))
    return np.asarray(radii), np.asarray(psf_t)


def library(extra, nboot=500):
    rows = json.load(open(BASE / 'napoli_veto_scan.json'))
    controls = json.load(open(BASE / 'napoli_clean_controls.json'))
    extra_rows = json.load(open(extra))
    for name in ('fixtures-core.json', 'fixtures-interp.json'):
        additional = Path(extra).with_name(name)
        if additional.exists():
            extra_rows += json.load(open(additional))
    clean = {(r['hlr'], tuple(r['psf_shear'])): r for r in controls}
    entries = {}
    for row in rows:
        fill = 'interpolate' if row['fill'] == 'production' and row['kind'] in INTERPOLATED else 'noise'
        r = dict(row, fill=fill, n_epochs=1)
        entries[(r['kind'], fill, r['d'], r['hlr'], tuple(r['psf_shear']))] = r
    for r in extra_rows:
        if r['n_epochs'] == 1 and r['kind'] != 'clean':
            entries[(r['kind'], r['fill'], r['d'], r['hlr'], tuple(r['psf_shear']))] = r
    rng = np.random.RandomState(91)
    draws = rng.randint(6, size=(nboot, 6))
    cboot = {k: bootstrap_response(r, draws) for k, r in clean.items()}
    lib, q, qb, invalid = [], [], [], []
    for r in entries.values():
        if 'error' in r:
            invalid.append(r)
            continue
        c = clean[(r['hlr'], tuple(r['psf_shear']))]
        lib.append(r)
        q.append(vector(r) - vector(c))
        qb.append(bootstrap_response(r, draws) - cboot[(r['hlr'], tuple(r['psf_shear']))])
    # Store N * delta(q_N) as an effective single-epoch coefficient.
    # Multiplication by an epoch's alpha then reproduces delta(q_N) when
    # alpha=1/N. For unequal weights, effective N=1/alpha is a transport
    # approximation, explicitly compared with uncorrected 1/N below.
    clean_n = {(r['hlr'], tuple(r['psf_shear']), r['n_epochs']): r
               for r in extra_rows if r['kind'] == 'clean' and 'error' not in r}
    multi = {(r['kind'], r['fill'], r['d'], r['hlr'], tuple(r['psf_shear']), r['n_epochs']): r
             for r in extra_rows if r['n_epochs'] > 1 and r['kind'] != 'clean'
             and 'error' not in r}
    for r in multi.values():
        ne = r['n_epochs']
        control = clean_n[(r['hlr'], tuple(r['psf_shear']), ne)]
        lib.append(r)
        q.append(ne * (vector(r) - vector(control)))
        qb.append(ne * (bootstrap_response(r, draws) - bootstrap_response(control, draws)))
    return lib, np.asarray(q), np.asarray(qb), invalid, extra_rows


def nearest_components(components):
    """At most one component per fill per epoch; retain all for envelopes."""
    order = np.lexsort((components['distance'], components['fill'], components['ei']))
    c = components[order]
    keep = np.ones(len(c), dtype=bool)
    keep[1:] = ((c['ei'][1:] != c['ei'][:-1]) | (c['fill'][1:] != c['fill'][:-1]))
    return c[keep]


def component_lookup(c, e, o, lib, q, radii, psf_t):
    # Coupled galaxy/PSF grid: compare observed HLR and sqrt(PSF T), in arcsec.
    obs_r = o['FLUX_RADIUS'][e['oi'][c['ei']]] * .1857
    obs_r = np.where(np.isfinite(obs_r) & (obs_r > 0), obs_r, radii[0])
    t = e['psf_T'][c['ei']]
    t = np.where(np.isfinite(t) & (t > 0), t, psf_t[0])
    cost = ((obs_r[:, None] - radii) ** 2
            + (np.sqrt(t)[:, None] - np.sqrt(psf_t)) ** 2)
    sizes = cost.argmin(axis=1)
    g = np.hypot(e['psf_g1'][c['ei']], e['psf_g2'][c['ei']])
    elliptical = g > np.hypot(.05, .02) / 2
    indices = np.empty(len(c), dtype=int)
    unsupported_shape = np.zeros(len(c), dtype=bool)
    for i, part in enumerate(c):
        kind = str(part['kind'])
        if kind == 'blob' or (part['fill'] == 'interpolate' and kind not in INTERPOLATED):
            kind = 'finite_bleed' if part['fill'] == 'interpolate' else 'cluster4'
            unsupported_shape[i] = True
        h = GALAXIES[sizes[i]][0]
        ell = (.05, .02) if elliptical[i] else (0., 0.)
        candidates = [j for j, r in enumerate(lib) if r['fill'] == part['fill']
                      and r['kind'] == kind and r['hlr'] == h and tuple(r['psf_shear']) == ell
                      and r.get('n_epochs', 1) == 1]
        indices[i] = min(candidates, key=lambda j: abs(lib[j]['d'] - min(float(part['distance']), 20.)))
    # Empirical scenario envelope: all tested sizes, PSFs and geometries at
    # neighbouring distance bins, summed over ALL components (no cancellation).
    envelopes = np.empty((len(c), 4))
    for fill in ('interpolate', 'noise'):
        for distance in range(21):
            target = (c['fill'] == fill) & (np.minimum(np.rint(c['distance']), 20) == distance)
            if not target.any():
                continue
            pool = [j for j, r in enumerate(lib) if r['fill'] == fill]
            nearest = min(abs(lib[j]['d'] - distance) for j in pool)
            pool = [j for j in pool if abs(lib[j]['d'] - distance) <= max(1., nearest)]
            envelopes[target] = np.max(np.abs(q[pool]), axis=0)
    return indices, envelopes, dict(unsupported_shape=int(unsupported_shape.sum()),
                                   outside_observed_size_grid=int(((obs_r < radii.min()) | (obs_r > radii.max())).sum()),
                                   beyond_20px=int((c['distance'] > 20).sum()))


def matched_size_envelope(c, ci, lib, q):
    """Keep measured size/footprint correlation; vary geometry, PSF and N.

    This still assumes the nearest convolved-size bin is appropriate; the
    unrestricted envelope deliberately drops even that assumption.
    """
    cache = {}
    output = []
    for part, index in zip(c, ci):
        h = lib[index]['hlr']
        d = min(int(np.rint(part['distance'])), 20)
        key = (str(part['fill']), h, d)
        if key not in cache:
            pool = [i for i, row in enumerate(lib) if row['fill'] == key[0] and row['hlr'] == h]
            nearest = min(abs(lib[i]['d'] - d) for i in pool)
            pool = [i for i in pool if abs(lib[i]['d'] - d) <= max(1., nearest)]
            cache[key] = np.max(np.abs(q[pool]), axis=0)
        output.append(cache[key])
    return np.asarray(output)


def build_keep(e, nobjects):
    masks = {k: np.asarray([epoch_rules(r)[k] for r in e]) & e['eligible']
             for k in ARMS if k != 'C'}
    rejected = np.zeros(nobjects, dtype=bool)
    flag = e['eligible'] & ((e['iso_bits'] != 0) | e['interp_hit0'] | ~e['known'])
    rejected[e['oi'][flag]] = True
    masks['C'] = e['eligible'] & ~rejected[e['oi']]
    return masks


def combine(e, o, c, ci, q, qb, envelopes, masks, sample, lib, ivar=True, use_multi=True):
    n = len(o)
    multi_index = {}
    for j, r in enumerate(lib):
        if r.get('n_epochs', 1) > 1:
            key = (r['kind'], r['fill'], r['d'], r['hlr'], tuple(r['psf_shear']))
            multi_index.setdefault(key, []).append(j)
    weight = e['ivar'].copy() if ivar else np.ones(len(e))
    weight = np.where(np.isfinite(weight) & (weight > 0), weight, 0.)
    blocks = (np.clip(((o['X_IMAGE'] - o['X_IMAGE'].min()) /
                      np.ptp(o['X_IMAGE']) * 6).astype(int), 0, 5)
              + 6 * np.clip(((o['Y_IMAGE'] - o['Y_IMAGE'].min()) /
                             np.ptp(o['Y_IMAGE']) * 6).astype(int), 0, 5))
    block_draws = np.random.default_rng(774).multinomial(36, np.ones(36)/36, size=500)
    baseline = masks['B'] & sample[e['oi']]
    nbase = int((np.bincount(e['oi'][baseline], minlength=n) > 0).sum())
    totalweight_b = np.bincount(e['oi'][baseline], weights=weight[baseline], minlength=n)
    kish_b = totalweight_b.sum() ** 2 / (totalweight_b ** 2).sum()
    result = {}
    for arm in ARMS:
        keep = masks[arm] & sample[e['oi']]
        wsum = np.bincount(e['oi'][keep], weights=weight[keep], minlength=n)
        alive = (wsum > 0) & sample
        alpha = np.zeros(len(e))
        alpha[keep] = weight[keep] / wsum[e['oi'][keep]]
        selected_ci = ci.copy()
        if use_multi:
            for ic in np.flatnonzero(alpha[c['ei']] > 0):
                row = lib[ci[ic]]
                key = (row['kind'], row['fill'], int(np.rint(c['distance'][ic])),
                       row['hlr'], tuple(row['psf_shear']))
                candidates = multi_index.get(key, [])
                neff = 1 / alpha[c['ei'][ic]]
                if candidates and neff >= 1.5:
                    selected_ci[ic] = min(candidates, key=lambda j: abs(lib[j]['n_epochs'] - neff))
        epoch_q = np.zeros((len(e), 4))
        np.add.at(epoch_q, c['ei'], q[selected_ci])
        oq = np.column_stack([np.bincount(e['oi'], weights=alpha * epoch_q[:, j], minlength=n)
                              for j in range(4)])
        mean = oq[alive].mean(axis=0)
        coefficient = np.bincount(selected_ci, weights=alpha[c['ei']] / alive.sum(), minlength=len(q))
        calibration = np.einsum('l,lbq->bq', coefficient, qb)
        fixture_sigma = calibration.std(axis=0)
        bn = np.bincount(blocks[alive], minlength=36)
        bs = np.array([np.bincount(blocks[alive], weights=oq[alive, j], minlength=36)
                       for j in range(4)]).T
        denominator = block_draws @ bn
        nonempty = denominator > 0
        bmean = (block_draws[nonempty] @ bs) / denominator[nonempty, None]
        block_sigma = bmean.std(axis=0)
        # Column defects correlate different square blocks along a stripe.
        # Keep both sky axes coherent and use the larger strip/block SD,
        # rather than treating 36 squares as independent by assumption.
        strip_sigmas = []
        for axis in (blocks % 6, blocks // 6):
            counts = np.bincount(axis[alive], minlength=6)
            sums = np.array([np.bincount(axis[alive], weights=oq[alive, j], minlength=6)
                             for j in range(4)]).T
            draws = np.random.default_rng(775).multinomial(6, np.ones(6)/6, size=500)
            denominator_s = draws @ counts
            valid = denominator_s > 0
            strip_sigmas.append(((draws[valid] @ sums) / denominator_s[valid, None]).std(axis=0))
        strip_sigma = np.max(strip_sigmas, axis=0)
        spatial_sigma = np.maximum(block_sigma, strip_sigma)
        # Summing absolute component predictions (not signed cancellation).
        env = (alpha[c['ei'], None] * envelopes).sum(axis=0) / alive.sum()
        result[arm] = dict(mean=mean.tolist(), fixture_sigma=fixture_sigma.tolist(),
                           spatial_sigma=spatial_sigma.tolist(),
                           spatial_block_sigma=block_sigma.tolist(),
                           spatial_strip_sigma=strip_sigma.tolist(),
                           sigma=np.hypot(fixture_sigma, spatial_sigma).tolist(),
                           fixture_envelope=env.tolist(), objects=int(alive.sum()),
                           objects_lost=nbase-int(alive.sum()),
                           object_loss_fraction=1-int(alive.sum())/nbase,
                           epochs=int(keep.sum()), epochs_lost=int(baseline.sum()-keep.sum()),
                           epoch_loss_fraction=1-keep.sum()/baseline.sum(),
                           ivar_information_loss=1-wsum.sum()/totalweight_b.sum(),
                           kish_proxy_ratio=float(wsum.sum()**2/(wsum**2).sum()/kish_b),
                           baseline_objects=nbase, baseline_epochs=int(baseline.sum()),
                           multi_corrected_components=int(((selected_ci != ci) & (alpha[c['ei']] > 0)).sum()))
        # Contribution by fill and central-distance bin, especially D's
        # retained defects OUTSIDE the footprint, not removed defects.
        parts = {}
        for fill in ('interpolate', 'noise'):
            for lo, hi in ((0, 7), (7, 13), (13, 99)):
                sel = (c['fill'] == fill) & (c['distance'] >= lo) & (c['distance'] < hi)
                contrib = (alpha[c['ei'][sel], None] * q[selected_ci[sel]]).sum(axis=0) / alive.sum()
                parts[f'{fill}_{lo}_{hi}'] = contrib.tolist()
        result[arm]['contributions'] = parts
    return result


def rates(o, e):
    n = len(o)
    eligible = e['eligible']
    baseline = np.bincount(e['oi'][eligible], minlength=n) > 0
    out = dict(all_objects=n, candidate_epochs=len(e), common_eligible_epochs=int(eligible.sum()),
               baseline_objects=int(baseline.sum()),
               epoch_N_distribution={str(i): int(v) for i, v in enumerate(np.bincount(np.bincount(e['oi'][eligible], minlength=n)))},
               invalid_flux_radius_objects=int((baseline & (~np.isfinite(o['FLUX_RADIUS']) | (o['FLUX_RADIUS'] <= 0))).sum()))
    for label, mask in (('candidate', np.ones(len(e), bool)), ('eligible', eligible)):
        near = mask & ((e['d_interp'] < 13) | (e['d_noise'] < 13))
        out[label + '_any_defect_13px_objects'] = int(len(np.unique(e['oi'][near])))
    out['exposures'] = sorted(set(x.split('-')[0] for x in e['exposure']))
    unknown = np.zeros(n, bool)
    unknown[e['oi'][eligible & ~e['known']]] = True
    positive = np.zeros(n, bool)
    positive[e['oi'][eligible & ((e['iso_bits'] != 0) | e['interp_hit0'])]] = True
    out['unknown_footprint_objects'] = int((unknown & baseline).sum())
    out['positive_iso_flag_objects'] = int((positive & baseline).sum())
    centre_missing = eligible & e['known'] & ~e['centre_in_foot']
    out['known_footprint_without_centre_objects'] = int(len(np.unique(e['oi'][centre_missing])))
    retained_centre = eligible & e['known'] & ~e['hit0'] & (e['d_interp'] == 0)
    out['D0_retained_central_interpolation_epochs'] = int(retained_centre.sum())
    out['psf_T_percentiles'] = np.percentile(e['psf_T'][eligible & (e['psf_T'] > 0)], [10, 50, 90]).tolist()
    blind = eligible & e['hit0'] & ~e['interp_hit0'] & (e['iso_bits'] == 0) & e['known']
    out['noise_iso_no_973_flag_epochs'] = int(blind.sum())
    out['noise_iso_no_973_flag_objects'] = int(len(np.unique(e['oi'][blind])))
    out['class_distance'] = []
    for fill, key in (('interpolate', 'd_interp'), ('noise', 'd_noise')):
        for lo, hi in ((0, 3), (3, 7), (7, 13), (13, 99)):
            m = eligible & (e[key] >= lo) & (e[key] < hi)
            out['class_distance'].append(dict(fill=fill, lo=lo, hi=hi, epochs=int(m.sum()),
                                             objects=int(len(np.unique(e['oi'][m])))))
    out['size_hits'] = []
    for lo, hi in ((0, 2), (2, 3), (3, 4), (4, 6), (6, 1000)):
        sample = baseline & (o['FLUX_RADIUS'] >= lo) & (o['FLUX_RADIUS'] < hi)
        known = np.zeros(n, bool)
        known[e['oi'][eligible & e['known']]] = True
        sample &= known
        r = dict(lo=lo, hi=hi, objects=int(sample.sum()))
        for k in ('hit0', 'hit2', 'interp_hit0', 'interp_hit2'):
            hit = np.zeros(n, bool)
            hit[e['oi'][eligible & e[k] & e['known']]] = True
            r[k] = int((hit & sample).sum())
        out['size_hits'].append(r)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', default='/scratch/cdaley/veto-survey')
    args = p.parse_args()
    root = Path(args.root)
    z = np.load(root / 'prevalence.npz')
    o, e, all_c = [z[k] for k in ('objects', 'epochs', 'components')]
    all_c = all_c[e['eligible'][all_c['ei']]]
    c = nearest_components(all_c)
    lib, q, qb, invalid, extra = library(root / 'fixtures.json')
    radii, psf_t = fixture_geometry()
    ci, env, diagnostics = component_lookup(c, e, o, lib, q, radii, psf_t)
    ai, aenv, all_diagnostics = component_lookup(all_c, e, o, lib, q, radii, psf_t)
    masks = build_keep(e, len(o))
    all_sample = np.ones(len(o), dtype=bool)
    # A transparent sensitivity sample, NOT the survey's actual shear cuts.
    galaxy_proxy = ((o['SNR_WIN'] >= 10) & (o['FLAGS'] == 0)
                    & (o['FLUX_RADIUS'] * .1857 >= .55)
                    & (o['FLUX_RADIUS'] * .1857 <= 2.))
    projected = combine(e, o, c, ci, q, qb, env, masks, all_sample, lib)
    single_epoch = combine(e, o, c, ci, q, qb, env, masks, all_sample, lib, use_multi=False)
    # Conservative all-component envelope, separate from nearest-component mean.
    all_components = combine(e, o, all_c, ai, q, qb, aenv, masks, all_sample, lib)
    sized_envelope = matched_size_envelope(all_c, ai, lib, q)
    matched = combine(e, o, all_c, ai, q, qb, sized_envelope, masks, all_sample, lib)
    uniform = combine(e, o, c, ci, q, qb, env, masks, all_sample, lib, ivar=False)
    proxy = combine(e, o, c, ci, q, qb, env, masks, galaxy_proxy, lib)
    near = c['distance'] <= 20
    no_tail = combine(e, o, c[near], ci[near], q, qb, env[near], masks, all_sample, lib)
    # ISO sentinel sensitivity: literal C flag!=0 removes unknowns; the
    # alternative cuts positive flags only. D's default remains fixed-veto
    # fallback on unknown footprints, with an explicitly optimistic variant.
    positive = np.zeros(len(o), bool)
    positive[e['oi'][e['eligible'] & ((e['iso_bits'] != 0) | e['interp_hit0'])]] = True
    optimistic_masks = dict(masks)
    optimistic_masks['C'] = e['eligible'] & ~positive[e['oi']]
    for arm, key in (('D0', 'hit0'), ('D2', 'hit2')):
        optimistic_masks[arm] = e['eligible'] & ~e[key]
    optimistic = combine(e, o, c, ci, q, qb, env, optimistic_masks, all_sample, lib)
    for arm in ARMS:
        projected[arm]['all_component_envelope'] = all_components[arm]['fixture_envelope']
        projected[arm]['all_component_mean'] = all_components[arm]['mean']
        projected[arm]['matched_size_envelope'] = matched[arm]['fixture_envelope']
    # Direct N=1 equivalence to the independent original helper.
    original = json.load(open(BASE / 'napoli_veto_scan.json'))
    direct = next(r for r in extra if r['kind'] == 'bleed' and r['fill'] == 'interpolate'
                  and r['hlr'] == .3 and r['d'] == 4 and r['n_epochs'] == 1)
    prior = next(r for r in original if r['kind'] == 'bleed' and r['fill'] == 'production'
                 and r['hlr'] == .3 and r['d'] == 4 and r['psf_shear'] == [0., 0.])
    assert np.allclose(vector(direct), vector(prior), atol=1e-12, rtol=0)
    result = dict(arms=projected, equal_epoch_weights=uniform, galaxy_proxy=proxy,
                  uncorrected_single_epoch=single_epoch, zero_beyond_20px=no_tail,
                  optimistic_unknown_footprints=optimistic,
                  prevalence=rates(o, e), mapping=diagnostics, all_mapping=all_diagnostics,
                  observed_fixture_hlr_arcsec=radii.tolist(), fixture_psf_T=psf_t.tolist(),
                  invalid_fixture_cases=[{k: r[k] for k in ('kind', 'fill', 'd', 'hlr', 'psf_shear')}
                                         for r in invalid],
                  library_rows=len(lib), nearest_components=len(c), all_components=len(all_c),
                  eligible_epochs_multiple_components=int((np.bincount(all_c['ei'], minlength=len(e)) > 1).sum()),
                  n1_equivalence=True)
    (root / 'combined.json').write_text(json.dumps(result, indent=2))
    np.savez_compressed(root / 'mapping.npz', components=c, library_index=ci, library_q=q)
    (root / 'library_keys.json').write_text(json.dumps([{k: r[k] for k in ('kind', 'fill', 'd', 'hlr', 'psf_shear', 'n_epochs')}
                                                      for r in lib]))
    for arm, r in projected.items():
        print(arm, 'q', r['mean'], 'sigma', r['sigma'], 'env', r['all_component_envelope'],
              'loss', r['objects_lost'], r['epochs_lost'])


if __name__ == '__main__':
    main()
