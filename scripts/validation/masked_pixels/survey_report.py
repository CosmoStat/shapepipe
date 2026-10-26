"""Self-contained HTML report for measured prevalence + modeled survey bias."""
import argparse
import base64
import csv
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from survey_combine import ARMS, bootstrap_response, vector

NAMES = dict(A='A · fixed 7/13', B='B · keep all', C='C · ISO object cut',
             D0='D · k=0', D2='D · k=2', E0='E · k=0 hybrid', E2='E · k=2 hybrid')
COLORS = ['#64748b', '#b45309', '#9333ea', '#0284c7', '#059669', '#475569', '#0f766e']


def figure(fig, path):
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return '<img alt="' + path.stem.replace('_', ' ') + '" src="data:image/png;base64,' + base64.b64encode(path.read_bytes()).decode() + '">'


def table(headers, rows):
    return '<div class="scroll"><table><thead><tr>' + ''.join(f'<th>{s}</th>' for s in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join(f'<td>{v}</td>' for v in r) + '</tr>' for r in rows) + '</tbody></table></div>'


def pm(value, error, scale=1):
    return f'{value*scale:+.3f} ± {error*scale:.3f}'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', default='/scratch/cdaley/veto-survey')
    args = p.parse_args()
    root = Path(args.root)
    j = json.load(open(root / 'combined.json'))
    arms, prevalence = j['arms'], j['prevalence']
    extra = []
    for name in ('fixtures.json', 'fixtures-core.json', 'fixtures-interp.json'):
        extra += json.load(open(root / name))
    stability = json.load(open(root / 'fixtures-stability.json'))
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
    images = {}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for axis, selection, bound, title in ((axes[0], slice(0,2), .01, 'Multiplicative bias'),
                                         (axes[1], slice(2,4), .0005, 'Additive bias')):
        for i, arm in enumerate(ARMS):
            r = arms[arm]
            v, s = np.abs(r['mean'][selection]), np.asarray(r['sigma'][selection])
            worse = int(np.argmax(v))
            axis.scatter(i-.12, v[worse]/bound, color=COLORS[i], marker='o', zorder=3)
            axis.scatter(i, (v[worse]+s[worse])/bound, color=COLORS[i], marker='+', zorder=3)
            axis.scatter(i+.08, max(r['matched_size_envelope'][selection])/bound,
                         color=COLORS[i], marker='D', s=20, zorder=3)
            axis.scatter(i+.22, max(r['all_component_envelope'][selection])/bound,
                         facecolors='none', edgecolors=COLORS[i], marker='D', zorder=3)
        axis.axhline(1, color='#b91c1c', ls='--', lw=1, label='requirement')
        axis.set_yscale('log')
        axis.set_ylabel('maximum absolute axis / requirement')
        axis.set_title(title + f' (bound {bound:g})')
        axis.set_xticks(range(len(ARMS)), ARMS)
        axis.grid(axis='y', alpha=.2)
    fig.suptitle('One-tile catalogue projection: dots nominal; + nominal + 1σ\n'
                 'Diamonds: filled = size-matched envelope; open = unrestricted-size stress')
    fig.tight_layout()
    images['budget'] = figure(fig, root / 'bias_budget.png')

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    bins = ['0–3', '3–7', '7–13', '≥13']
    for ax, field, denom, title in ((axes[0], 'epochs', prevalence['common_eligible_epochs'], 'Eligible object–epochs'),
                                     (axes[1], 'objects', prevalence['baseline_objects'], 'Objects with ≥1 such epoch')):
        for k, fill in enumerate(('interpolate', 'noise')):
            rows = [r for r in prevalence['class_distance'] if r['fill']==fill]
            ax.bar(np.arange(4)+(k-.5)*.35, [100*r[field]/denom for r in rows], .35,
                   label=fill, color=('#0284c7','#b45309')[k])
        ax.set_xticks(range(4), bins)
        ax.set_xlabel('nearest defect pixel [px]')
        ax.set_ylabel('% of denominator')
        ax.set_title(title)
        ax.grid(axis='y', alpha=.2)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    images['prevalence'] = figure(fig, root / 'prevalence_distance.png')

    fig, ax = plt.subplots(figsize=(9, 4))
    rows = prevalence['size_hits']
    for k, key in enumerate(('hit0','hit2')):
        ax.bar(np.arange(5)+(k-.5)*.35, [100*r[key]/r['objects'] for r in rows], .35,
               color=('#0284c7','#059669')[k], label=('ISO', 'ISO dilated 2 px')[k])
    ax.set_xticks(range(5), ['<2','2–3','3–4','4–6','≥6'])
    ax.set_xlabel('SExtractor FLUX_RADIUS [px], observed (not intrinsic HLR)')
    ax.set_ylabel('objects with ≥1 eligible footprint hit [%]')
    ax.set_title('Known footprints only; 0.1857 arcsec / pixel')
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=.2)
    fig.tight_layout()
    images['size'] = figure(fig, root / 'footprint_size.png')

    fig, axes = plt.subplots(1, 2, figsize=(12,4))
    lookup = {(r['kind'],r['fill'],r['d'],r['hlr'],tuple(r['psf_shear']),r['n_epochs']):r
              for r in extra if 'error' not in r}
    draws = np.random.RandomState(91).randint(6,size=(500,6))
    dilution_rows = []
    for ax, fill, h, field, label in ((axes[0],'interpolate',.3,2,'N × Δc1 [10⁻³]'),
                                      (axes[1],'noise',.7,0,'N × Δm11 [%]')):
        distances = (4,5,6) if fill=='interpolate' else (6,8,10,12)
        scale = 1000 if field==2 else 100
        for distance in distances:
            ns, vals, errs = [], [], []
            for ne in (1,2,3,5,7):
                r = lookup[('bleed',fill,distance,h,(0.,0.),ne)]
                c = lookup[('clean','interpolate',0,h,(0.,0.),ne)]
                delta = vector(r)-vector(c)
                sigma = (bootstrap_response(r,draws)-bootstrap_response(c,draws)).std(axis=0)
                ns.append(ne); vals.append(ne*delta[field]*scale);errs.append(ne*sigma[field]*scale)
                dilution_rows.append([fill,h,distance,ne,*delta.tolist(),*sigma.tolist()])
            ax.errorbar(ns,vals,yerr=errs,marker='o',ms=3,lw=1,label=f'{distance} px')
            ax.axhline(vals[0],color=ax.lines[-1].get_color(),ls=':',alpha=.5)
        ax.set_xlabel('N epochs, exactly one defective')
        ax.set_ylabel(label)
        ax.set_xticks([1,2,3,5,7])
        ax.set_title(('0.3″ galaxy / 0.7″ PSF' if fill=='interpolate' else '0.7″ galaxy / 0.9″ PSF')+'; round PSF')
        ax.grid(alpha=.2)
        ax.legend(frameon=False,ncol=2)
    fig.suptitle('A constant curve would mean exact 1/N dilution; dotted lines are N=1 predictions')
    fig.tight_layout()
    images['dilution'] = figure(fig, root / 'multi_epoch_dilution.png')

    fig, ax = plt.subplots(figsize=(10,3.8))
    keys = list(arms['D0']['contributions'])
    for k, arm in enumerate(('D0','D2')):
        ax.bar(np.arange(6)+(k-.5)*.35, [1e6*arms[arm]['contributions'][key][2] for key in keys],.35,
               label=NAMES[arm], color=('#0284c7','#059669')[k])
    ax.set_xticks(range(6),['interp\n0–7','interp\n7–13','interp\n≥13','noise\n0–7','noise\n7–13','noise\n≥13'])
    ax.set_ylabel('contribution to catalogue Δc1 [10⁻⁶]')
    ax.set_title('Residual bias of adaptive rules comes from retained defects outside the footprint')
    ax.axhline(0,color='#64748b',lw=.8)
    ax.legend(frameon=False)
    fig.tight_layout()
    images['outside'] = figure(fig, root / 'outside_footprint.png')

    primary_rows = []
    csv_rows = []
    for a in ARMS:
        r = arms[a]
        primary_rows.append([NAMES[a],*[pm(r['mean'][i],r['sigma'][i],1e5) for i in range(4)],
                             f"{r['objects_lost']:,} ({100*r['object_loss_fraction']:.3f}%)",
                             f"{r['epochs_lost']:,} ({100*r['epoch_loss_fraction']:.3f}%)",
                             f"{100*r['ivar_information_loss']:.3f}%"])
        csv_rows.append(dict(arm=a, **{k+suffix:r[field][i] for i,k in enumerate(('m11','m22','c1','c2'))
                                      for suffix,field in (('', 'mean'),('_sigma','sigma'),('_envelope','all_component_envelope'),
                                                          ('_matched_size_envelope','matched_size_envelope'))},
                             **{k:r[k] for k in ('objects','objects_lost','object_loss_fraction','epochs','epochs_lost',
                                                'epoch_loss_fraction','ivar_information_loss','kish_proxy_ratio')}))
    with open(root/'arms.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(csv_rows[0]));w.writeheader();w.writerows(csv_rows)
    with open(root/'dilution.csv','w') as f:
        w=csv.writer(f);w.writerow(['fill','hlr','distance','N','delta_m11','delta_m22','delta_c1','delta_c2',
                                  'sigma_m11','sigma_m22','sigma_c1','sigma_c2']);w.writerows(dilution_rows)
    rate_rows = [[r['fill'],f"[{r['lo']},{r['hi']})" if r['hi']<99 else '≥13',f"{r['epochs']:,}",
                  f"{100*r['epochs']/prevalence['common_eligible_epochs']:.3f}%",f"{r['objects']:,}",
                  f"{100*r['objects']/prevalence['baseline_objects']:.3f}%"] for r in prevalence['class_distance']]
    size_rows = [[f"{r['lo']}–{r['hi']}" if r['hi']<1000 else '≥6',r['objects'],
                  f"{100*r['hit0']/r['objects']:.2f}%",f"{100*r['hit2']/r['objects']:.2f}%"] for r in prevalence['size_hits']]
    env_rows = [[NAMES[a],f"{100*max(arms[a]['matched_size_envelope'][:2]):.3f}%",
                 f"{1e6*max(arms[a]['matched_size_envelope'][2:]):.1f}",
                 f"{100*max(arms[a]['all_component_envelope'][:2]):.3f}%",
                 f"{1e6*max(arms[a]['all_component_envelope'][2:]):.1f}",
                 f"{100*(arms[a]['kish_proxy_ratio']-1):+.3f}%"] for a in ARMS]
    sensitivity_rows=[]
    for name in ('equal_epoch_weights','galaxy_proxy','uncorrected_single_epoch','zero_beyond_20px','optimistic_unknown_footprints'):
        for a in ('B','D0','D2'):
            r=j[name][a]
            sensitivity_rows.append([name,a,*[f'{x*1e5:+.3f}' for x in r['mean']],
                                     f"{100*r['object_loss_fraction']:.3f}%",f"{100*r['epoch_loss_fraction']:.3f}%"])
    for a in ('B', 'D0', 'D2'):
        r = arms[a]
        sensitivity_rows.append(['all_components', a, *[f'{x*1e5:+.3f}' for x in r['all_component_mean']],
                                 f"{100*r['object_loss_fraction']:.3f}%", f"{100*r['epoch_loss_fraction']:.3f}%"])
    checks_file = root / 'checks_summary.json'
    checks = json.load(open(checks_file)) if checks_file.exists() else {'status': 'See adjacent check logs.'}
    checks_html = table(['Check', 'Observed result'], [[k, v] for k, v in checks.items()])
    # Stability follow-up is independent seeds 6..29, never a reason to
    # delete the original zero-flag metacal-response outlier.
    unstable = lookup[('bleed','interpolate',4,.7,(0.,0.),2)]
    repeat_bad=next(r for r in stability if r['kind']=='bleed')
    repeat_clean=next(r for r in stability if r['kind']=='clean')
    repeat_delta=vector(repeat_bad)-vector(repeat_clean)
    orig_delta=vector(unstable)-vector(lookup[('clean','interpolate',0,.7,(0.,0.),2)])
    d0central=arms['D0']['contributions']['interpolate_0_7'][2]
    d2central=arms['D2']['contributions']['interpolate_0_7'][2]
    html = f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Replacement rules for the ShapePipe central defect veto</title>
<style>
:root{{color-scheme:light;font-family:system-ui,sans-serif;color:#172334;background:#f6f8fa}}
body{{max-width:1200px;margin:0 auto;padding:32px 24px 80px;line-height:1.6}} h1,h2,h3{{line-height:1.2;color:#12243a}}
h1{{font-size:2.2rem;max-width:950px}} h2{{margin-top:2.5rem;border-top:1px solid #cbd5e1;padding-top:1.2rem}}
p{{max-width:1050px}} .lede{{font-size:1.15rem}} .box{{background:#e9f3f5;border-left:5px solid #0284c7;padding:14px 22px;margin:24px 0}}
.warn{{background:#fff5e6;border-color:#b45309}} table{{border-collapse:collapse;width:100%;font-size:.88rem;background:white}}
th,td{{padding:9px 10px;border-bottom:1px solid #dbe3ed;text-align:right}} th:first-child,td:first-child{{text-align:left}}
th{{background:#e8eef4;white-space:nowrap}} .scroll{{overflow:auto}} img{{display:block;width:100%;height:auto;margin:20px 0;background:white}}
code{{background:#e7edf2;padding:2px 5px;font-size:.88em}} pre{{overflow:auto;background:#162536;color:#e2edf7;padding:20px}}
.small{{font-size:.88rem;color:#475569}} summary{{cursor:pointer;font-weight:650}} a{{color:#0369a1}}
</style>
<header><p class="small">ShapePipe · 188.307 · conditional survey-bias projection · 25 October 2026</p>
<h1>The fixed veto has plausible adaptive replacements; a single-epoch counterexample does not decide the catalogue mean.</h1>
<p class="lede">Under the stated budget and nearest-fixture transport, I recommend <strong>D(k=0) for the next validation run</strong>, not a claim of survey-certified calibration.
It retains 1,015 more epochs than the fixed veto, with nominal defect-induced bias far below both bounds.
Its size-matched fixture envelope is also inside the budget by factors of about 25 in m and 5 in c.
D(k=2) buys additional margin at the cost of another 731 epochs.</p></header>
<div class="box"><strong>Measured:</strong> object–epoch masks, fill classification, footprint intersections, relative inverse-variance proxy and rule losses on seven rebuilt exposures.
<strong>Simulated:</strong> conditional m/c, including N=3,5,7 joint fits with one defective epoch.
<strong>Assumed:</strong> transport to real galaxy profiles, nearest fixture/PSF, relative epoch information, equal object weights and additive superposition of multiple defects.
No detector-column orientation averaging is used.</div>
<h2>Catalogue projection and measured loss</h2>
<p>All four bias columns below are in <strong>10⁻⁵</strong>; ± is the quadrature of a paired-six-seed fixture bootstrap and the largest spatial bootstrap SD from 6×6 blocks or six row/column strips of this tile.
These are <strong>defect-induced</strong> biases relative to the matching clean fixture, not a calibration of every other survey systematic.
The comparison bounds are |m11|,|m22| &lt; 0.01 (=1000 in this table) and |c1|,|c2| &lt; 5×10⁻⁴ (=50).
These are bootstrap standard deviations, not guaranteed Gaussian coverage intervals.
The statistical errors do not include population/geometry transport uncertainty; see the separate envelope and sensitivity tests.</p>
{table(['Arm','m11 ×10⁵','m22 ×10⁵','c1 ×10⁵','c2 ×10⁵','Objects lost','Epochs lost','Ivar-information lost'],primary_rows)}
<p>Losses are relative to B's <strong>35,694 objects / 115,299 eligible epochs</strong>.
All arms retain the existing 10% masked-fraction cut: it removes 5,122 of 120,421 candidate epochs and leaves 357 of 36,051 detections with no eligible epoch before any replacement rule.
“Objects lost” means additional all-epoch loss, or C's explicit object cut.
No new galaxy fits were run on the real catalogue; these are mask-based retention counts, not measured fit-success or science-sample counts.</p>
{images['budget']}
<h3>Definitions and unknown footprints</h3>
<ul><li><strong>A:</strong> strict &lt;7 px interpolation / &lt;13 px noise epoch veto.</li>
<li><strong>B:</strong> no central veto, retain all otherwise eligible epochs.</li>
<li><strong>C:</strong> no central veto; remove the object if any eligible epoch has instrument ISO bits or interpolated ISO pixels.
A literal nonzero cut also removes the −1 unknown-footprint sentinel.
There are {prevalence['positive_iso_flag_objects']:,} positive-flag objects and {prevalence['unknown_footprint_objects']:,} unknown-footprint objects; their union loses 2,042.
A positive-only variant is reported below.
Unflagged noise-filled ISO pixels are not captured by these two #973 flags, but none occurred in known eligible footprints here.</li>
<li><strong>D:</strong> reject an epoch if any production defect intersects the object's segmentation footprint, dilated by an Euclidean disk of radius k=0 or 2 px.</li>
<li><strong>E:</strong> retain the 13 px noise veto; apply D only to interpolation.</li></ul>
<p>For missing/edge-clipped footprints D/E fall back to A, rather than treating missing information as a clean footprint.
The optimistic unknown-footprint variant is shown separately.
The coadd segmentation is MegaCam-flipped into each epoch with no reprojection, matching the existing #973 overlay and its limitations.</p>
<h2>Prevalence and source size</h2>
<p><strong>{prevalence['candidate_any_defect_13px_objects']:,}/36,051 = {100*prevalence['candidate_any_defect_13px_objects']/36051:.2f}%</strong> of detections have a candidate defect epoch within 13 px.
After the common fraction cut, this is <strong>{prevalence['eligible_any_defect_13px_objects']:,}/35,694 = {100*prevalence['eligible_any_defect_13px_objects']/35694:.2f}%</strong>.
Only <strong>48/115,299 = 0.0416%</strong> of eligible object–epochs have noise-filled defects within 13 px; interpolation has 5,760/115,299 = 5.00% within 13 px.
The two classes and object distance bins can overlap; do not sum object fractions as disjoint populations.</p>
{table(['Class','Distance [px]','Epochs','% eligible epochs','Objects','% baseline objects'],rate_rows)}
{images['prevalence']}
{table(['FLUX_RADIUS [px]','Known objects','ISO hit','ISO +2 px hit'],size_rows)}
{images['size']}
<p>Eligible epoch counts per object are N=1:307, N=2:5,302, N=3:16,903, N=4:12,231, N=5:951.
The seven exposures are {', '.join(prevalence['exposures'])}; the private rebuild does not cover all nine exposures available for this tile.
This is one tile's detection population, not the final galaxy selection or a survey-wide prevalence estimate.</p>
<h2>Multi-epoch dilution is measured, not assumed to be exact</h2>
<p>The new joint fits use N=1,2,3,4,5,7 with exactly one defective epoch, the same galaxy in every epoch, equal RMS/PSF, independent subpixel shifts and six paired seeds.
Interpolation at 4–6 px is measured on the 0.3″/0.7″ grid and across the other grid points/elliptical PSFs; noise-filled bleeds at 6,8,10,12 px include the decisive 0.7″/0.9″ grid.
Additional core-noise cases resolve single-epoch ill-conditioning that would otherwise give B a misleadingly unstable 1/N error bar.</p>
{images['dilution']}
<p>For example, the 0.3″/0.7″ interpolated bleed at 4 px has NΔc1≈4.45×10⁻³ at N=1 but ≈2.17×10⁻³ at N=3 and ≈1.90×10⁻³ at N=7.
The 0.7″/0.9″ noise-filled bleed at 10 px has NΔm11≈−6.50% at N=1 and ≈−6.74% at N=3.
Thus 1/N can overpredict the first case's additive bias by about a factor two, while working well for the second.</p>
<p>The projection uses measured joint-fit NΔq where a matched multi-epoch fixture exists, choosing the nearest measured effective N=1/α; α is the retained epoch's relative inverse-RMS² information proxy.
Elsewhere it uses αΔq_single, explicitly the 1/N/linear-information assumption.
Unequal-epoch and multiple-defect interactions are not directly calibrated; equal epoch weights and all-component superposition are sensitivity checks below.</p>
<div class="box warn"><strong>A numerical response outlier is retained, not hidden.</strong>
The 0.7″/0.9″ interpolated bleed at 4 px, N=2, has Δm11={orig_delta[0]:+.3f} with a 0.349 single-case bootstrap error.
One of six seeds has an unflagged metacal response element near 12 instead of 0.94; the other five have near-normal responses.
An independent 24-seed repetition (seeds 6–29) gives Δm11={repeat_delta[0]:+.5f} and Δc1={repeat_delta[2]:+.3g}.
This diagnoses a sparse, zero-fit-flag response failure, not a reproducible −66% physical dilution law.
I nevertheless retain the original case in the broad envelope; removing it without a production-quality response criterion would be unjustified.</div>
<h2>Adaptive residuals are outside the footprint</h2>
<p>Yes: D's retained known-footprint defects are outside that footprint by construction.
More specifically, the <strong>interpolated 0–7 px outside-footprint</strong> contribution to D(k=0)'s catalogue c1 is <strong>{d0central*1e6:+.2f}×10⁻⁶</strong>, versus its total {arms['D0']['mean'][2]*1e6:+.2f}×10⁻⁶.
For k=2, the corresponding contribution is {d2central*1e6:+.2f}×10⁻⁶ and the total is {arms['D2']['mean'][2]*1e6:+.2f}×10⁻⁶.
These signed totals have partial cancellation between distance/classes; the envelope below explicitly does not rely on that cancellation.</p>
{images['outside']}
<p>{prevalence['known_footprint_without_centre_objects']:,} known footprints do not contain the stamp centre.
D(k=0) retains {prevalence['D0_retained_central_interpolation_epochs']} epochs with an interpolated central pixel outside the ISO footprint.
An adaptive footprint rule is not equivalent to a smaller circular radius.
The retained bias here is primarily from near-footprint defects, not a large population of missing central labels.</p>
<h2>Model uncertainty and effective density</h2>
<p>Both envelopes maximize absolute effective bias over tested geometries, both PSF states, measured N values and neighbouring distance bins, then sum over <strong>every component</strong> with the measured epoch-information fractions.
The <strong>size-matched envelope preserves each component's nearest observed-size bin</strong>, so it retains the measured size–footprint correlation.
The unrestricted envelope also maximizes over all three sizes; it deliberately discards that correlation.
Both use no orientation or cross-case sign cancellation and retain the sparse N=2 response outlier.
For example, only three large-grid 3–5 px interpolation components survive D(k=0), versus 49 in B; none survive D(k=2).
Applying the large-grid outlier to every small-footprint object would overstate that particular correlated scenario.
These are <strong>empirical fixture scenario envelopes, not confidence intervals or mathematical bounds for arbitrary galaxies/defects</strong>.</p>
{table(['Arm','Size-matched max |m|','Size-matched max |c| [10⁻⁶]','Unrestricted max |m|','Unrestricted max |c| [10⁻⁶]','Kish proxy Δ'],env_rows)}
<p>The density column is (ΣW_i)²/ΣW_i² with W_i the object's retained epoch inverse-RMS² sum, relative to B.
It is <strong>not actual shear n_eff</strong>: the science selection, per-object shear weights, shape dispersion and noisy fitted responses are unavailable here.
In the equal-shape-weight limit the n_eff loss is simply the measured object-loss fraction; the information-loss column gives a complementary depth proxy.</p>
<details><summary>Weight, population, tail and flag-sentinel sensitivities</summary>
<p>All bias entries remain ×10⁻⁵.
The galaxy proxy requires SNR_WIN≥10, FLAGS=0, and 0.55″≤observed FLUX_RADIUS≤2″; it is a sensitivity selection, not the collaboration's actual shear cuts.
“Uncorrected single epoch” exposes the 1/N approximation before direct N fits.
“Zero beyond 20 px” replaces the nominal nearest-20-px tail with zero.
“Optimistic unknown footprints” cuts positive #973 flags only in C and trusts partial footprints in D.</p>
{table(['Scenario','Arm','m11','m22','c1','c2','Object loss','Epoch loss'],sensitivity_rows)}
</details>
<div class="box warn"><strong>Limits on transport:</strong>
{j['mapping']['outside_observed_size_grid']:,}/{j['nearest_components']:,} mapped components have observed FLUX_RADIUS outside the three convolved fixture radii {', '.join(f'{x:.3f}″' for x in j['observed_fixture_hlr_arcsec'])}.
The full detection sample contains unresolved sources; {prevalence['invalid_flux_radius_objects']} objects have invalid/nonpositive FLUX_RADIUS and use the smallest fixture.
{j['mapping']['unsupported_shape']} nearest components have blob/unsupported morphology and use a stated finite-bleed/cluster surrogate.
{j['mapping']['beyond_20px']:,} nearest components are beyond the extended 20-px support and use that last distance (with a zero-tail sensitivity).
Actual PSF orientation/profile, galaxy morphology/SNR and detection/shear-selection response are not calibrated by these fixtures.
The strip/bootstrap calculation retains coherent stripes; it does not rotate or average their orientation away.
A spatial bootstrap inside one tile cannot estimate field-to-field prevalence uncertainty.</div>
<h2>Recommendation and decision threshold</h2>
<p><strong>Pilot D(k=0), keeping A as the unchanged production default until a selected-sample multi-tile injection test.</strong>
With the nearest-fixture transport specified here, D(k=0)'s nominal mean and size-matched envelope are well inside the stated budget while it rejects 39.6% fewer epochs than A and avoids C's 5.72% object loss.
This recommendation uses catalogue dilution and the measured size–footprint relation, not a per-case acceptance threshold.
D(k=2) is a more conservative alternative if the extra profile/segmentation margin is needed: its size-matched envelope is smaller, but it saves only 284 epochs versus A and loses five more objects than A.
B also has a small nominal mean, so the earlier single-epoch scan did not demonstrate unacceptable survey bias.
Its size-matched no-cancellation additive envelope still crosses the target; the unrestricted envelope crosses both targets.
The hybrid provides only a modest envelope reduction for 30 extra epoch drops at k=0; interpolation dominates the adaptive residual and I would not retain a fixed noise radius solely on this evidence.</p>
<p>I would choose D(k=2) or the hybrid if actual shear weights, real-mask response tails or field-to-field prevalence consume D(k=0)'s margin.
I would choose B if a matched-profile test establishes that the uncut selected-catalogue mean and flags-only policy meet the allocated systematic budget across tiles.
A tighter defect budget near 10⁻⁴ in c would also favor k=2.
The decisive missing measurement is the selected catalogue's response with real masks, not another comparison of a per-object bias against a catalogue requirement.</p>
<h2>Checks</h2>
{checks_html}
<h2>Reproduction and provenance</h2>
<p>Allocation 25005956, account def-waerbeke, ≤16 CPUs, ≤2 h; all computation in
<code>/project/def-mjhudson/cdaley/containers/shapepipe-develop-240b37e4.sif</code> with apptainer/1.4.5 and single-thread BLAS/OpenMP.
All writes are under <code>/scratch/cdaley/wt-veto</code> or <code>/scratch/cdaley/veto-survey</code>.
The private input stores are read-only/immutable; the census reproduces all 36,051 cached fixed-veto counts and candidate counts exactly.
N=1 recovery matches the original helper to 10⁻¹².
No production source, other worktree, campaign params or PR was changed.</p>
<pre>sbatch scripts/validation/masked_pixels/survey_run.sbatch
# Outputs: /scratch/cdaley/veto-survey/reproduction/report.html</pre>
<p class="small">Machine-readable companions: combined.json, arms.csv, dilution.csv, prevalence.npz, prevalence_checks.json,
fixtures.json, fixtures-core.json, fixtures-interp.json and fixtures-stability.json.
The branch stores copies under scripts/validation/masked_pixels/survey_results/.
The report is self-contained: every figure is inlined as base64; it needs no network assets.</p>
</html>'''
    (root/'report.html').write_text(html)
    summary = '# Survey-level replacement-rule projection\n\n'
    summary += 'One tile, modeled defect-induced bias, coherent orientation; not a measured survey calibration.\n'
    summary += 'Recommendation: pilot D(k=0); retain production A pending selected-sample multi-tile validation.\n\n'
    summary += '| Arm | max abs(m) ×10⁴ ± σ | max abs(c) ×10⁵ ± σ | Objects lost | Epochs lost |\n|---|---:|---:|---:|---:|\n'
    for a in ARMS:
        r=arms[a];ii=int(np.argmax(np.abs(r['mean'][:2])));ic=2+int(np.argmax(np.abs(r['mean'][2:])))
        summary += f"| {a} | {abs(r['mean'][ii])*1e4:.3f} ± {r['sigma'][ii]*1e4:.3f} | {abs(r['mean'][ic])*1e5:.3f} ± {r['sigma'][ic]*1e5:.3f} | {r['objects_lost']} ({100*r['object_loss_fraction']:.3f}%) | {r['epochs_lost']} ({100*r['epoch_loss_fraction']:.3f}%) |\n"
    summary += '\nRequirements in these units: max |m| < 100, max |c| < 50.\n'
    summary += 'Errors exclude transport uncertainty; report.html includes empirical stress envelopes and sensitivity samples.\n'
    (root/'summary.md').write_text(summary)
    print(root/'report.html')


if __name__=='__main__':
    main()
