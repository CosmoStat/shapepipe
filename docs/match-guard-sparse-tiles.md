# Sparse DR6 tiles and the catalogue-pairing guard

The external-only guard rejected localized bright-star halo/spike detections on four sparse DR6 tiles, although essentially every ShapePipe detection had a DR6 partner.
The fix keeps the SExtractor-side guard unchanged and adds a bounded, empirical 50-object allowance to the external-side fractional loss budget.
It changes acceptance, not matching, measurements, segmentation labels, or membership of an already accepted catalogue.

## Guard purpose and history

[PR #933](https://github.com/CosmoStat/shapepipe/pull/933) replaced the external-catalogue converter with native tile SExtractor measurements followed by a mutual-nearest positional join to DR6 IDs.
`git log -S min_fraction` and `git log -S tolerated_unpaired` identify the two-sided guard (`7ac98844`) and the 2%-beyond-20 policy (`59832081`), subsequently merged in `e3ef2523`.
The SExtractor side detects wrong pixels: a DR5 image against its DR6 catalogue paired only about 95% of detections in the original validation.
The external side detects a detection configuration that finds substantially fewer objects than the reference catalogue.
Neither side is a shear-eligibility cut.
External-only objects retain their external photometry but acquire no ShapePipe shape.

## Four-tile diagnosis

All diagnostics ran in the campaign SIF, `shapepipe-develop-240b37e4.sif`, on Nibi allocation 25000049, with four cores and 16 GiB.
The rerun logs identify Source Extractor 2.25.0; none reports a pixel-stack overflow or truncation warning.
Production inputs and snapshots were read-only.
The failed node-local SExtractor catalogues had already been reclaimed, so detection was reproduced with the logged command, identical binary, persisted image, weight, configuration, parameter file and filter.
Only output paths changed.
The full parameter list matters: a preliminary scalar-only run changed three detections on 577.246 and is not the diagnostic reference.
The full run reproduces all four production count tuples exactly.
Matching alone was then run on scratch copies, with the old and new external allowance.

| Tile | External | Detections | Pairs | External only | At zero weight | In star/halo masks 1 OR 2 OR 4 | Largest 100-pixel-linked group |
|---|---:|---:|---:|---:|---:|---:|---:|
| 576.248 | 1968 | 1922 | 1919 | 49 | 2 | 8 | 36 |
| 577.246 | 2490 | 2431 | 2431 | 59 | 1 | 59 | 43 |
| 603.227 | 1680 | 1625 | 1625 | 55 | 0 | 55 | 50 |
| 618.220 | 2861 | 2760 | 2760 | 101 | 9 | 101 | 38 |

The usable footprints are thin edge strips; 86.6–93.8% of each full tile has zero weight.
Nevertheless, 252 of the 264 external-only objects have positive weight at the nearest centroid pixel, and their median weights are nearly the same as the paired population's.
The missing objects concentrate in visible bright-star halos and spikes, rather than throughout the usable strip.
A 100-pixel linking length is only a diagnostic summary, not a new pipeline criterion.
The mask queries use external RA/Dec and the production `dr6-2026-09-29` maps.
The 41 unmasked external-only objects on 576.248 show why even a star-mask exemption would not explain every loss.

| Tile | Median valid MAG_AUTO, paired / external only | Median valid flux S/N, paired / external only | Invalid external photometry |
|---|---:|---:|---:|
| 576.248 | 23.27 / 23.83 | 8.50 / 5.20 | 1 |
| 577.246 | 23.21 / 24.62 | 9.05 / 3.85 | 15 |
| 603.227 | 23.17 / 25.17 | 9.23 / 2.29 | 14 |
| 618.220 | 22.94 / 24.28 | 9.17 / 4.31 | 33 |

S/N is estimated as `1.085736 / MAGERR_AUTO`; invalid photometry means nonfinite values, magnitude at least 90, or magnitude uncertainty outside `(0, 90)`.
Histograms omit those sentinels and label their counts.
There are 104 external-only objects in `20 < MAG_AUTO < 24.5`; magnitude alone does not establish that they are usable shear galaxies.
These are DR6 detections associated with halo/spike structure, not a demonstrated population of missing real galaxies.

Paired-position p99 separations are 6.56–6.69e-5 pixels, with only 3, 0, 0 and 0 unmatched ShapePipe detections.
A different stack or a WCS mismatch is therefore strongly disfavoured.
The detection command uses the committed 1-sigma threshold, three-pixel minimum area, MAP_WEIGHT treatment and zero-weight interpolation.
There is no evidence here for a global threshold or weight-treatment regression.
Local deblending/build differences are consistent with the morphology and with the unseeded SExtractor deblender discussed in #933; this diagnosis does not identify the precise binary-level cause.

## Retained-log population

The recursive census read 4,925 `process-*.log` files under all 218 existing `smk-b[0-9]*` batch directories and found 366 matching summaries.
Combining that scan with the initial matching-log scan gives 408 unique tile summaries: 400 passing and eight failing.
There is one summary per tile, not repeated attempts.
The union retains 42 summaries whose logs disappeared between scans while active production cleanup continued; newer logs also appeared during collection.
This is a collection window, not an atomic filesystem snapshot.
The other four failures are 416.238 (99/4018 external only), 048.252 (134/6619), 035.238 (67/3221) and 036.238 (72/3126), all with zero dropped SExtractor rows.
The new external budget also accepts those count tuples; their images were not part of this four-tile diagnosis.
For the 333 passing tiles with more than 30,000 reference objects, the median external-only count is 128 and the median fraction is 0.341%; the maximum fraction is 0.958%.
For the 65 tiles with 10,000–30,000 objects, the median count is 113 and the median fraction is 0.426%; the maximum fraction is 1.250%.
Only two passing summaries have fewer than 3,000 objects, with 2 and 6 external-only rows.

The four requested failures are fraction outliers relative to the retained passing sample, but their absolute losses of 49–101 objects are ordinary compared with dense tiles.
A localized halo can dominate a small footprint's denominator.
They are not established Poisson fluctuations about the dense-tile mean: the losses are spatially correlated, and the passing sample is censored by the old guard.
Successful tile cleanup deletes the module logs (`workflow/scripts/clean_tile.py`), so retained log summaries are not a complete census of all published tiles.
No false-rejection probability can be inferred from this sample.
The final recursive census and its scope are recorded beside the figures.

## Selected policy

The SExtractor side still fails when its unpaired count exceeds both 20 and `0.02 * n_sextractor`.
The external side now fails when its unpaired count exceeds both 20 and `50 + 0.02 * n_external`.
The additive allowance requires at least one pair; an empty detection catalogue does not gain it.
`MATCH_EXTERNAL_UNPAIRED_ALLOWANCE = 0` restores the old external guard.
The public function parameter follows `w_log`, preserving existing positional callers.
Older INI files receive the default through the runner's fallback.

The 50-object allowance is a rounded empirical margin above the largest observed excess over the old fractional budget, 43.78 objects on 618.220.
The new budgets for the four tiles are 89.36, 99.80, 83.60 and 107.22 objects.
It is deliberately not described as a binomial confidence bound.
The cost is weaker detection-config-drift sensitivity externally on small tiles: at 2,000 objects, losses up to 90 rather than 40 can pass.
The wrong-pixels check on the SExtractor side is not relaxed.

Rejected alternatives:

- An independent-object binomial/Poisson bound would assign unsupported significance to correlated halo children, with too few uncensored sparse tiles to calibrate it.
- Excluding zero-weight external-only objects would leave all four failures above the old 2% boundary.
- Excluding masked or clustered objects would add spatial and mask-selection policies and could hide a localized input/configuration defect.
- Raising the shared absolute tolerance would unnecessarily weaken the SExtractor-side wrong-pixels check.

## Workflow and published products

No allocation/workflow code changes are needed for tile-failure isolation.
The Nibi profile already sets `keep-going: true`; `--in-allocation` keeps that profile setting while selecting the local executor.
A failed tile writes a failed stage log but no success manifest, so its dependants stop and a later invocation retries the missing stage.
The examined latest b493 and b499 compute DAGs contain only one unfinished tile's chain, not 25 tile chains.
At inspection, b493 and b494 each have 24 final tile catalogues for 25 currently declared tiles; b499 has 6 for 7, and b512 has 13 for 14.
Their `LAUNCH.md` records declare 7 and 14 tiles respectively; they are smaller production batches, not 25-tile batches.
The remaining batch-level merge cannot finish while one required tile fails, but existing sibling tile products survive.

After an approved integration backport, resume the four campaigns through their existing allocation launchers; do not delete outputs or restart all tiles from scratch.
`sp run --in-allocation` refreshes the code snapshot and already repairs reclaimed node-local stores on resume.
Some unfinished work may be rebuilt after those stores disappear, but published per-tile catalogues need no recomputation.
No already-published tile's matching result changes: the pair selection and output rewrite are untouched, and a byte-identity regression compares the old and new policies on an accepted catalogue, including segmentation labels.
Batch-level merged catalogues will gain rows from the newly completed tiles.
The g19 branch has not been changed; deployment is awaiting review of the diagnostics.

## Evidence and verification

Canonical diagnostics are under `/scratch/cdaley/match-guard-diag/`:

- `<tile>/positions.png`: paired and external-only positions over the weight footprint, with a science-image close-up.
- `<tile>/magnitudes.png`: magnitude and S/N histograms.
- `<tile>/masks.png`: external-only positions split by downstream star/halo masks.
- `population.png`: external-only counts and fractions versus reference count, with old and new boundaries.
- `refined-stats.json`, `spatial-summary.json`, `photometry-summary.json`, `masks.json`, `full-log-census.json`, `population-observed.json`, `population-summary-final.json`, and `batch-products.json`: numerical evidence.
- `refine.py`, `masks_census.py`, `census.py`, `summarize.py`, `finish_figures.py`, `photometry.py`, and `verify_real_match.py`: reproduction scripts.

The four synthetic production-count fixtures were seen to fail under the old guard before the implementation commit.
Bad-match tests cover shifted sparse/dense grids, DR5-like 5% mismatches on either side, external-only budget boundaries, empty detections and a disabled allowance.
Three private-source mutations restore the old external guard, disable the external guard and disable the SExtractor guard; all are killed by the intended tests.
Inside the SIF, `tests/module tests/unit` passed (1036 passed, four skipped), and `tests/workflow` passed (39 tests).
A final focused run including the byte-identity regression and ASTRA integrity checks passed all 92 tests.
The three mutations produce respectively four, three and one intended failures.
The real four-tile scratch matching check also passes with exactly the production pair counts; restoring the old guard rejects each file without modifying it.
An independent read-only review found no blocking defects.
The numerical diagnostics are not production outputs and must not be copied into a running campaign.
