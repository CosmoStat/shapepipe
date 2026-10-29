# Tile-object defect diagnostics and D0 epoch selection

D0 rejects an epoch if any raw defect pixel touches the object's own ISO footprint, `SEG_VIGNET == NUMBER`, without dilation.
The footprint is the coadd segmentation stamp, MegaCam-flipped into the exposure's orientation, without reprojection.
The 10% masked-fraction cut runs first.
Both interpolated and noise-filled defects count, including unflagged zero-weight pixels, invalid background RMS and off-tile pixels.
Neighbour markers and clean quarter-turn copies introduced by defect weighting do not count as defects.

A footprint touching the stamp border is tested within the stamp, the same pixel domain used for the fit.
An empty or missing footprint rejects the epoch because overlap cannot be established.
A segmentation/defect shape mismatch raises an error.
No distance-to-centre threshold enters this selection.
The scientific rationale, alternatives and one-tile limitations live in `astra.yaml`, under `shape_measurement.footprint_defect_veto`.

## Diagnostic columns

The diagnostic records describe epochs handed to the fit.
They impose no additional cut.

- `NGMIX_IMAFLAGS_ISO` is the bitwise OR of exposure instrument flag pixels on the own footprint, OR'ed across kept epochs.
- `NGMIX_N_EPOCH_INTERP_ISO` counts kept epochs containing an interpolated defect on that footprint.
- Complete kept footprints have zero hits under D0.
  Both footprint diagnostics are `-1` for border-touching footprints because the full ISO flag OR is unknown, even when the in-stamp overlap test passes.
  The record API also gives `-1` for missing or empty footprints; production rejects those epochs before recording them.
  Never-fit final-catalogue rows carry `0`, so downstream consumers must check fit status as well as the sentinel.
- `NGMIX_N_EPOCH_INTERP`, `NGMIX_MIN_DIST_INTERP` and `NGMIX_MIN_DIST_NOISEFILL` describe defects anywhere in kept stamps, including off-footprint defects close to the centre.
  They support stricter downstream cuts.
  An absent distance is `-1`.

A zero-weight pixel with flag `0` contributes no instrument bit.
Integral floating-point flags are accepted; non-finite, fractional, negative or int64-overflow values on the footprint raise rather than silently corrupting the OR.
FITS shear extensions and the final catalogue retain the diagnostic columns and their integer/sentinel semantics.

## Reproduction

Run validation only inside a compute allocation and the campaign SIF, with `PYTHONPATH=src:.`.
`d0_prevalence.py` invokes the real `prepare_postage_stamps` path on retained survey stores, without running galaxy fits or writing to inputs.
It reads SQLite values through immutable read-only connections and reports object/epoch edge-case frequencies, rejection counts, and exact agreement with the survey's cached raw-overlap decisions.
It also reports differences from the survey's clipped-footprint fallback policy.

```bash
python scripts/validation/masked_pixels/d0_prevalence.py \
  --inputs /scratch/cdaley/veto-survey/prevalence_checks.json \
  --survey /scratch/cdaley/veto-survey/prevalence.npz \
  --out /scratch/cdaley/veto-d0/188.307.json --nproc 4
```

`tile_defect_flags.py` measures kept-epoch diagnostics and can overlay them on an exposure-flag map.
Its `--plot-only` mode regenerates the plot from saved arrays.
The independent-pixel Poisson reference uses `1 - exp(-0.0028 * A * N)` with each object's footprint area and kept-epoch count.
Spatially correlated defects, non-random detection and conditioning on epoch selection violate that reference's independence assumptions; it is not a calibration.

## Regression coverage

The unit tests cover own-footprint overlap, adjacent non-overlap, another object's footprint, raw defect sources, border touch and unavailable footprints.
The module tests cover both fills, flipped/unflipped CCDs, cut ordering and counters, off-tile defects, neighbour masks, retained stamp-wide diagnostics, footprint-record OR/sentinel semantics, and catalogue round trips.
The workflow test resolves detection for data and image simulations under each blend handling and checks that `SEG_VIGNET` reaches `Tile_cat` and the D0 predicate.

`mutate_d0.py` restores a fixed 7/13-px rule in a private source copy and requires assertion failures from the footprint tests.
It never edits the checkout.

```bash
python scripts/validation/masked_pixels/mutate_d0.py --out /scratch/cdaley/veto-d0/mutant.log
```
