# Tile-object defect flags

We propagate exposure instrument bits onto each tile object's coadd segmentation footprint, using only epochs that pass the masked-fraction cut and central veto.
The implementation does not change fits or impose a catalogue cut.
This change stacks on #925 (`feat/uberseg-seg-vignet`), which puts the int32 `SEG_VIGNET` column in the committed tile catalogue under every blend handling and makes `Tile_cat` read it.
Stacking only on #922 would leave the committed workflow without the footprint input and produce `-1` for fitted objects.

- `NGMIX_IMAFLAGS_ISO` is the bitwise OR of exposure flag pixels on `tile_seg == obj_id`, OR'ed across kept epochs, with every instrument bit preserved.
- `NGMIX_N_EPOCH_INTERP_ISO` counts kept epochs containing at least one interpolated defect pixel on that footprint.
- Both are `-1` when segmentation is missing, the footprint is empty, or it touches a stamp edge.
  We cannot establish completeness in those cases.
  Never-fit final-catalogue rows retain the existing count convention of `0`.
- A zero-weight/off-tile pixel with flag `0` contributes no instrument bit, but contributes to the epoch count if interpolated.
  Clean quarter-turn copies made by defect weighting do not count as defects.

The coadd segmentation is MegaCam-flipped to each exposure's stamp orientation, without reprojection, matching the existing segmentation overlay.
The columns supplement the stamp-wide diagnostics rather than replacing the tile SExtractor's `IMAFLAGS_ISO`.

## Fixture validation

The synthetic path runs `prepare_postage_stamps` with multiple epochs and both flipped and unflipped CCDs.
It checks inside/outside-footprint defects, clean objects, unflagged zero-weight defects, rejected epochs, missing/empty/clipped footprints, safe float-to-integer conversion, integer FITS round trips in every shear extension, and final-catalogue values and never-fit sentinels.

| Mutation | Observed failure | Result |
| --- | --- | --- |
| Reduce/count over the entire stamp | Outside-footprint case becomes `1028 / 2`, expected `0 / 0` | 2 assertion failures |
| Overwrite the OR on each epoch | Final clean epoch erases `1028`, leaving `0` | 2 assertion failures |
| Accumulate before either epoch cut | Rejected epochs contribute bit `64` and interpolation counts | 8 assertion failures |
| Drop the make_cat write mapping | Final-catalogue value stays `0`, expected `1028` | 1 assertion failure |
| Remove the workflow's SEG_VIGNET export | Committed detection disables the footprint | 6 assertion failures |
| Disable SEG_VIGNET in the committed ini | Committed detection disables the footprint | 6 assertion failures |
| Drop Tile_cat's SEG_VIGNET read | Ordinary clean objects receive `-1` instead of `0` | 6 assertion failures |

All seven mutations were restored before running the suite.
`test_committed_detection_supplies_defect_flag_footprint` resolves the actual workflow and ini for data and image simulations under default, noisefill and uberseg handling, then checks the `SEG_VIGNET` load and a clean object's zero-valued diagnostics.
The runtime SIF lacks pytest; validation uses `shapepipe-develop-240b37e4.sif` on Slurm allocations under `def-mjhudson`.
The rebased branch's ngmix, make_cat, unit and PSF property suite passes 830 tests with one skipped; the workflow/DAG/params-pin suite passes 31 tests with one skipped.
The seven mutations were rerun on the rebased branch on allocation `24920131`, with source edits and pytest on the same compute node and bytecode caches cleared between mutations.
The params-pin test runs inside that image with the existing `/project/6001537/cdaley/snakemake-env/lib/python3.12/site-packages` prepended and pytest plugin autoload disabled, because the image lacks Snakemake.
Regeneration leaves the pin unchanged.
`completeness.py` counts files, not columns; this change adds no product file.

## Flag dtype

`split_exp.py` casts flag images to `int16`, and `get_stamps` preserves the input dtype.
On real exposure `2254826p.flag.fits.fz`, CCD HDU 1 is `uint8`, with values `0, 1, 3, 8, 10, 11`; applying the pipeline's split cast and stamp extraction produces `int16` stamps.
This checks the extraction path, **not a retained production SQLite flag vignet**, which was not available.
The fixture also exercises integral `float32` flag stamps.
The new reduction validates finite, nonnegative, integral, int64-representable values before casting to `int64`, so floats cannot silently truncate or overflow.

## Real-tile validation remains pending

**No real-tile fractions or spatial plot have been measured for this change.**
The retained final catalogues and `nflagged` maps are available, but a complete tile with its image, flag, weight, PSF and segmentation stores has not been located.
For example, the `smk-ab925-us` tile `181.308` manifest points to `/local/scratch/sp-181.308`, and the shared tile directory contains `cleaned.json` and only the exposure-finder output.
`fir-b0001.yaml` likewise points to node-local inputs/stores and requests `clean_tiles: true`.
The `ngmix-fixtest` inputs for `196.305` are broken symlinks to removed `p3-first10` stores.
A bounded search of the shared products, scratch run directories and test runs did not locate a complete alternative; this is not proof that no retained tile exists elsewhere.

A private regeneration config is saved at `/scratch/cdaley/defectflag-check/run/config.yaml`, with separate work, products, index and tile-store directories under that prefix and both cleaners disabled.
It selects tile `181.308`, imports its retained PSF tar directories through read-only symlinks, and reads the shared raw mirrors.
The dry run targets only `manifests/tile_vignets.json`; its command and log are `run/dryrun.sh` and `run/dryrun.log` beneath the same prefix.
The plan contains **35 jobs: 9 exposure retrievals, 9 exposure splits, 9 exposure PSF chains, and 8 tile-level jobs**.
All nine source exposure output directories are absent, consistent with the recorded cleanup; the PSF tarballs do not include the split images, headers or SExtractor background maps needed for the vignets.
The run stopped at the dry run, as required when any exposure-level job is planned.
No split, detection, PSF fitting, vignette generation or ngmix job was executed.

`tile_defect_flags.py` opens stores read-only and calls the real `prepare_postage_stamps` path without fitting.
Its module docstring gives the full command.
It reads `SEG_VIGNET` through the rebased branch's `Tile_cat`, without a separate segmentation path or a script-only fallback.
It writes per-object arrays, counts, flag-bit breakdowns, explicit denominators, the edge-touching and empty-footprint counts, the unknown-footprint fraction among objects with kept epochs, and two Poisson reference estimates under the requested scratch output prefix.
With `--nflagged` and `--plot`, it overlays flagged objects on the exposure defect map.
A synthetic one-object end-to-end smoke test produced bit `4` and one interpolated epoch; those are fixture results, not survey measurements.

For the real check, compare the measured instrument-bit fraction with the mean of `1 - exp(-0.0028 * A * N)` over known footprints, using each object's area and kept-epoch count.
That expression ignores spatial correlations and the central veto: every kept interpolated defect is at least 7 px from the stamp centre, so compact footprints can be protected completely.
The script also reports a reference using only footprint area outside that radius; even this is not a prediction calibrated for the epoch-selection correlations or invalid-RMS/noise-filled defects.
Explain the discrepancy rather than treating agreement with an independent-pixel approximation as a requirement.

Back-filling would require retained or regenerated stores, but would not require shear refits.
No back-fill has been attempted.
