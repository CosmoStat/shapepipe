---
name: verify
description: >
  Run the ShapePipe stage a change touches through workflow/bin/sp on a small
  real-data fixture (1 tile, 3 exposures) and compare the output (columns,
  counts, star/galaxy metrics, and a figure) against a known-good baseline.
  Use before committing or opening a PR on any science or workflow change.
---

# Verify a ShapePipe change

Before calling a change done, run the stage it touches through the real launcher
(`workflow/bin/sp`) on the small fixture and inspect what came out beside the
known-good baseline. Passing unit tests confirm contracts still hold; running
the fixture confirms the pipeline still produces what you meant.

## One command

Run inside a SLURM allocation (never on a login node):

```bash
uv run workflow/scripts/verify.py --stage <stage> --baseline <baseline-dir> --out <out-dir>
```

The script has no harness dependence — a human or any agent runs the same
command. It writes `verify_summary.json` and a comparison figure
(`verify_<stage>.png`) to `--out`, and prints the counts, column diffs, and
star/galaxy metrics to stdout.

## Which stage to run

Pick the earliest stage that exercises your change:

| Stage | Rules run through `workflow/bin/sp` | Covers | Measured wall time (16 cores) |
|---|---|---|---|
| `prepare` | `tile_get_images` → `tile_uncompress` → `tile_find_exposures` | Tile ingestion, exposure discovery | ~20 s |
| `detect` | `prepare` + `tile_get_catalogue` → `tile_detect` | SExtractor detection, external catalogue join, `SEG_VIGNET` | ~1–5 min |
| `psf` | `prepare` + `exp_get_images` → `exp_split` → `exp_psf` → `exp_persist` → `star_cat_merge` | Star selection, masks, PSFEx fit, PSF validation, star catalogue | ~2.7 min (0.7 core-h) |
| `vignets` | `psf` + `detect` + `tile_exp_forest` → `tile_merge_headers` → `tile_vignets` | PSF interpolation at galaxy positions, multi-epoch postage stamps | ~4–5 min |
| `all` | Full `sp run` through `tile_ngmix` → `tile_merge_cats` → `tile_make_cat` → `final_cat_merge` + `exp_maps` | Everything, including ngmix metacalibration and final catalogue | ~14–20 min (3.7–5.3 core-h) |

Pass `--strict` on a pure refactor to exit non-zero on any numeric or schema
difference from the baseline.

## The fixture

Default tile: **`592.232`** (from the DR6 patch testbed), which touches **3
exposures** (`2514827` at 0.86″ seeing, `2780860` at 0.49″, `2780862` at 0.49″)
— ~1.6 GB of raw inputs, ~33,800 detected objects (~2,080 v2.0.0-selected
galaxies), and ~5,080 good reserve stars. One tile exercises both sharp
sub-0.55″ seeing and typical seeing, multi-epoch overlap (`N_EPOCH >= 2`), and
all exposure and tile rules.

Pixels are never committed to the repository.

- **On Fir:** inputs are already staged at
  `/project/rrg-kyi/cdaley/unions/analyses/shapepipe-patch-testbed/inputs`,
  sky masks at `/project/def-waerbeke/cdaley/fir-assets/masks/dr6-2026-09-29`,
  the container at
  `/project/def-waerbeke/cdaley/fir-assets/shapepipe-develop-240b37e4.sif`,
  and the `develop` baseline at
  `/project/rrg-kyi/cdaley/unions/analyses/shapepipe-verify-skill/baseline_develop`.
  `verify.py` picks these up automatically when present; override with
  `--inputs`, `--masks`, `--container`, `--baseline` (or `SP_VERIFY_INPUTS`,
  `SP_VERIFY_MASKS`, `SP_CONTAINER`, `SP_VERIFY_BASELINE`; pass `--no-baseline`
  when creating a fresh baseline).
- **Elsewhere:** stage tile `592.232` and its three exposures (`2514827p`,
  `2780860p`, `2780862p`) from CADC/VOSpace into `<root>/{tiles,exposures,catalogues}`
  (tile image + weight, external `.cat` + `.seg.fits.fz`, and each exposure's
  `.fits.fz`, `.flag.fits.fz`, `.weight.fits.fz`), then point `--inputs <root>`
  at it and run once on an unmodified `develop` checkout to create your
  `--baseline` directory.

## What to look at before committing

1. **Read the printed summary and open `verify_<stage>.png`.**
   - For a **refactor**: verdict should be `NO_DIFFERENCE` (note that `tile_ngmix`
     has ≤1 object in ~34k with convergence jitter on rare runs; verify any
     non-zero diff is single-object ngmix jitter and not a systematic shift).
   - For a **science change**: verify the output moves *where intended* (e.g.
     star validation residuals move while galaxy shapes stay bit-identical for a
     validation-only change) and *nowhere else*.
2. **State how the change was verified** in the commit message or PR description
   (stage run, baseline compared against, key numbers, and attach
   `verify_<stage>.png` to PRs that move a scientific quantity).
