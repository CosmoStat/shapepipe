# The Snakemake workflow

Production ShapePipe runs — a tile list, not a tile — are orchestrated by the
Snakemake workflow in `workflow/`. It is the current and only orchestration for
real data: the CANFAR submission layers that preceded it were retired at
`2ef07e45`, while the bit-coded bash job scripts they drove
(`run_job_sp_canfar_v2.0.bash`, `job_sp_canfar_v2.0.bash`) survive only as the
entry point of sp_validation's image-simulation chain.

**Module code is untouched by it.** Every rule ultimately calls
`shapepipe_run -c <config>` on the config chains described in
[Configuration](configuration.md); the workflow supplies the scheduling, the
per-unit isolation, the bookkeeping and the reclamation around them. Reading
[Basic execution](basic_execution.md) first is still the right order.

This page is the orientation. The deep reference — design rationale, every
rule, the profile, the two-roots argument — is
[`workflow/README.md`](https://github.com/CosmoStat/shapepipe/blob/develop/workflow/README.md)
in the repository, and the living design record is
[issue #848](https://github.com/CosmoStat/shapepipe/issues/848).

## Quick start

```bash
# One-time, on a SHARED filesystem: the SLURM executor re-invokes this python
# inside every job, so a node-local path (/tmp) will not do.
uv venv /project/<alloc>/<user>/snakemake-env --python 3.12
source /project/<alloc>/<user>/snakemake-env/bin/activate
uv pip install 'snakemake>=9,<10' 'snakemake-executor-plugin-slurm>=2.7,<3'

# Edit workflow/config.yaml: tile_list, run_dir, products_dir, container,
# star_cats.

workflow/bin/sp run       # bring products on disk up to date with the tile list
workflow/bin/sp report    # the success/failure tables, any time, mid-run is fine
workflow/bin/sp container status   # which image the jobs will run inside
workflow/bin/sp cancel <name>      # scancel this workflow's jobs
```

`workflow/bin/sp` is the entry point for everything. It loads the apptainer
module and the venv, snapshots the code it is about to launch, and sets the
state directory, the SLURM profile and `SP_PHASE`. **Bare `snakemake` outside
`sp` is unsupported** — it skips all of that, including the phase variable the
Snakefile needs to build its index at parse time. Any verb `sp` does not
recognise passes straight through to snakemake with the right profile and state
directory, which is the escape hatch for `sp --unlock`, `sp --dag` or
`sp exp_psf ...`.

## The two static invocations

The exposure job set is *derived from data the run itself produces* — a tile's
`find_exposures` output — so it cannot live in the same static DAG that produces
it. `sp run` is therefore two snakemake invocations over one Snakefile:

1. **PREPARE** (`snakemake prepare_all_tiles`) — the per-tile static chain
   `tile_get_images → tile_uncompress → tile_find_exposures`, with `keep-going`
   so tile failures stay independent. A nonzero exit is a warning, not a fatal
   error: tiles that lost their exposure list are dropped at the compute parse,
   and `SP_MISSING_THRESHOLD` is the real gate on how many may be lost.
2. **COMPUTE** (`snakemake all`) — this invocation's *parse* builds the
   tile↔exposure index (`run_index.sqlite`) and the DAG then runs the full
   exposure and tile chains:

   | level | chain |
   | --- | --- |
   | exposure | `exp_get_images → exp_star_cat → exp_split → exp_mask → exp_psf → exp_persist → exp_footprint → clean_exposure` |
   | tile | `tile_exp_forest → tile_merge_headers → tile_detect → tile_vignets → tile_ngmix → tile_merge_cats → tile_make_cat → clean_tile` |
   | campaign | `star_catalogue`, `coverage_map` |

`sp run` chains both and fails if either phase did. The index **accumulates**
across invocations, so appending tiles to `tile_list` later changes which jobs
exist without invalidating completed work.

Every rule instance is its own SLURM job (`executor: slurm` in
`profiles/nibi/config.yaml`) carrying that rule's own attempt-scaled resources,
and runs inside the container via snakemake's software-deployment method — the
workflow never calls apptainer itself.

## Manifests are the currency

The atom is **one rule == one `shapepipe_run` on one unit** (one tile, or one
exposure), in its own sharded store: `<run_dir>/tiles/<2-char prefix>/<ID>/` or
`<run_dir>/exp/<prefix>/<base>/`.

A rule's single declared output is not its product files but that unit's
**manifest**, `<unit dir>/manifests/<stage>.json`. Per-CCD declaration would
mean millions of paths at DR6 scale, and a missing CCD is frequently
legitimate.

Manifests are *success-only*. After each run, `completeness.py check` writes its
full verdict — per-runner counts against floors, scraped failure reasons, the
`shapepipe_run` exit status — to the rule's `log:`
(`<unit dir>/logs/<stage>.json`) **every** time, and additionally to the
declared manifest **only** when that verdict is a success. So a manifest on disk
means "this stage succeeded", and a resume after an unclean death can never
schedule downstream work on top of a failure.

Completeness is a **count floor, not a taxonomy**: a stage is a real failure iff
a mandatory runner produced fewer products than its floor. Per-CCD attrition
between the floor and the expected count is tolerated and recorded. There is no
error-signature whitelist.

## `sp report`

`sp report` reads the index and those two directories — the manifest for
success, the log for failure, a unit with neither ran nothing — and prints the
per-stage tables plus a `run_report.json` beside the index. It is deliberately
**not** a DAG node: a report rule declaring every tile's output would be a
descendant of every job, so one hard failure under `--keep-going` would poison
its cone and the report would never run in exactly the situation it exists for.
It is a plain script, runnable at any time including mid-run, and the COMPUTE
invocation's `onsuccess`/`onerror` hooks fire it automatically so every run ends
with one.

## `sp container`

`sp container` owns which image the jobs run inside, in one resolution order
shared by the CLI and the workflow: **sandbox → your cached SIF → the
`container:` path in `config.yaml`**.

```bash
sp container status                  # layers present, active one, revision vs HEAD
sp container pull                    # ghcr.io/cosmostat/shapepipe:develop-runtime
sp container sandbox                 # unpack the SIF writable (opt-in)
sp container exec --writable pip install <pkg>
sp container resolve                 # just the path the workflow will run
```

The cached SIF is a private pull, so nobody else's refresh moves the ground
under your running jobs; the sandbox is the escape hatch for work needing a
package the image does not carry yet. `status` places the image's
`org.opencontainers.image.revision` label against this checkout's HEAD.
`pull` needs the network — run it on a login node, never from a batch job.

See [Container workflow](container.md) for how the images are built.

## Durable products

Two roots. `run_dir` is scratch: the bulk per-unit stores, sized to finish
inside the purge window. `products_dir` is persistent and holds the low-volume
things worth keeping:

- **Final catalogues** — `<products_dir>/tiles/<prefix>/<ID>/final_cat-<ID>.fits`,
  mirroring the scratch shard structure. This file is also the *tile-finished
  marker* that lets a completed tile stop pinning its exposures, which is the
  second reason it cannot live on scratch.
- **PSF products** — `exp_persist` packs the files named by `persist_exp:` in
  `config.yaml` (by default the psfex_interp `validation_psf-*.fits`, the rho/tau
  statistics input) into one uncompressed tar per exposure, plus a manifest
  recording the members. It answers the scratch purge, not the workflow's own
  reclamation, so it runs whether or not `clean:` is on.
- **The index and the report** — `run_index.sqlite`, `missing.json`,
  `run_report.json`.
- **Coverage** — see below.

Exposure stores are reclaimed by the in-DAG `clean_exposure` rule rather than by
`temp()`: exposures overlap tiles, so `temp()` would cascade destructive reruns
whenever a tile is appended. Reclamation leaves a `cleaned.json` tombstone
holding the absorbed manifests, so a reclaimed exposure still reports its
per-runner counts and blocks no tile.

## Coverage

Sky coverage is a workflow product, built from records the DAG already writes
rather than by scraping a finished campaign.

`exp_footprint` writes one JSON per exposure to
`<products_dir>/exp/<prefix>/<base>/manifests/exp_footprint.json` giving the
four sky corners of every CCD that got a PSF model. The valid-PSF CCD set comes
off `exp_persist.json`'s tar members — exact, because `psfex_interp` returns
*without* writing `validation_psf-*.fits` on NOT_ENOUGH_STARS, BAD_CHI2 or
FILE_NOT_FOUND — and the WCS off the `headers-<exp>.npy` written by `exp_split`.
That makes `validation_psf-*` in `persist_exp:` a precondition of the chain, not
a preference.

Set `coverage: {enabled: true}` and one further job, `coverage_map`, stamps every
footprint into `<products_dir>/coverage/coverage.hsp` — a HealSparse map counting,
per sky pixel, the exposures with a valid PSF there. It is **campaign-cumulative**:
its declared inputs are the in-scope footprints, but the script reads every record
on the products root, reclaimed exposures included, so appending tiles grows the map
instead of replacing it. `nside` is set in `config.yaml` to the production
128/131072 pair, chosen to align pixel-wise with the UNIONS bit masks — nothing
defaults to it, and a coarser map would look plausible and silently fail to align.

Plotting stays out of the DAG, as a human act on a durable product:

```bash
plot_coverage_map -i <products_dir>/coverage/coverage.hsp ...
```

with the sky windows under `coverage.plot` in `config.yaml`.
