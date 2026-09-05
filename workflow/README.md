# ShapePipe Snakemake orchestration

Snakemake workflow that orchestrates real-data ShapePipe runs. It replaces the
`curl_canfar_local.sh → run_job_sp_canfar_v2.0.bash → job_sp_canfar_v2.0.bash`
bash layers and the per-site sbatch reimplementations. **Module code is
untouched**: rules call `shapepipe_run -c <config>` on the existing config
chains. Design and rationale:
[CosmoStat/shapepipe#848](https://github.com/CosmoStat/shapepipe/issues/848)
(the design document).

Use `workflow/bin/sp` for everything. Bare `snakemake all` outside `sp` is
unsupported: `sp` sets the state directory, the SLURM profile, and
`SP_PHASE`, which the Snakefile needs to build the tile/exposure index at
parse time. Running snakemake directly skips all of that.

## Quick start (nibi)

```bash
# One-time: a snakemake env on a SHARED filesystem (/project — NOT /tmp, which
# is node-local; the SLURM executor re-invokes this python inside every job).
uv venv /project/def-mjhudson/cdaley/snakemake-env --python 3.12
source /project/def-mjhudson/cdaley/snakemake-env/bin/activate
uv pip install 'snakemake>=9,<10' 'snakemake-executor-plugin-slurm>=2.7,<3'

# Edit workflow/config.yaml: tile_list, inputs.tiles/exposures, outputs.run_dir,
# outputs.products_dir/index_db, and container.

# `psf_model` is `psfex` or `mccd`; mccd is wired but unvalidated here, while psfex is exercised by smk-g4 through smk-g6.

# The committed launcher loads apptainer/1.4.5 + the /project venv, so a
# fresh shell always has the right state.
workflow/bin/sp run       # bring products on disk up to date with the tile list
workflow/bin/sp report    # emit run_report.json now (mid-run is fine)
workflow/bin/sp cancel <run-name-substring>   # scancel this workflow's jobs
workflow/bin/sp container status              # which image the jobs will run
```

Installed and pinned versions on nibi (`/project/def-mjhudson/cdaley/snakemake-env`,
queried 2026-07-30): `snakemake==9.23.1`, `snakemake-executor-plugin-slurm==2.7.1`.
Pin range: `snakemake>=9,<10`, `snakemake-executor-plugin-slurm>=2.7,<3`. The
v8→v9 breaks matter here: `--use-singularity` became `--sdm`, executors became
plugins, and full `rerun-triggers` became the default.

Anything other than `run`, `report`, `container`, `cancel` passes straight through to
snakemake with the workflow's profile and state dir — the direct command path for
`sp --unlock`, `sp --dag`, `sp exp_psf ...`.

## The container image

`sp container` owns which image the jobs run inside. Two layers, and the second
only exists if you ask for one:

* your **cached SIF** (`~/.cache/shapepipe/shapepipe.sif`, `SP_CACHE_DIR` or
  `SP_CONTAINER` to move it) — a read-only pull of the published image, private
  to you, so nobody else's refresh moves the ground under your running jobs;
* an optional **sandbox** (`~/.cache/shapepipe/sandbox/`, `SP_SANDBOX`) — the
  same image unpacked writable, so a `pip install` into it sticks. The escape
  hatch for work needing a package the image does not carry yet.

The Snakefile's `container:` is that resolution, in one order shared by the CLI
and the workflow: **sandbox → cached SIF → the `container:` path in
`config.yaml`**. With an empty cache — the normal case — that lands on the
shared `/project` `.sif` the workflow has always used, so this changes nothing
until you opt in.

```bash
sp container status                      # layers present, active one, revision vs HEAD
sp container pull                        # ghcr.io/cosmostat/shapepipe:develop-runtime
sp container pull --tag docker://...     # some other image
sp container sandbox                     # unpack the SIF writable (opt-in)
sp container exec --writable pip install <pkg>
sp container exec python -c 'import shapepipe'
sp container resolve                     # just the path the workflow will run
```

`status` reads the image's OCI labels and places its
`org.opencontainers.image.revision` against this checkout's HEAD:
in-sync / behind / ahead / diverged, or unknown when the image carries no label
or the commit was never fetched here.

**`pull` needs the network.** Compute nodes on Alliance clusters generally have
none, so run it on a login node or inside an `salloc` allocation — never from a
batch job. `pull` and `sandbox` both stage to a sibling path and swap it in, so
an in-flight job never sees a half-written image and a failed rebuild leaves the
one you had intact.

## Execution: two static invocations

The exposure job set is data-derived from the tiles' `find_exposures` output,
so it cannot live in the same static DAG that produces it. `sp run` is
therefore two snakemake invocations over one Snakefile:

1. **PREPARE** — `snakemake prepare_all_tiles`: per-tile static DAG
   (`Git_vos → Uz → Fe`), `keep-going` so tile failures are independent. A
   nonzero exit here is not fatal to the run — tiles that lost their
   exposure list are dropped at the compute parse — but it is a warning:
   `SP_MISSING_THRESHOLD` (default 0.0) is the real gate.
2. **COMPUTE** — `snakemake all`: this invocation's *parse* builds the
   tile↔exposure index (`build_index.py`, imported at parse time, not a DAG
   node) and runs the full tile/exposure compute chain. The index
   accumulates across invocations, so appending tiles later changes which
   jobs exist without invalidating completed work.

`sp run` chains both so the UX is one command; both exit codes are checked
and the run fails if either phase failed.

## The launch code snapshot

`sp run` copies the code it is about to launch — `workflow/` (config symlinks
dereferenced), `src/`, the repo's `scripts/` (`final_cat_merge` loads
`scripts/python/create_final_cat.py` by path) and the profile — into
`<state dir>/code`, records HEAD
plus a dirty flag in `<state dir>/code/snapshot.json`, and runs the campaign
entirely out of that copy. It matters because a campaign is not one process: the
SLURM executor re-invokes snakemake on every job's node, so jobs re-parse the
Snakefile and read `workflow/scripts/*` and the ini chain hours after launch.
**Editing the checkout while a campaign runs is therefore harmless; a change
takes effect on the next `sp run`.**

Everything workflow-internal hangs off `workflow.basedir`, which *is* the
snapshot, including the profile's `PYTHONPATH`
pin, which YAML cannot interpolate: `sp run` rewrites that single path in the
snapshot's copy of the profile and launches `--profile` at the copy. The snapshot
is refreshed wholesale on every `sp run` — a new `sp run` *is* the relaunch — and
the other verbs run out of the existing snapshot, `sp container` excepted (it is
about the image you are working with now, not about a campaign). The mechanism
and its rationale live in one place: `bin/sp`.

## Execution mode: one SLURM job per rule

The profile (`profiles/nibi/config.yaml`) sets `executor: slurm`. Every rule
instance becomes its own SLURM job carrying that rule's own attempt-scaled
resources (`cpus_per_task = threads`, `mem_mb`, `runtime`) and runs inside
`apptainer exec` via the profile's software-deployment method — the workflow
never calls apptainer directly. Snakemake feeds the queue as jobs finish, so
the full campaign's job count never needs to queue at once, and multi-node
scaling is inherent. `jobs:` in the profile caps concurrent submissions at a
fraction of the cluster's per-user submit limit (queried, not invented — see
the comment in the profile file for the query and date).

The profile intentionally sets no `set-resources` / `set-threads` overrides:
those replace a rule's own values wholesale, which would kill the
attempt-scaled `mem_mb = lambda wc, attempt: ...` OOM retries and the tuned
ngmix thread count. Rules own their own resources; the profile only supplies
defaults for rules that state nothing.

`group:` labels that fuse short rules (uncompress, merges) into their chunky
neighbours (queue-latency amortization, per the PRD) are **not yet wired**:
they require labels in `workflow/rules/*.smk`, out of scope for this
profile-only pass.

## Layout

```
workflow/
  Snakefile              parse-time index load; global container:; onsuccess/onerror report hooks
  config.yaml            the run: tile list, input/output paths, container, chunk count
  bin/sp                 committed launcher (module load + /project venv + launch code snapshot + run/report/container/cancel)
  rules/
    prepare.smk          tile get_images/uncompress/find_exposures
    exposure.smk         per-exposure: get_images, split, psf, persist, footprint (no temp()); campaign star_cat_merge
    tile.smk             per-tile: exp forest, merge_headers, detect, vignets, ngmix, merge, make_cat; campaign final_cat_merge
    coverage.smk         campaign-level: the HealSparse nexp mask from the exposure footprints
  scripts/
    build_index.py       prepare-phase run_index.sqlite builder (plain script)
    build_forest.py      per-tile exposure symlink forest (group-compatible shell)
    completeness.py      the ported count table (shared by sp_rule + run_report)
    run_report.py        standalone report (NOT a DAG node; run_report hooks call it)
    container.py         image layers + the resolution order behind `sp container` (stdlib-only)
    persist_exp.py       ONE exposure's keepable PSF products -> one tar on products_dir (the exp_persist rule)
    exp_footprint.py     ONE exposure's per-CCD sky corners, for the CCDs with a PSF (the exp_footprint rule)
    coverage_map.py      every exposure footprint on products_dir -> coverage.hsp (the coverage_map rule)
    hdf5_reconcile.py    bring an hdf5 catalogue into agreement with a campaign (shared by both merges)
    merge_star_cat.py    ALL exposures' validation_psf, out of the tars -> full_starcat_<campaign>.hdf5
    merge_final_cat.py   ALL tiles' final_cat -> final_cat_<campaign>.hdf5 (the final_cat_merge rule)
    clean_exposure.py    ONE exposure's store + manifests + logs -> tombstone (the clean_exposure rule)
profiles/nibi/config.yaml  SLURM executor; apptainer SDM; per-user jobs cap; keep-going
```

## How it works

- **The atom is one rule == one `shapepipe_run` on one unit.** Its single
  declared output is that unit's manifest
  (`<unit dir>/manifests/<stage>.json`), not its product files — a missing
  CCD is often legitimate, and at DR6 scale per-CCD declaration means
  millions of paths.
- **Manifests are the DAG's currency, and they are success-only.**
  `completeness.py check` writes its full verdict — per-runner counts against
  expected counts, scraped failure reasons, the `shapepipe_run` exit status when nonzero
  — to the rule's `log:` (`<unit dir>/logs/<stage>.json`) on *every* run, and
  additionally to the declared `<stage>.json` manifest only when that verdict is
  a success. So `<stage>.json` on disk means "this stage succeeded", and a
  resume after an unclean death cannot schedule downstream work on top of a
  failure. Nothing unlinks anything: Snakemake deletes a failed job's declared
  output natively and never touches its log, which is why the profile runs
  *without* `keep-incomplete`. `sp report` reads both dirs — the manifest for
  success, the log for failure — and a unit with neither ran nothing.
- **Completeness is an exact-count check, not a taxonomy.** After a run,
  `sp_rule.py` counts products per mandatory runner against
  `completeness.py`'s expected count and exits nonzero below it. A shortfall
  below `expect` fails unless `warn` is set. No 3-class taxonomy, no
  error-signature whitelist. `--keep-going` isolates a failure to its own
  DAG cone.
- **Stores are sharded.** Every tile/exposure runs its own `shapepipe_run`
  in `tiles/<2-char prefix>/<ID>/` or `exp/<prefix>/<base>/`. Configs are
  committed under `workflow/config/cfis/` and version with the rules that set
  the env vars they interpolate — there is no `config_src` knob, and no per-unit
  config symlink; `$SP_CONFIG` points straight at the committed directory.
- **There is no masking stage, on either side.** ShapePipe generates no masks
  (PR #847). The one mask that reaches pixels is the instrument flag image
  delivered with each exposure, which `exp_split` splits per CCD beside image
  and weight and SExtractor reads as `IMAFLAGS_ISO`. Everything else — star
  halos, manual masks, per-band coverage, MaxiMask — is supplied as sky-fixed
  healsparse maps and QUERIED once per object: the `mask_query` module writes a
  `MASK_EXT` column onto each CCD's detection catalogue (inside `exp_psf`), and
  `make_cat` writes one `MASK_<band>` column per band onto the final catalogue
  (inside `tile_make_cat`). Neither is cut on in the pipeline: instrument flags
  mark corrupted measurements and are the only masks that reject anything here,
  while the healsparse masks are location flags and every decision about them is
  downstream. On exposures the query ships off — `MASK_PATHS` is commented out,
  making `mask_query` a no-op pass-through — so turning it on is a config edit,
  not a chain edit. Map paths are config, not
  code, so regenerated products cost a config edit. Nothing is fetched from a
  catalogue server, staged, or rasterized, which is why the old
  `star_catalogue` / `exp_star_cat` / `exp_mask` rules and their cache root are
  gone.
- **The index is parse-time data, never a rule input.** Appending tiles
  changes which jobs exist without invalidating completed work.
- **Exposure products are not `temp()`.** Exposures overlap tiles, so
  `temp()` would cascade destructive reruns when a tile is appended later.
  Reclamation is the in-DAG `clean_exposure` rule instead: one job per
  exposure, taking every consuming tile's `tile_vignets` manifest as input
  (the campaign-wide consumer set comes from the accumulating index), which
  deletes the store *and* the exposure's `manifests/` and `logs/`, and leaves a
  `cleaned.json` tombstone. Deleting the manifests is what makes a late append
  correct: the
  appended tile finds an unbuilt chain and regenerates it. The `clean:` flag
  in `config.yaml` gates it; flipping it on later reclaims retroactively,
  since the missing tombstones schedule exactly the outstanding clean jobs.
  The tombstone is written *before* anything is deleted, so a crash can cost
  disk but never the record.
- **A finished tile declares no reclaimed exposures.** Deleting an exposure's
  manifests would otherwise rerun every other tile that reads it, and those
  reruns spread across the exposure-overlap component. So a tile whose
  `final_cat` exists drops the exposure manifests that are gone from its input
  list, and holds the rest through `ancient()`. This is why the profile runs
  with `rerun-triggers: [mtime, params, code, software-env]`: the `input`
  trigger reads that cut as a reason to rerun the very tiles it protects.
  Know the consequence — `--forcerun` on a tile whose `final_cat` exists will
  not rebuild its reclaimed exposures. Delete the `final_cat` first.
- **PSF products leave scratch before the purge does.** `exp_persist` packs
  the products named by `persist_exp:` in `config.yaml` from the exposure's
  scratch store into ONE uncompressed tar,
  `<products_dir>/exp/<prefix>/<base>/psf/<base>.tar` (inodes, not bytes, bind
  on /project), and writes ONE manifest beside it recording the patterns, the
  members and their sizes. The
  threat it answers is the /scratch purge, not `clean_exposure` — the store goes
  in 60 days whether or not the workflow reclaimed it — so it runs even with
  `clean: false`, requested directly by `rule all`. `clean_exposure` takes its
  manifest as an input, so reclamation can never overtake the copy. It is a
  rule of its own rather than a `cp` on the end of `exp_psf` because the keep
  list rides on `params`: adding a pattern reruns seconds of packing, not four
  hours of PSF fitting per exposure. A pattern that matches nothing is a
  recorded warning (setools rejects sparse CCDs); matching nothing at all is a
  failure. A `localrule`, by the same arithmetic as `clean_exposure`.
- **The star catalogue's inputs are always kept; `persist_exp:` is what you
  keep on top.** `exp_persist` packs `psf_validation` — the psfex_interp
  validation catalogue, one per CCD — for every exposure whatever the config
  says, because `star_cat_merge` stacks exactly those into the campaign's
  `full_starcat`. They are that catalogue's provenance, and they are what keeps
  appending a tile next month cheap rather than a rebuild from VOS. About 2 MB
  per exposure: ~40 GB and ~40k inodes at DR6 scale, against a ~1 M-inode group
  quota. `persist_exp:` is purely additive, and an empty list is legal — the tar
  then holds the merge's inputs and nothing else.
- **The keep list names products, not globs.** Entries are names from a
  catalogue in `workflow/scripts/persist_exp.py`, which is the single source of
  truth for what each one means and what keeping it buys
  ([#844](https://github.com/CosmoStat/shapepipe/issues/844)); `config.yaml`'s
  block is that catalogue rendered, and `persist_exp.py --list-products` prints
  it. Sizes are per exposure, 40 CCDs, measured on smk-m2.

  | product | glob | per exposure | what it buys |
  |---|---|---|---|
  | `psf_model` | `*.psf` | 2.8 MB | re-interpolate the PSF anywhere later, no rebuild |
  | `psfex_cat` | `psfex_cat-*.cat` | unmeasured | which stars PSFEx's outlier rejection clipped |
  | `star_selection` | `star_selection-*.fits` | 24.5 MB | which stars the selection cuts rejected, and why |
  | `star_train` | `star_split_ratio_80-*.fits` | 19.9 MB | the 80% sample PSFEx fitted |
  | `star_test` | `star_split_ratio_20-*.fits` | 7.1 MB | the 20% sample `psf_validation` corresponds to |
  | `star_stats` | `star_stat-*.txt` | unmeasured | setools' per-CCD counts, density and FWHM cuts |

  The default is `psf_model`. `psf_validation` is in the catalogue too but needs
  no naming; naming it anyway is harmless. **Retention is additive**: an
  existing tar is a floor, so shrinking the list adds nothing and removes
  nothing. Dropping a product is a deliberate act on `products_dir`, not a
  config edit — otherwise editing a config would delete products from the
  backed-up filesystem whose scratch originals are long gone. A raw glob is still accepted as an
  escape hatch — anything with a glob metacharacter or a dot is read as one —
  and an unknown *name* is a parse-time error listing the valid ones. The list
  is exposure-side only; tile-side retention is #844 follow-up.
- **The campaign ends in two merged catalogues, and the workflow makes both.**
  Everything above is per unit; the two products downstream analysis actually
  opens are per *campaign*, and until these rules existed each was a manual pass
  after the run.
  `star_cat_merge` collects every exposure's every CCD's `psf_validation` into
  `<products_dir>/full_starcat_<campaign>.hdf5`, one dataset per exposure at
  `exposures/<exp>` — the rho/tau statistics input. It reads the members
  straight out of the per-exposure tars (`tarfile`; unpacking ~800k files to
  merge them would defeat the tar's whole purpose), keeps their native dtypes,
  and stores `CCD_NB` as an int. sp_validation still opens the old flat FITS
  name, `full_starcat-0000000.fits`; its readers move to this file under
  [sp_validation#340](https://github.com/CosmoStat/sp_validation/issues/340),
  the same migration that retires the `patches/` key on the tile side. The rule
  exists whenever the campaign has a persisted exposure.
  **Two writers, one schema.** The module runner still emits the flat FITS
  table through `MergeStarCatPSFEX`, and this rule emits the hdf5; they are
  separate implementations on purpose, because only one of them reads tars,
  keeps native dtypes and reconciles. Their 16 COLUMN NAMES must not drift
  apart, and nothing else would notice if they did — a column added to one
  writer would just be missing from the other's product. `tests/unit/`
  `test_star_cat_columns.py` is what holds them together.
  `final_cat_merge` collects every ready tile's `final_cat-<ID>.fits` into
  `<products_dir>/final_cat_<campaign>.hdf5`: one dataset per tile under a group
  named for the campaign, the `final_cat.param` columns, an `n_tiles` attribute.
  That schema is what sp_validation's reader opens, so it is fixed; the column
  extraction reuses `scripts/python/create_final_cat.py` while the file is
  written here, because that script's own discovery walks a directory layout
  this workflow does not have. `campaign:` in `config.yaml` names the group and
  defaults to the persistent root's basename.
  BOTH RECONCILE, through one shared module (`hdf5_reconcile.py`) so the
  campaign's two products cannot disagree about what an output owes its inputs.
  Each adds the units that have no dataset, drops datasets whose unit left the
  campaign, re-reads one whose source changed (every dataset records its
  source's size and mtime) or whose column set moved (a digest on the file's
  root), and leaves the rest unread — because re-reading a campaign to add one
  unit is ~800 GB of IO at DR6 scale. The *content* is still a function of the
  input set; the byte layout is not, and a no-op leaves the file untouched
  rather than rewritten.
  Both rerun when the set changes: the unit ids' fingerprint rides on `params`.
  Neither is a `localrule` — one job over ~20k units is real work — and neither
  puts its input paths in its shell, which is not fastidiousness: ~20k paths is
  an order of magnitude over Linux's 128 KiB `MAX_ARG_STRLEN` for a single argv
  entry, so each job is handed the tile list and the run index and derives the
  same set from them.
- **Coverage is a workflow product, built from records the DAG already
  writes.** `exp_footprint` writes one JSON per exposure to
  `<products_dir>/exp/<prefix>/<base>/manifests/exp_footprint.json`, giving the
  four sky corners of every CCD that got a PSF model. It reads the valid-PSF CCD
  set off `exp_persist.json`'s tar members — exact, because `psfex_interp`
  returns *without* writing `validation_psf-*.fits` on NOT_ENOUGH_STARS,
  BAD_CHI2 or FILE_NOT_FOUND — and the WCS off `headers-<exp>.npy`, written by
  `exp_split`. It needs nothing of `persist_exp:`: those catalogues are the
  `psf_validation` product `exp_persist` packs for every exposure whatever the
  keep list says. Like `exp_persist` it runs whatever `clean:` and `coverage:`
  say, because its input is on /scratch and the purge takes it; `clean_exposure`
  takes its manifest as an input. Set
  `coverage: {enabled: true}` and one further job, `coverage_map`, stamps every
  footprint into `<products_dir>/coverage/coverage.hsp` — a HealSparse map
  counting, per sky pixel, the exposures with a valid PSF there. That job is
  **campaign-cumulative**: its declared inputs are the in-scope footprints, but
  the script reads *every* record on the products root, reclaimed exposures
  included, so appending tiles grows the map instead of replacing it. `nside` is
  set in `config.yaml` to the production 128/131072 pair, ~0.1"/pixel, chosen to
  align pixel-wise with the UNIONS bit masks; nothing defaults to it, and a
  coarser map would look plausible and not align. Plotting stays out of the DAG: run `plot_coverage_map -i
  <products_dir>/coverage/coverage.hsp ...` by hand, with the sky windows under
  `coverage.plot` in `config.yaml`. sp_validation consumes the map in
  `notebooks/demo_apply_hsp_masks.py`.
- **A dead tile can be told to stop pinning exposures.** An exposure is
  cleanable only once every consuming tile has its vignets, so one
  permanently-failed tile holds its ~80 exposures for the life of the
  campaign. List it under `clean_ignore_tiles:` in `config.yaml` and it leaves
  the consumer sets. Retrying an ignored tile later is legal and expensive:
  its exposure chains are gone and rebuild from scratch.
- **A reclaimed exposure reports as `cleaned`.** `run_report.py` reads the
  absorbed manifests out of `cleaned.json`, so a reclaimed exposure keeps its
  per-runner counts and blocks no tile. The logs go with the manifests — a log
  claiming `complete` for a store that is gone would contradict the unbuilt
  chain the DAG must now see, and its content duplicates the manifest anyway.
  The `exp_psf` benchmark tsv lives beside both dirs, not inside either, so
  reclamation does not eat the memory-sizing data.
- **Failure is a report, not a gate.** `run_report.py` disk-scans the trees
  against the count table and enumerates shortfalls (whole-unit absence vs
  per-CCD attrition). It runs standalone — a DAG report node would itself be
  poisoned by the failures it must enumerate — and fires automatically from
  the COMPUTE invocation's `onsuccess`/`onerror` hooks, or on demand via
  `sp report`.
