"""Tile post-chain — per tile: gather exposures, then detect / PSF / shape / catalogue.

    tile_exp_forest
    tile_merge_headers -> tile_detect -> tile_vignets -> tile_ngmix x N
                                                      -> tile_merge_cats -> tile_make_cat

Exposure manifests, looked up through the index, provide the tile-to-exposure
DAG edges. The per-tile exposure "forest" is a symlink view for ShapePipe's
``$SP_EXP`` globs, not the exposure dependency. See build_forest.py and
``exp_utils.get_exp_output_files`` for the required sharded layout.

The post-chain uses shell rules with no mid-chain localrules. Each connected
component within a group label is one tile, so each group job covers one tile.
See prepare.smk's docstring for Snakemake's group resource composition:

* ``group: "tile_gather"`` — tile_exp_forest (2 GB, 20 min) and
  tile_merge_headers (median 0:38; 8 GB, 4 threads, 120 min). The group asks max
  mem_mb = 8000*attempt, max threads = 4, sum runtime = 140. tile_detect also
  consumes the forest, but that edge does not pull it into the group.
* ``group: TILE_GROUP`` ("tile_shape") — tile_vignets -> 8 x tile_ngmix ->
  tile_merge_cats -> tile_make_cat in one allocation, so the vignette and WCS
  stores can stay node-local (see tile_local() below).
  First-attempt composition, verified against snakemake 9.23.1 and sbatch:
  cpus = max over levels of summed siblings = max(8, 8x1, 8, 8) = 8;
  mem_mb = max(32000, 8x5000, 16000, 16000) = 40000;
  runtime = sum of each level's maximum = 20 + 120 + 10 + 15 = 165 min.
  See tile_ngmix's resource notes for the measurements behind these budgets.
  The 165-minute total fits ``cpubase_bycore_b1``'s 180-minute ceiling.

tile_detect stays out: the shape chain does not need it co-scheduled, and
folding it in would add its runtime to a sum that has no room.

tile_detect runs SExtractor on the tile image (sextractor_runner,
config_tile_Sx.ini) for data and image simulations alike, so both get the same
detection, windowed positions, VIGNET neighbour marking and columns.
``tile_detection: unions_catalogue`` (config.yaml, the data default) adds one
rule and one step: tile_get_catalogue fetches the UNIONS per-tile catalogue
(get_images_runner, config_tile_Gic.ini), and tile_detect joins its SExtractor
rows to it before multi-epoch post-processing (MATCH_CATALOGUE, set through
SP_MATCH_CATALOGUE). See config_tile_Sx.ini and sextractor_runner's catalogue
matching implementation for the join contract and shared object identifiers.

There is no `tile_mask` rule. Tiles have no instrument flag image of their own,
so tile_detect runs SExtractor with FLAG_IMAGE = False against
default_noimaflags.param; see config_tile_Sx.ini for the detection settings.
Sky-fixed mask columns are configured in config_tile_Mc.ini and written by
make_cat for downstream selection.
"""

# --- node-local tile stores -----------------------------------------------
# The shape group writes its vignette store on one node and keeps every reader
# in that allocation. This avoids staging and ~163 GB of random NFS reads per
# tile, and saves ~5.6 GB of shared scratch per tile (~190 GiB for 34 tiles).
# Measured on nibi with remote stores: median 7h34m elapsed versus 52 min CPU,
# with ~20 GB read per chunk at 0.78 MB/s. Per-read latency dominates; the
# sequential local/shared bandwidth ratio is only ~3x.
#
# Manifests, verdict logs, chunk directories and final_cat stay on shared
# storage as DAG outputs. The node-local paths are wired through
# config_tile_PiViVi_<psf_model>.ini, config_tile_Ng_template.ini and
# config_tile_Mc.ini. tile_local() also stages the WCS sqlite (see its docstring).
# A group retry rebuilds the stores; they cannot be inherited by another job.
#
# EDITING PROLOGUE OUTPUT CAN DELETE FINISHED CATALOGUES ON RESUME.
# tile_local() contributes to params.pre. See workflow/CONTRACTS,
# unit-pre-changes-at-campaign-boundary, for the rerun hazard and pinned params.
# Change emitted prologues only at a campaign boundary on a fresh root.
# Per-campaign prefix of the node-local store name; see tile_local().
LOCAL_TAG = hashlib.sha1(str(RUN_DIR).encode()).hexdigest()[:8] + "-"


def tile_local(tile):
    """The node-local prologue, as bash, for one tile.

    Derive a literal path from the tile wildcard and LOCAL_TAG (the run-dir
    hash), so concurrent campaigns over the same tile do not share a store.
    See profiles/nibi/config.yaml's /local bind note for container path handling
    and bin/sp's tile_store_root handling for the host-side bind.
    On nibi the container's /tmp is RAM-backed, so these stores go on the
    /local/scratch bind instead.

    The WCS store is small (11.3 MB), but ngmix reads it once per object per
    epoch. With local vignettes and remote WCS, 3,200 thread-state samples
    across eight chunks on nibi show 56% CPU and 44% rpc_wait_bit_killable;
    147 of 176 D-state samples are in NFS RPC. Staging WCS avoids those
    per-read round trips, with a measured ~1.8x gain for ngmix.

    Copy WCS unconditionally to a temporary name, then rename it atomically.
    A symlink would still read NFS; cp -u could preserve a truncated file with
    a newer mtime. Concurrent group members can safely publish identical
    copies: an open reader keeps its inode. tile_merge_headers is upstream,
    so a missing source is an error rather than a reason to read another tree.

    tile_vignets clears its stale vignette run directory before writing, and
    tile_make_cat removes the local root on exit. The age-based sweep below
    removes this user's matching directories older than one day; SIGKILL can
    leave stores behind because it bypasses the exit trap.
    """
    # unit_num() supplies the leading dash and dashed tile ID. Use an exact
    # WCS filename so a missing source produces a precise error, not an empty glob.
    tile_num = unit_num(tile)
    return f'''
if [ ! -d /local/scratch ]; then
  echo "tile_shape: node-local storage unavailable." >&2
  echo "  /local/scratch is not a directory inside the container." >&2
  echo "  profiles/nibi/config.yaml apptainer-args must carry --bind /local" >&2
  exit 1
fi
export SP_LOCAL="/local/scratch/sp-{LOCAL_TAG}{tile}"
export SP_VIGNET_OUT="$SP_LOCAL/output"
export NGMIX_VIGNET_DIR="$SP_LOCAL/output/run_sp_tile_PiViVi"
export SP_WCS_DIR="$SP_LOCAL/wcs"
mkdir -p "$SP_VIGNET_OUT" "$SP_WCS_DIR" || {{
  echo "tile_shape: cannot create $SP_LOCAL on this node." >&2
  exit 1
}}
# the WCS store -- see the docstring. Copy to a temp name and rename, because
# `cp` is not atomic and the eight chunks share this destination.
sp_wcs_src="$SP_RUN/output/run_sp_tile_Mh_exp/merge_headers_runner/output/log_exp_headers{tile_num}.sqlite"
cp "$sp_wcs_src" "$SP_WCS_DIR/.log_exp_headers.$$" && \
mv -f "$SP_WCS_DIR/.log_exp_headers.$$" \
      "$SP_WCS_DIR/log_exp_headers{tile_num}.sqlite" || {{
  rm -f "$SP_WCS_DIR/.log_exp_headers.$$"
  echo "tile_shape: could not stage the WCS store." >&2
  echo "  source: $sp_wcs_src" >&2
  echo "  dest:   $SP_WCS_DIR" >&2
  echo "  check that tile_merge_headers ran, and that /local/scratch has room." >&2
  exit 1
}}
find /local/scratch -maxdepth 1 -name 'sp-*.*' -user "$(id -u)" -mmin +1440 \
     -exec rm -rf {{}} + 2>/dev/null || true
'''


# All group members must use the same string resource (Snakemake's
# resources.py::_is_string_resource), so each uses this constant.
# Slurm's --tmp selects nodes by configured TmpDisk; it does not reserve or
# decrement disk space per job. 16 GB leaves room above the measured 5.6 GB
# store and excludes nodes configured without sufficient temporary disk.
TILE_SLURM_EXTRA = "--tmp=16000"

# One label for the whole shape chain; see this file's docstring for resource
# composition. The summed runtime determines eligible partitions on nibi.
TILE_GROUP = "tile_shape"

# Set cleanup only on tile_make_cat, the last store reader. Earlier members'
# shells exit while other members still need the store. See tile_local() for
# cleanup after a hard kill.
TILE_CLEAN = r"""
trap 'rm -rf "$SP_LOCAL"' EXIT
"""

# Only tile_vignets may clear the vignette run directory before writing.
# Putting this in tile_local() would make readers delete their own input.
TILE_VIGNET_FRESH = r"""
rm -rf "$SP_VIGNET_OUT/run_sp_tile_PiViVi"
"""

# tile_vignets writes the partition once, before all chunks; each chunk reads
# its row rather than recomputing boundaries. A measured mid-run splitter edit
# with independent chunk splitting orphaned 520 objects and measured 831 twice;
# merge_sep_cats would accept that inconsistent partition without error.
# The input sexcat is on the shared root; the ranges are group-local plumbing,
# not a rule output. See ngmix_range.py for the partition algorithm.
NGMIX_RANGES = "$SP_LOCAL/ngmix_ranges.json"
TILE_NGMIX_RANGES = (
    f'python {SCRIPTS}/ngmix_range.py --run-dir "$SP_RUN" '
    f'--n-chunks {NGMIX_CHUNKS} --write "{NGMIX_RANGES}" || exit 1'
)

# Snakemake tracks the manifest, not the node-local store. After a partial
# group failure, a resume can omit tile_vignets because its manifest exists,
# leaving readers without a store. Delete tile_vignets.json and resume to put
# the writer back in the group. An in-flight retry keeps the original member set.
# Keep the manifest non-temp: clean_exposure uses it to establish eligibility.
TILE_VIGNET_REQUIRED = r"""
if [ ! -d "$NGMIX_VIGNET_DIR/vignetmaker_runner_run_2/output" ]; then
  echo "tile_shape: the node-local vignette store is missing." >&2
  echo "  expected: $NGMIX_VIGNET_DIR" >&2
  echo "  This means tile_vignets did NOT run in this group job — its manifest" >&2
  echo "  was already satisfied, most likely because a previous fused job for" >&2
  echo "  this tile died after tile_vignets succeeded. The store lives and dies" >&2
  echo "  with the job, so it cannot be inherited." >&2
  echo "  FIX: rm \"\$SP_RUN/manifests/tile_vignets.json\" and resume." >&2
  exit 1
fi
"""




# Preserve run and module logs before temporary or node-local directories go.
# This includes ngmix's epoch-cut counts and per-object failure messages.
# Each rule copies logs to the tile's shared logs/modules/<run name>/ even on
# failure; a copy failure does not change the rule's exit status. clean_tile
# reclaims these logs with the rest of the tile store.
def keep_module_logs(run_root, run_name):
    dest = f'"$SP_RUN/logs/modules/{run_name}"'
    return (
        f'rm -rf {dest} && mkdir -p {dest} && '
        f'(cd "{run_root}" && find . -path "*/logs/*" -type f '
        f'-exec cp -p --parents -t {dest} {{{{}}}} +) || '
        f'echo "could not keep the module logs of {run_root}" >&2\n'
    )


def tile_exp(wc):
    return tile_exposures(wc.tile)

# --- tile-to-exposure edges after reclamation -----------------------------
# A tile with final_cat needs no reclaimed exposure store. Drop only missing
# exposure manifests for finished tiles, so rebuilding a shared exposure for
# one tile does not propagate reruns through the overlap component. Unfinished
# tiles keep every edge and rebuild the exposures they need. ancient() suppresses
# timestamp changes, but cannot suppress reruns propagated from upstream jobs.
#
# The input rerun trigger must be off (see the profiles' trigger lists): it
# treats the deliberate edge removal as a reason to rerun. On a four-tile
# fixture with one damaged tile, job counts are 82 without these protections,
# 70 with edge removal but the input trigger on, and 28 with both protections.
#
# Use final_cat as the marker, not tile_vignets.json: the manifest can survive
# without the node-local store. To rebuild reclaimed exposures with --forcerun,
# delete the tile's final_cat first so its exposure edges return to the DAG.
# See clean_exposure in exposure.smk for reclamation and consumer-set tracking.
def tile_finished(tile):
    return Path(final_cat(tile)).exists()


def exp_manifests(wc, stage):
    paths = [exp_manifest(e, stage) for e in tile_exp(wc)]
    if tile_finished(wc.tile):
        paths = [p for p in paths if Path(p).exists()]
    return [ancient(p) for p in paths]

def tile_exp_split(wc):  return exp_manifests(wc, "exp_split")
def tile_exp_psf(wc):    return exp_manifests(wc, "exp_psf")
def tile_exp_all(wc):    return tile_exp_split(wc) + tile_exp_psf(wc)


# Build the per-tile symlink forest. Declaring the exposure manifests as input
# makes this wait on its exposures; the forest itself is only the $SP_EXP view.
# Its output is a directory() (it has no ShapePipe run dir and no manifest —
# it is not a shapepipe_run at all).
rule tile_exp_forest:
    group: "tile_gather"
    input:
        tile_exp_all
    output:
        forest = directory(f"{TILE_DIR}/exp_forest")
    params:
        cmd = lambda wc: (f"python {SCRIPTS}/build_forest.py --tile {wc.tile} "
                          f"--run-dir {RUN_DIR} --index {INDEX_DB}"),
        # build_forest.py's own content hash rides here and nowhere else.
        script_hash = FOREST_HASH
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 20
    shell:
        # --forest {output} lives in the shell string: snakemake formats shell
        # once, so an {output} placeholder inside params.cmd would survive
        # literally and every forest job would race one './{output}'.
        "{params.cmd} --forest {output.forest}"

# Merge single-exposure WCS headers into the tile-level sqlite
# (log_exp_headers-<IDra>-<IDdec>.sqlite, which Sx / PiViVi / ngmix consume).
# Reads headers-*.npy through the forest -> the split manifests are the edge.
rule tile_merge_headers:
    group: "tile_gather"
    input:
        forest = rules.tile_exp_forest.output.forest,
        split  = tile_exp_split,
        # config_tile_Mh_exp.ini reads the prepare phase's run_sp_tile_Fe output.
        # Its manifest provides a regeneration edge in the compute DAG.
        fe     = f"{TILE_DIR}/manifests/tile_find_exposures.json",
    output:
        manifest = f"{TILE_DIR}/manifests/tile_merge_headers.json"
    log:
        f"{TILE_DIR}/logs/tile_merge_headers.json"
    params:
        pre = lambda wc: unit_pre("tile_merge_headers", wc.tile,
                                  forest=forest_dir(wc.tile)),
        script_hash = SCRIPT_HASH
    threads: 4
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 120
    shell:
        sp_shell("tile_merge_headers", "config_tile_Mh_exp.ini")

# The UNIONS per-tile catalogue (CFIS.<tile>.r.cat), from a local mirror or
# from vos, the way tile_get_images fetches the image; tile_detect joins its
# SExtractor detections to it. Reads only tile_numbers.txt, which unit_pre
# writes; the edge on the image manifest is what puts it after the prepare
# phase.
if TILE_DETECTION == "unions_catalogue":

    rule tile_get_catalogue:
        input:
            git = f"{TILE_DIR}/manifests/tile_get_images.json",
        output:
            manifest = f"{TILE_DIR}/manifests/tile_get_catalogue.json"
        log:
            f"{TILE_DIR}/logs/tile_get_catalogue.json"
        params:
            pre = lambda wc: unit_pre(
                "tile_get_catalogue", wc.tile,
                env={"SP_INPUT_CATALOGUES": CATALOGUES,
                     "SP_RETRIEVE_CATALOGUES": CATALOGUE_RETRIEVE}),
            script_hash = SCRIPT_HASH
        threads: 1
        retries: 2
        resources:
            mem_mb = lambda wc, attempt: 4000 * attempt,
            runtime = 60
        shell:
            sp_shell("tile_get_catalogue", "config_tile_Gic.ini")


# tile_detect's inputs: the catalogue manifest too when it joins one.
DETECT_INPUTS = {"uz": f"{TILE_DIR}/manifests/tile_uncompress.json",
                 "mh": rules.tile_merge_headers.output.manifest}
if TILE_DETECTION == "unions_catalogue":
    DETECT_INPUTS["cat"] = rules.tile_get_catalogue.output.manifest


def detect_env(tile):
    """tile_detect's prologue exports: the catalogue to join, if any.

    Empty under tile_detection: sextractor (image simulations), so that no
    value left exported in the submitting shell reaches the job.
    """
    if TILE_DETECTION != "unions_catalogue":
        return {"SP_MATCH_CATALOGUE": ""}
    gic = f"{tile_dir(tile)}/output/run_sp_tile_Gic/get_images_runner/output"
    return {"SP_MATCH_CATALOGUE": f"{gic}/CFIS_cat{unit_num(tile)}.cat"}


# SExtractor object detection on the tile; under unions_catalogue joined to
# the UNIONS catalogue (config_tile_Sx.ini, MATCH_CATALOGUE).
rule tile_detect:
    input:
        **DETECT_INPUTS
    output:
        manifest = f"{TILE_DIR}/manifests/tile_detect.json"
    log:
        f"{TILE_DIR}/logs/tile_detect.json"
    params:
        pre = lambda wc: unit_pre("tile_detect", wc.tile,
                                  env=detect_env(wc.tile)),
        script_hash = SCRIPT_HASH
    # One core: the container's SExtractor has no threading support, and the
    # join and post-processing are serial. See default_tile.sex for NTHREADS.
    # Across eight DR6 tiles on candide (1-10 exposures, including low latitude,
    # a cluster, a bright-star halo and survey edges), wall time is 23-160 s
    # and peak RSS is 1.84 GiB. 4000 MB gives ~2.1x memory margin; 20 min gives
    # ~7x wall-time margin for shared-storage I/O on nibi. Both scale on retry.
    threads: 1
    retries: 1
    resources:
        mem_mb = lambda wc, attempt: 4000 * attempt,
        runtime = lambda wc, attempt: 20 * attempt
    shell:
        sp_shell("tile_detect", "config_tile_Sx.ini")

# PSF interpolation and galaxy postage stamps: the last exposure-product reader.
# The vignette store is node-local (see tile_local()).
rule tile_vignets:
    group: TILE_GROUP
    input:
        sx     = rules.tile_detect.output.manifest,
        forest = rules.tile_exp_forest.output.forest,
        split  = tile_exp_split,
        psf    = tile_exp_psf,
        # config_tile_PiViVi_<psf_model>.ini reads run_sp_tile_Fe output — same reason as
        # tile_merge_headers above.
        fe     = f"{TILE_DIR}/manifests/tile_find_exposures.json",
    output:
        # The manifest is the only declared output. The store is node-local,
        # all readers are in this group, and TILE_CLEAN reclaims it regardless
        # of --notemp.
        manifest = f"{TILE_DIR}/manifests/tile_vignets.json",
    log:
        f"{TILE_DIR}/logs/tile_vignets.json"
    params:
        pre = lambda wc: unit_pre("tile_vignets", wc.tile,
                                  forest=forest_dir(wc.tile),
                                  pre_run=[tile_local(wc.tile), TILE_VIGNET_FRESH,
                                           TILE_NGMIX_RANGES]),
        script_hash = SCRIPT_HASH,
        # Fingerprint the partition writer and readers together; see the
        # tile_ngmix range_hash note and the prologue hazard above.
        range_hash = NGMIX_RANGE_HASH
    # Eight threads match the ngmix wave's width, without making this stage
    # request a wider group allocation. ShapePipe's batch parallelism is over
    # input file sets; a tile supplies one set (see tile_ngmix's thread note).
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 32000 * attempt,
        # Measured median 3m29s, max 5:25. The 20-minute budget includes
        # margin and contributes to the group's summed runtime.
        runtime = 20,
        slurm_extra = TILE_SLURM_EXTRA
    shell:
        # The completeness check is pointed at the node-local run root; see
        # sp_shell's check_args for what the two flags do.
        sp_shell("tile_vignets", f"config_tile_PiViVi_{PSF_MODEL}.ini",
                 check_args=' --run-dir "$SP_LOCAL" --unit {wildcards.tile}',
                 post=keep_module_logs("$NGMIX_VIGNET_DIR", "run_sp_tile_PiViVi"))

# ngmix shape measurement: N chunks per tile, each reading its closed row range
# from TILE_NGMIX_RANGES. The tile's sexcat supplies the bounds at execution time,
# so a params function cannot determine them. ngmix treats ID_OBJ_MAX <= 0 as
# unbounded; a closed upper bound prevents the last chunk from overrunning.
#
# Chunks write separate directories: each has its own run_sp_tile_ngmix_Ng<k>u, and
# merge_sep_cats — DAG-serialised after all chunks — is the gather.
rule tile_ngmix:
    group: TILE_GROUP
    input:
        # The manifest is the DAG edge; see tile_vignets for store lifetime.
        vignets  = rules.tile_vignets.output.manifest,
        sx       = rules.tile_detect.output.manifest,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_ngmix_{{chunk}}.json",
        # Shared DAG output (~300 KB/chunk), reclaimed after tile_merge_cats.
        # merge_sep_cats derives chunks 2..N from chunk 1's path.
        chunkdir = temp(directory(f"{TILE_DIR}/output/run_sp_tile_ngmix_Ng{{chunk}}u")),
    log:
        f"{TILE_DIR}/logs/tile_ngmix_{{chunk}}.json"
    params:
        pre = lambda wc: unit_pre("tile_ngmix", wc.tile,
            env={"SP_NGMIX_CHUNK": wc.chunk, "NGMIX_N_CHUNKS": NGMIX_CHUNKS},
            # tile_local exports SP_LOCAL before the range lookup. Capture and
            # check the lookup before eval: eval "$(...)" would discard its
            # failure status and let ShapePipe run with unset bounds.
            pre_run=[tile_local(wc.tile), TILE_VIGNET_REQUIRED,
                     f'ngmix_range_out=$(python {SCRIPTS}/ngmix_range.py '
                     f'--read "{NGMIX_RANGES}" --chunk {wc.chunk}) || exit 1',
                     'eval "$ngmix_range_out"']),
        script_hash = SCRIPT_HASH,
        # Hash the script body, not just its invocation, for resume reproducibility.
        # Keep this fingerprint on tile_vignets too: a readers-only params change
        # can omit the writer, trip TILE_VIGNET_REQUIRED, and make failed-group
        # cleanup delete an existing final_cat (observed in a nibi group job).
        # See the prologue hazard above and workflow/CONTRACTS,
        # unit-pre-changes-at-campaign-boundary, before changing params.
        range_hash = NGMIX_RANGE_HASH
    # One core per chunk: ShapePipe's -b controls parallelism over input file
    # sets (pipeline/args.py), and a chunk is one catalogue. ngmix has no internal
    # worker pool; the prologue and profile pin OMP/BLAS to one thread.
    threads: 1
    retries: 2
    benchmark:
        f"{TILE_DIR}/manifests/tile_ngmix_{{chunk}}.benchmark.tsv"
    resources:
        # Across 31 nibi tiles, cgroup high-water is 27.40 GiB max, 22.14 GiB
        # mean. Eight sibling chunks sum to 40000 MiB = 39.1 GiB, a 1.43x
        # margin over the measured maximum. nibi's cgroup-based sacct MaxRSS
        # includes reclaimable page cache, not just process memory.
        # Prefer the larger per-worker memory estimate in config_tile_Ng_template.ini's
        # SAVE_BATCH note over the ~1.25 GiB psutil benchmark: 30-second sampling
        # can miss flush peaks. That leaves about 20 GiB for the local store's cache.
        #
        # nibi bills max(cores, mem_GB/4): 40 GB costs 10 core-equivalents for
        # eight cores. Reducing chunks to 4000 would make tile_vignets' 32000
        # the group limit, only 1.14x the measured maximum; measure cache
        # behaviour under pressure before using that budget. Retries scale memory.
        mem_mb = lambda wc, attempt: 5000 * attempt,
        # Measured on nibi: 1.29 s/object for one local chunk, 1.54 s/object
        # for eight concurrent chunks (~19% slower). At 4678 objects/chunk,
        # the latter projects to ~120 min. Slurm enforces the group's total,
        # not this member budget: attempt 1 gets 165 min, attempt 2 gets 285.
        # Scale retries so a timeout receives more wall time rather than
        # repeating the same allocation; attempt 1 fits cpubase_bycore_b1.
        runtime = lambda wc, attempt: 120 * attempt,
        slurm_extra = TILE_SLURM_EXTRA
    shell:
        sp_shell("tile_ngmix", "config_tile_Ng_template.ini",
                 post=keep_module_logs("{output.chunkdir}",
                                       "run_sp_tile_ngmix_Ng{wildcards.chunk}u"))


def ngmix_manifests(wc):
    return [f"{tile_dir(wc.tile)}/manifests/tile_ngmix_{k}.json"
            for k in range(1, NGMIX_CHUNKS + 1)]

def ngmix_chunkdirs(wc):
    return [f"{tile_dir(wc.tile)}/output/run_sp_tile_ngmix_Ng{k}u"
            for k in range(1, NGMIX_CHUNKS + 1)]

# The gather: merge the N chunk catalogues. N_SPLIT_MAX comes from the workflow's
# own chunk count via $NGMIX_N_CHUNKS (env-expanded by the module).
rule tile_merge_cats:
    group: TILE_GROUP
    input:
        manifests = ngmix_manifests,
        chunkdirs = ngmix_chunkdirs,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_merge_cats.json"
    log:
        f"{TILE_DIR}/logs/tile_merge_cats.json"
    params:
        pre = lambda wc: unit_pre("tile_merge_cats", wc.tile,
                                  env={"NGMIX_N_CHUNKS": NGMIX_CHUNKS}),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        # Measured median 0:15, max 0:42. A term in the group's runtime sum.
        runtime = 10,
        slurm_extra = TILE_SLURM_EXTRA
    shell:
        sp_shell("tile_merge_cats", "config_tile_Ms.ini")

# The run's science product. make_cat also reads the vignette store's
# configured PSF-interpolation output, so it — not ngmix — is the store's last reader.
#
# No protected(): the profile's rerun triggers govern catalogue regeneration.
rule tile_make_cat:
    group: TILE_GROUP
    input:
        # See tile_vignets for store lifetime, and config_tile_Mc.ini for
        # make_cat's PSF-interpolation input through NGMIX_VIGNET_DIR.
        ms    = rules.tile_merge_cats.output.manifest,
    output:
        manifest  = f"{TILE_DIR}/manifests/tile_make_cat.json",
        final_cat = f"{PROD_TILE_DIR}/final_cat-{{tile}}.hdf5",
    log:
        f"{TILE_DIR}/logs/tile_make_cat.json"
    params:
        pre = lambda wc: unit_pre("tile_make_cat", wc.tile,
                                  pre_run=[tile_local(wc.tile), TILE_VIGNET_REQUIRED,
                                           TILE_CLEAN]),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        # Measured median 1:33, max 1:50. A term in the group's runtime sum.
        runtime = 15,
        slurm_extra = TILE_SLURM_EXTRA
    shell:
        # The plain body plus the catalogue publish: final_cat is a real file, so
        # it is a real declared output (and it persists — never temp()). The
        # publish is guarded on rc, so a job whose shapepipe_run died never
        # publishes a catalogue for a manifest snakemake is about to delete.
        sp_shell("tile_make_cat", "config_tile_Mc.ini",
                 post="if [ $rc -eq 0 ]; then\n"
                      '  cp -f "$(ls -1 "$SP_RUN"/output/run_sp_tile_Mc/make_cat_runner'
                      '/output/final_cat*.hdf5 | head -1)" {output.final_cat}\n'
                      "fi\n")


# --- tile reclamation -----------------------------------------------------
# Reclaim each ready tile after its persistent final_cat lands, keeping shared
# scratch proportional to tiles in flight. The input is the tile-finished marker
# defined by final_cat() in the Snakefile. See clean_tile_targets() there for
# scope, and clean_tile.py for the measured footprint and survivor set.
#
# Rebuilding a cleaned tile must also rebuild its node-local stores. Deleting
# its upstream outputs makes their reruns propagate to tile_vignets, the store's
# writer. Preserving enough upstream outputs to leave tile_vignets up to date
# would instead schedule its readers without it. See clean_tile.py before
# changing the survivor set.
# This rule is a local DAG leaf (declared in the Snakefile), outside all groups.
rule clean_tile:
    input:
        # The tile-finished marker, on the persistent root. A lambda rather than
        # the PROD_TILE_DIR pattern so there is exactly one definition of this
        # path (final_cat() in the Snakefile), the same way clean_exposure keys
        # off wildcards.exp alone.
        lambda wc: final_cat(wc.tile)
    output:
        tombstone = f"{TILE_DIR}/cleaned.json"
    params:
        # clean_tile.py is external to the shell string, so the `code`
        # rerun-trigger cannot see it — same reason SCRIPT_HASH exists.
        script_hash = CLEAN_TILE_HASH
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    shell:
        f"python {SCRIPTS}/clean_tile.py"
        " --tile-dir $(dirname {output.tombstone}) --tile {wildcards.tile}"
        " --tombstone {output.tombstone}"


# --- campaign shear catalogue --------------------------------------------
# One job merges the same ready-tile final_cat list requested by rule all.
# Pass the tile list and index, not every path (MAX_ARG_STRLEN), and fingerprint
# that input set so an append triggers reconciliation. Globbing products_dir
# could include tiles outside the campaign. See merge_final_cat.py for the
# reconciliation algorithm and output schema, and workflow/CONTRACTS for
# campaign-name-is-run and final-cat-param-is-exact-allow-list.
#
# Run on a compute node: the first build reads ~32-46 MB per tile (~800 GB for
# DR6's 23k tiles). Memory holds one tile plus the HDF5 buffer; runtime budgets
# the full first build. Reconciliation reads only new or changed source catalogues.
rule final_cat_merge:
    input:
        lambda wc: [final_cat(t) for t in TILES_READY]
    output:
        merged = final_cat_hdf5()
    params:
        products_dir = str(PRODUCTS_DIR),
        tile_list    = str(config["tile_list"]),
        index_db     = str(INDEX_DB),
        param_file   = str(CONFIG_DIR / "final_cat.param"),
        campaign     = CAMPAIGN,
        snapshot     = str(SNAPSHOT_JSON),
        inputs       = unit_fingerprint(TILES_READY),
        script_hash  = MERGE_FINAL_HASH
    threads: 1
    resources:
        # See the Snakefile's sizing block for the largest-tile memory model.
        mem_mb = lambda wc, attempt: capped_mem(attempt * (
            FINAL_MEM_BASE_MB
            + FINAL_MEM_FACTOR * final_cat_max_bytes() // 1_000_000),
            "final_cat_merge"),
        # First-build budget: ~1 min per 10 tiles measured, with 3x margin
        # above a fixed floor. Reconciliation can read fewer tiles.
        runtime = lambda wc, attempt: attempt * (30 + len(TILES_READY) // 3)
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/merge_final_cat.py"
        " --products-dir '{params.products_dir}'"
        " --tile-list '{params.tile_list}' --index-db '{params.index_db}'"
        " --output {output.merged}"
        " --campaign '{params.campaign}'"
        " --param-file '{params.param_file}'"
        " --snapshot-json '{params.snapshot}'"
