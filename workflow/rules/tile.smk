"""Tile post-chain — per tile: gather exposures, then detect / PSF / shape / catalogue.

    tile_exp_forest
    tile_merge_headers -> tile_detect -> tile_vignets -> tile_ngmix x N
                                                      -> tile_merge_cats -> tile_make_cat

The DAG edge to the exposures is always the exposures' MANIFESTS, looked up
through the index (TILE_EXP). The per-tile exposure "forest" — a symlink view of
exactly this tile's exposures' products — exists only so the ShapePipe configs
have one deterministic ``$SP_EXP`` to glob; it is NEVER the edge. Its 2-char
shard level is not cosmetic: ``exp_utils.get_exp_output_files`` hardwires
``<SP_EXP>/<prefix>/<base>/output/run_sp_*`` into its glob, so a flat forest
makes every tile gather stage fail "No split_exp_runner output found".

All rules are group-compatible (shell only, no mid-chain localrules), and the two
short regions are grouped, per the composition rules in prepare.smk's docstring.
Distinct tiles share no edge, so each is one group job per tile:

* ``group: "tile_gather"`` — tile_exp_forest (2 GB, 20 min) and
  tile_merge_headers (median 0:38; 8 GB, 4 threads, 120 min). The group asks max
  mem_mb = 8000*attempt, max threads = 4, sum runtime = 140. tile_detect also
  consumes the forest, but a consumer OUTSIDE the group is just an ordinary DAG
  edge on the group job — it does not pull tile_detect (16 GB) in.
* ``group: TILE_GROUP`` ("tile_shape") — tile_vignets -> 8 x tile_ngmix ->
  tile_merge_cats -> tile_make_cat, THE WHOLE SHAPE CHAIN AS ONE JOB. This is
  not a latency optimisation; it is what lets the 5.6 GB vignette store live on
  node-local NVMe and never touch /scratch (see tile_local() below).
  Composition, verified against snakemake 9.23.1 and a real sbatch:
  cpus = max over levels of summed siblings = max(8, 8x1, 8, 8) = 8;
  mem_mb = max(32000, 8x5000, 16000, 16000) = 40000, so nibi's
  max(cores, mem_GB/4) bills the group 10 core-equivalents rather than the 28
  it billed at 8x14000 (see the mem_mb note on tile_ngmix for the 31-tile
  measurement that sized it); tile_vignets' own 32000 is now the second term
  and becomes binding if the chunks ever go below 4000;
  runtime = sum along the chain of each level's MAX = 20 + 120 + 10 + 15 = 165.
  165 min is under the 180 min ceiling of ``cpubase_bycore_b1``, so the fused
  job reaches the widest partition set (plus cpubackfill) — which is why each
  member's runtime is measured p99 plus margin and not the old ceiling.

The heavy middle (tile_detect) stays out: it is a 16 GB / 8 thread SExtractor
run that the shape chain does not need co-scheduled, and folding it in would add
its runtime to a sum that has no room.

Note there is no `tile_mask` rule: the committed config chain is the
"sx_nomask" tile_detect variant (config_tile_Sx.ini reads Git + Uz + Mh, no mask
run), and no tile-mask config was committed in the S2 sweep. Adding the masked
variant is a config + one rule, at the config selector the PRD describes.

That rule also needs a tile-side analogue of ``exp_star_cat``: tile star cats key
on TILE id, so they are a separate cache namespace and a separate node, and the
earliest point it can run is after ``tile_uncompress`` (create_star_cat.py's
``-k tile`` mode reads the uncompressed tile image's primary header). The mask
config would then read a real per-unit ``$SP_RUN/star_cat_tiles`` directory,
built the same way and for the same reason (the file handler intersects numbers
across INPUT_DIRs, so a shared pool cannot be symlinked in wholesale).
"""

# --- the node-local tile root (the I/O + scratch fix) ----------------------
#
# MEASURED PROBLEM. tile_ngmix was I/O-bound, not CPU-bound: median 7h34m
# elapsed against 52 min of CPU (3% efficiency, 0.12 of 4 reserved cores).
# Each chunk read ~20 GB at 0.78 MB/s and all 8 chunks of a tile re-read the
# SAME 5.6 GB vignette store -- ~163 GB of small random sqlite preads per tile
# against /scratch, which on nibi is NFS (VAST, 4.5 PB, 95% full). What hurts
# is per-read LATENCY, not bandwidth: the sequential local-vs-scratch gap is
# only ~3x and is NOT the argument for this change.
#
# THE FIX IS THE FUSE. tile_vignets -> 8 x tile_ngmix -> tile_merge_cats ->
# tile_make_cat run as ONE group job on ONE node, and the vignette store is
# WRITTEN to that node's NVMe and never lands on /scratch at all. So the store
# is not copied, it is simply never remote: zero staging cost, zero NFS random
# reads, and -- the reason this was chosen over per-chunk staging -- the
# per-tile scratch high-water drops by the whole 5.6 GB store (~190 GiB across
# 34 concurrent tiles). Scratch quota, not speed, is what caps batch size.
#
# WHAT STAYS ON SHARED STORAGE, and it is everything that is DAG currency:
# every manifest (the success sentinels), every verdict log, the ngmix chunk
# dirs merge_sep_cats gathers, and final_cat on the PERSISTENT root. Only the
# bulk intra-tile intermediate is node-local. That split is possible because
# ShapePipe's configs set input and output paths independently -- see
# config_tile_PiViVi.ini (OUTPUT_DIR = $SP_VIGNET_OUT) and
# config_tile_Ng_template.ini and config_tile_Mc.ini ($NGMIX_VIGNET_DIR), and
# config_tile_Ng_template.ini alone for $SP_WCS_DIR.
#
# TWO node-local stores, not one, and the second is small on purpose. The
# vignette store is ~5-8 GB of bulk; $SP_WCS_DIR holds a single 11.3 MB file,
# log_exp_headers-<tile>.sqlite, which ngmix reads once per object per epoch.
# Moving 5.6 GB bought 4.3x and left 44% of every chunk blocked in NFS RPC on
# that one file. What costs on a network filesystem is the operation, not the
# byte, and the two stores are the same fix applied at the two ends of that
# distribution.
#
# WHY IT NEEDS THE GROUP. The store is node-local, so it exists only on the
# machine that wrote it and only while that job lives. Unfused, tile_vignets and
# tile_ngmix are different jobs on (usually) different nodes, and the store would
# have to be on shared storage -- which is the whole problem. Fused, all members
# run in one allocation on one node, so the store the first member writes is
# simply there for the rest. (Verified live: group members share the directory.)
#
# HOW THE PATH GETS IN HERE. Not through the environment. profiles/nibi passes
# --bind /local, and tile_local() below builds the path from the tile wildcard,
# which snakemake substitutes at DAG time. $SLURM_TMPDIR was the obvious choice
# and does not work: snakemake escapes the `$` when it splices apptainer-args
# into the job command, so the container gets the literal string and every fused
# job dies on the guard (campaign 20798193, cancelled after 8 minutes). The guard
# stays, because its job now is to catch a missing --bind rather than a missing
# variable -- and either way the failure must be loud, never a silent fall back
# onto NFS or into /tmp, which in this container is a tmpfs charged to the job's
# memory cgroup.
#
# THE COST WE ACCEPT: a failure anywhere in the tile re-runs the WHOLE tile,
# not one chunk, because the store dies with the job. At ~1 h per fused tile
# that is a cheap trade for the scratch it buys. Noted, not engineered around.
#
# It lives in each rule's `pre_run`, i.e. in THAT RULE's params.pre. The shared
# prologue (unit_pre) and completeness.py are untouched by design: both are
# fingerprinted into every rule, and changing either would invalidate the whole
# campaign, including the 117 finished exposure chains.
#
# EDITING THIS FUNCTION INVALIDATES FINISHED TILES. READ THIS BEFORE YOU DO.
# The returned string lands in `params.pre`, `params` is an active rerun trigger
# in profiles/nibi/config.yaml, and `params.pre` is a non-derived param (its
# lambda takes only `wc`), so snakemake records it and compares it. Any edit
# here therefore reschedules tile_vignets, all eight tile_ngmix chunks and
# tile_make_cat for EVERY tile the campaign already finished.
#
# On a fresh campaign that is free. On a RESUME it is destructive, and the
# mechanism is worth spelling out because it is not obvious:
#   * clean_exposure has already reclaimed the exposure stores those tiles read,
#     and tile_finished() deliberately drops the manifests of a finished tile's
#     exposures, so the rerun is UNSATISFIABLE -- PiViVi runs against dangling
#     symlinks and the group fails;
#   * a failed group job's postprocess(error=True) fans out over every member in
#     every toposort level and removes each member's EXISTING outputs. One of
#     those is tile_make_cat's final_cat on the persistent root. So the failure
#     deletes the science product of a tile that was finished and correct;
#   * with final_cat gone, tile_finished() flips and the tile re-declares its
#     whole exposure edge set, rebuilding the reclaimed chains from VOS. That is
#     the rerun avalanche D5 exists to prevent, arriving through the params
#     trigger instead of through inputs.
# So: land changes to this function BETWEEN campaigns, never mid-campaign. If a
# resume is genuinely needed after an edit, run it once with
# `--rerun-triggers mtime code software-env`, accepting that clean_exposure's
# consumer-set staleness detection (which rides on params) is off for that
# invocation.
def tile_local(tile):
    """The node-local prologue, as bash, for one tile.

    A FUNCTION of the tile rather than a constant, because the path has to be
    literal by the time apptainer sees it. `$SLURM_TMPDIR` cannot be used:
    snakemake escapes the `$` when it splices `apptainer-args` into the job
    command, so the container receives the eight characters `$SLURM_TMPDIR`
    and every fused job dies (observed live, campaign 20798193). Passing it
    through `--env` fails the same way, and `APPTAINERENV_*` would have to be
    exported by the job wrapper, which is snakemake's and not ours.

    So the path is derived instead of communicated. `/local/scratch` on nibi is
    `drwxrwxrwt root:root` -- world-writable with the sticky bit, verified in
    the container (probe job 20798618) -- and the tile id is a wildcard
    snakemake substitutes at DAG time, so the shell string carries a concrete
    path with no `$` left for anything to escape. One group job per tile means
    the name cannot collide; the sticky bit means nobody else can remove it.

    What we give up is Slurm's own cleanup of `$SLURM_TMPDIR`. TILE_VIGNET_FRESH
    reclaims a stale directory on the next attempt for the same tile, the trap
    in the last group member removes it on the way out, and the sweep below
    catches what a hard kill leaves behind -- on a shared node that last one is
    manners, not housekeeping.

    NOT `/tmp`: inside the container that is a 378 GB tmpfs, i.e. RAM charged to
    the job's cgroup, so a 5.6 GB store there would be paid for twice.

    TWO STORES, AND THE SECOND IS 11 MB. `$SP_WCS_DIR` holds one file,
    `log_exp_headers-<tile>.sqlite`, and staging it is the rest of the same fix
    rather than a refinement of it. Moving the 5.6 GB vignette store off NFS
    bought 4.3x and did not solve the problem: thread-state sampling of the
    fused job (3,200 samples across 8 chunks, campaign smk-g4 job 20799387) put
    every chunk at 56% on CPU and 44% in `rpc_wait_bit_killable`, with 147 of
    176 D-state samples inside NFS RPC and exactly one hot NFS file still open.
    ngmix reads that file once per object PER EPOCH for the WCS, so what costs
    is the operation count and not the byte count -- a network round trip per
    read. The two stores are the same intervention at the two ends of the size
    distribution, and the small one is worth ~1.8x on the rule that is 99.3% of
    tile-side core-hours.

    A copy, never a symlink: a link resolves straight back to NFS. The prologue
    runs in `tile_vignets` and in each of the eight `tile_ngmix` chunks (not in
    `tile_merge_cats`, which has no pre_run), so the staging is attempted nine
    times per tile and eight of those are concurrent siblings in one toposort
    level. It is written as copy-to-a-temp-name plus `mv -f`, which is atomic
    within a filesystem, and it is UNCONDITIONAL -- no `cp -u`, no
    already-there test. `cp` is not atomic, so an interrupted copy leaves a
    truncated destination whose mtime is NEWER than the source; `-u` would then
    skip it forever, and nothing downstream would catch it, because
    TILE_VIGNET_FRESH and TILE_VIGNET_REQUIRED both look only at the vignette
    store. A silently truncated WCS store reaches ngmix as wrong astrometry
    rather than as an error, which is the worst failure mode available here.
    Nine unconditional copies of 11 MB per tile, once at job start, is not a
    cost worth reasoning about; a staleness rule would be. Renaming over a file
    a sibling chunk already has open is safe -- the open descriptor keeps the
    old inode, whose contents are identical.

    `tile_merge_headers` is upstream of the whole group, so the file exists by
    the time any member runs; failing loudly if it does not is correct, because
    the alternative is ngmix silently reading a different tree.
    """
    # `log_exp_headers-<IDra>-<IDdec>.sqlite`, named from the tile in ShapePipe's
    # dashed image-number form -- spelled out rather than globbed so a missing
    # file is one precise error instead of an empty expansion. unit_num() owns
    # that convention and supplies the SEPARATING DASH itself, so there is none
    # in the literals below.
    tile_num = unit_num(tile)
    return f'''
if [ ! -d /local/scratch ]; then
  echo "tile_shape: node-local storage unavailable." >&2
  echo "  /local/scratch is not a directory inside the container." >&2
  echo "  profiles/nibi/config.yaml apptainer-args must carry --bind /local" >&2
  exit 1
fi
export SP_LOCAL="/local/scratch/sp-{tile}"
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


# Every member of a group job must agree on a STRING resource or snakemake
# raises "Resource slurm_extra is a string but not all jobs in group require
# the same value" (resources.py::_is_string_resource). So the tmp-disk request
# is one constant on all four members, not a property of tile_ngmix alone.
#
# --tmp is a node SELECTION FLOOR, not a reservation: Slurm matches it against
# the node's configured TmpDisk (nibi: 1.67-12 TB, TmpFS=/local) and does not
# decrement it per job. It buys exactly one thing, and it is the thing worth
# having -- the fused job can no longer land on a node with no usable local
# disk, which is the one configuration in which the whole design fails. 16 GB
# against a measured 5.6 GB store leaves room for an unusually rich tile, and
# filters no nibi node today, so it costs no queue time. Resources are not a
# rerun trigger, so adding it invalidates nothing.
TILE_SLURM_EXTRA = "--tmp=16000"

# One label for the whole shape chain. Composition, verified empirically
# against snakemake 9.23.1 (resources.py::GroupResources.basic_layered, and a
# real sbatch): within a toposort level siblings SUM their mem_mb and cores and
# MAX their runtime; across levels the sums are MAXed and the runtimes SUMMED.
# So the eight chunks are one wave -- 8 x 1 core, 8 x mem, ONE chunk's runtime
# -- and the group's declared runtime is the sum along the chain. That sum is
# what gates partitions on nibi, so each member's runtime below is measured
# p99 plus margin, NOT the old defensive ceiling. See the tile_ngmix docstring.
TILE_GROUP = "tile_shape"

# Cleanup, on the LAST member only. $SLURM_TMPDIR is Slurm's to reclaim, but
# /local/scratch on nibi demonstrably carries stale directories from
# long-finished jobs, so the epilog is not on its own reliable. tile_make_cat is
# the last reader of the store, so its EXIT trap is the earliest moment the
# 5.6 GB can go. It must NOT be set by the earlier members: their shells exit
# while the store is still needed. Residual failure mode accepted: SIGKILL
# (OOM-kill, node death) skips the trap, and if the epilog also misses it the
# store leaks until the node drains -- 5.6 GB against a 3.3 TB disk.
TILE_CLEAN = r"""
trap 'rm -rf "$SP_LOCAL"' EXIT
"""

# tile_vignets ONLY. unit_pre clears each stage's own run dir on the shared root
# ("the job clears its run dir at start" — ShapePipe's FileHandler raises on an
# existing one); the node-local root needs the same treatment, and only the rule
# that WRITES the store may do it. Putting this in TILE_LOCAL would have
# tile_ngmix delete the store it is about to read.
TILE_VIGNET_FRESH = r"""
rm -rf "$SP_VIGNET_OUT/run_sp_tile_PiViVi"
"""

# THE FUSE'S ONE SHARP EDGE, made loud rather than mysterious.
#
# The vignette store is not a declared output any more, so snakemake cannot see
# it. tile_vignets' manifest CAN be up to date while the store does not exist on
# this node — concretely, after a fused group job dies part-way: the manifest is
# a success and survives, so a LATER invocation re-plans the group with only the
# surviving members and the store is simply not there. (An in-flight `retries:`
# resubmission is safe: it reruns the same member set, tile_vignets included.)
#
# Recovery is one line, and the message says it: delete the tile's
# tile_vignets.json and resume — that puts tile_vignets back in the group and
# the store comes back with it. The alternative fixes are worse: temp()-ing the
# vignets manifest would break clean_exposure, which keys reclamation
# eligibility on exactly that file.
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




def tile_exp(wc):
    return tile_exposures(wc.tile)

# --- the tile->exposure edge, and why it is cut for finished tiles ----------
#
# THE cascade fix. clean_exposure deletes the exposure's manifests on purpose:
# that is what makes a tile appended later rebuild the chain instead of running
# against an empty store. But those manifests are the tile side's inputs, and an
# exposure is read by ~7-10 tiles. So the moment ONE tile's chain rebuilt an
# exposure, every other tile reading it saw "input files updated by another job"
# and reran — and that rerun rebuilt ITS exposures, which reran THEIR other
# consumers, propagating across the whole exposure-overlap connected component.
# On fixture t4, asking for one damaged tile scheduled all four tiles' chains.
#
# Two mechanisms, and only the second one actually cuts it:
#
#  1. ancient() on every exposure-manifest edge. Correct on its own terms — a
#     tile has no business rerunning because an exposure manifest is NEWER — and
#     it is what keeps a pure-mtime disturbance (a re-touched manifest, a
#     restored backup) from waking finished tiles. But ancient() governs
#     TIMESTAMPS only. Snakemake propagates "my input is produced by a job that
#     will run in this DAG" separately, and ancient does not suppress it
#     (measured: t4 counts were identical with ancient alone).
#
#  2. Cutting the RECLAIMED edges of a FINISHED tile — the mechanism that works.
#     A tile whose final_cat is on disk needs nothing further from its
#     exposures: it has already extracted everything it will ever read. So for
#     such a tile the input list drops the manifests that are GONE, and the
#     propagation has nowhere to go. An UNfinished tile keeps its full edge set
#     and therefore still drags in — and rebuilds — every exposure it needs,
#     which is the accepted price of a late append, unchanged.
#
# Only the missing ones are dropped, never a manifest that still exists: a
# campaign that has cleaned nothing then declares exactly the edges it always
# did, and the cut cannot perturb it.
#
# The marker is final_cat, not the tile's own vignets manifest: on a tile whose
# catalogue was lost, the vignets manifest still exists while the vignette store
# (temp()) does not, so keying on vignets would cut the edge on exactly the tile
# that has to rerun, and run it against a deleted exposure store.
#
# THE CUT REQUIRES THE `input` RERUN-TRIGGER TO BE OFF (profiles/nibi sets the
# trigger list). Dropping an input is itself a change in the set of input files,
# which that trigger reads as a reason to rerun — reinstating the very cascade,
# now as "Set of input files has changed", and running finished tiles against a
# store that is gone. Measured on fixture t4, one damaged tile of four: 82 jobs
# with neither fix, 70 with the cut but the trigger on, 28 with both (= exactly
# the damaged tile's own chain, its two exposures, and the clean jobs).
#
# The cost, stated plainly: `--forcerun` on a tile whose final_cat exists will
# NOT rebuild its reclaimed exposures, because those edges are not in the DAG.
# Delete that tile's final_cat first and the whole chain comes back.
#
# What none of this weakens: clean_exposure's own inputs are neither ancient nor
# cut, so a tile that really did rebuild its vignets still reschedules the cleans
# of the exposures it read, and a grown consumer set still travels through
# params.consumers.
def tile_finished(tile):
    return Path(final_cat(tile)).exists()


def exp_manifests(wc, stage):
    paths = [exp_manifest(e, stage) for e in tile_exp(wc)]
    if tile_finished(wc.tile):
        paths = [p for p in paths if Path(p).exists()]
    return [ancient(p) for p in paths]

def tile_exp_split(wc):  return exp_manifests(wc, "exp_split")
def tile_exp_mask(wc):   return exp_manifests(wc, "exp_mask")
def tile_exp_psf(wc):    return exp_manifests(wc, "exp_psf")
def tile_exp_all(wc):    return tile_exp_split(wc) + tile_exp_mask(wc) + tile_exp_psf(wc)


# Build the per-tile symlink forest. Declaring the exposure manifests as input
# makes this wait on its exposures; the forest itself is only the $SP_EXP view.
# Its output stays a directory() (it has no ShapePipe run dir and no manifest —
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
        # ONCE, so an {output} placeholder inside params.cmd would survive
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
        # config_tile_Mh_exp.ini reads run_sp_tile_Fe output, which the PREPARE
        # phase produced. Declaring the Fe manifest gives the COMPUTE DAG a
        # regeneration path for it instead of a silent dependency on a
        # previous invocation (prepare.smk is included in every parse, so the
        # rule exists here too). Normally a satisfied no-op.
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

# SExtractor object detection on the tile.
rule tile_detect:
    input:
        uz = f"{TILE_DIR}/manifests/tile_uncompress.json",
        mh = rules.tile_merge_headers.output.manifest,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_detect.json"
    log:
        f"{TILE_DIR}/logs/tile_detect.json"
    params:
        pre = lambda wc: unit_pre("tile_detect", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        runtime = 180
    shell:
        sp_shell("tile_detect", "config_tile_Sx.ini")

# PSFEx interpolation to galaxies + vignet postage stamps: the last stage that
# reads exposure products, and the bulk intra-tile intermediate.
#
# The vignette store is declared as a SECOND, temp(directory()) output alongside
# the manifest. This is the ONE scoped exception to "no directory() outputs"
# (D5): the store is ~tens of GB per tile and must be reclaimed when its last
# intra-tile reader finishes, but it must not become DAG currency — so the
# manifest stays the edge, and the directory rides along purely so native temp()
# fires at the right moment. Its readers (tile_ngmix, tile_make_cat) declare
# BOTH. `--notemp` keeps it for debugging.
rule tile_vignets:
    group: TILE_GROUP
    input:
        sx     = rules.tile_detect.output.manifest,
        forest = rules.tile_exp_forest.output.forest,
        split  = tile_exp_split,
        psf    = tile_exp_psf,
        # config_tile_PiViVi.ini reads run_sp_tile_Fe output — same reason as
        # tile_merge_headers above.
        fe     = f"{TILE_DIR}/manifests/tile_find_exposures.json",
    output:
        # THE MANIFEST IS THE ONLY DECLARED OUTPUT NOW. The vignette store used
        # to ride along as a second temp(directory()) output purely so native
        # temp() would reclaim it; it lives on node-local storage under
        # $SLURM_TMPDIR, whose path is not knowable at DAG time, so it cannot
        # be declared at all — and does not need to be. It was never DAG
        # currency (the manifest always was the edge), its readers are now all
        # inside this group job, and Slurm reclaims it with the job. The one
        # affordance lost: `--notemp` can no longer keep it for debugging.
        manifest = f"{TILE_DIR}/manifests/tile_vignets.json",
    log:
        f"{TILE_DIR}/logs/tile_vignets.json"
    params:
        pre = lambda wc: unit_pre("tile_vignets", wc.tile,
                                  forest=forest_dir(wc.tile),
                                  pre_run=[tile_local(wc.tile), TILE_VIGNET_FRESH]),
        script_hash = SCRIPT_HASH
    # 8, not 16, for the same reason tile_ngmix is 1: `-b {threads}` is SMP
    # batch size over input FILE SETS, and a tile is one set -- this run's own
    # log says "Batch size: 16 / Total number of processes: 1". 16 was the
    # widest member and therefore set the whole GROUP's cpus_per_task; at 8 the
    # group asks exactly what the eight-chunk ngmix wave needs. Billing is
    # unchanged either way (nibi is MAX_TRES and the 112 GB memory term is
    # 28 core-equivalents, well above both), but the group now packs onto a
    # node in 8 cores instead of 16.
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 32000 * attempt,
        # Measured median 3m29s, max 5:25 (this branch: 2m38s on 198.305).
        # 20 min is p99 plus a wide margin, and it is a term in the GROUP's
        # runtime sum, so the old defensive 240 is not free any more.
        runtime = 20,
        slurm_extra = TILE_SLURM_EXTRA
    shell:
        # The completeness check is pointed at the NODE-LOCAL run root.
        # completeness.py is unmodified: --run-dir is its own documented override
        # for "$SP_RUN", and --unit keeps the manifest's unit field the tile ID
        # instead of the basename of the local root (it would otherwise default
        # to "sp-tile"). See sp_shell for the rest.
        sp_shell("tile_vignets", "config_tile_PiViVi.ini",
                 check_args=' --run-dir "$SP_LOCAL" --unit {wildcards.tile}')

# ngmix shape measurement — N chunks per tile (D4). Each chunk computes its own
# CLOSED object-ID range at EXECUTION time from this tile's own sexcat: a params
# function cannot, because params evaluate before the sexcat exists. Closed, not
# open-ended: `ID_OBJ_MAX = -1` on the last chunk was the 13-hour straggler's
# root cause (ngmix treats id_obj_max <= 0 as unbounded).
#
# Chunks write nothing shared: each has its own run_sp_tile_ngmix_Ng<k>u, and
# merge_sep_cats — DAG-serialised after all chunks — is the gather.
rule tile_ngmix:
    group: TILE_GROUP
    input:
        # The manifest, and only the manifest. The vignette store is no longer a
        # declared input because it is no longer a declared output: it is
        # node-local, produced by tile_vignets earlier in THIS SAME group job
        # (see TILE_LOCAL). The manifest was always the real edge.
        vignets  = rules.tile_vignets.output.manifest,
        sx       = rules.tile_detect.output.manifest,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_ngmix_{{chunk}}.json",
        # STAYS ON SCRATCH, unlike the vignette store. It is DAG currency:
        # tile_merge_cats reads it, and merge_sep_cats derives chunks 2..N from
        # chunk 1's path. It is also small (~300 KB/chunk), so it costs the
        # scratch high-water nothing worth chasing. temp() still reclaims it
        # once merge_cats has run.
        chunkdir = temp(directory(f"{TILE_DIR}/output/run_sp_tile_ngmix_Ng{{chunk}}u")),
    log:
        f"{TILE_DIR}/logs/tile_ngmix_{{chunk}}.json"
    params:
        pre = lambda wc: unit_pre("tile_ngmix", wc.tile,
            env={"SP_NGMIX_CHUNK": wc.chunk, "NGMIX_N_CHUNKS": NGMIX_CHUNKS},
            # Two steps, not `eval "$(...)"`: a command substitution inside eval
            # discards the script's exit status, so a missing sexcat would fall
            # through to shapepipe_run with an unset range and fail as something
            # else. Capture, check, then eval — the range script fails as itself.
            pre_run=[f'ngmix_range_out=$(python {SCRIPTS}/ngmix_range.py --run-dir '
                     f'"$SP_RUN" --chunk {wc.chunk} --n-chunks {NGMIX_CHUNKS}) '
                     f'|| exit 1',
                     'eval "$ngmix_range_out"',
                     tile_local(wc.tile), TILE_VIGNET_REQUIRED]),
        script_hash = SCRIPT_HASH,
        # The one rule whose fingerprint has to cover ngmix_range.py.
        # `pre_run` already puts the INVOCATION in params, but the invocation
        # is invariant to the script's body — and what that body decides is a
        # PARTITION. Resume a tile across an edit to the split without this and
        # the already-done chunks keep the old ranges while the reruns take the
        # new ones: some objects measured twice, others by nobody, the tile
        # green, no error anywhere.
        #
        # ADDING THIS PARAM IS ITSELF A PARAMS CHANGE, AND ON A RESUME IT DOES
        # THE VERY THING IT PREVENTS, BY ANOTHER ROUTE. Snakemake compares the
        # SET of recorded param values, and smk-g4's metadata for a finished
        # chunk (record_format_version 6) holds two of them -- the SCRIPT_HASH
        # digest and the pre string -- so a third entry is new and the job
        # replans with "Params have changed since last execution". On snakemake
        # 9.23.1 under profiles/nibi's trigger set, on a fixture mirroring this
        # group with everything up to date, that reschedules the eight chunks,
        # tile_merge_cats and tile_make_cat -- and NOT tile_vignets or
        # tile_detect, whose manifests exist, so missing_output never queues
        # them and nothing propagates down to them either.
        #
        # A group replanned without tile_vignets is precisely the state
        # TILE_VIGNET_REQUIRED exists to catch. Every chunk trips the guard and
        # the group fails; from there GroupJob.postprocess(error=True) fans out
        # over every member and removes its EXISTING outputs -- tile_make_cat's
        # final_cat on the persistent root among them. So the cost of getting
        # this wrong is a deleted science product, not a wasted hour.
        #
        # ALL OF THIS IS OBSERVED, end to end, on a real campaign root -- SLURM
        # job 20818649 against smk-g5, 2026-08-30. It had been derived twice from
        # snakemake's cleanup path and never seen, because "Group jobs: inactive
        # (local execution)" puts it out of reach of any login-node fixture. The
        # run took 32 seconds:
        #   * the dry run scheduled tile_ngmix x8 + tile_merge_cats +
        #     tile_make_cat, reason "params have changed since last execution",
        #     with tile_vignets and tile_detect ABSENT -- the group's own SLURM
        #     label came back as tile_shape_tile_make_cat_tile_merge_cats_tile_ngmix,
        #     with no tile_vignets in it;
        #   * every chunk tripped the guard below;
        #   * the group FAILED, and final_cat-186.307.fits was GONE from
        #     /project afterwards. The tile's directory on the persistent root
        #     was empty.
        # (Restored from a copy taken first; the design and the transcript are in
        # sp-products/smk-g5/EXPERIMENT_postprocess_deletion.md.)
        #
        # WORTH SITTING WITH: the guard below is what makes this loud rather than
        # silent, and the loud failure is precisely what triggers the deletion. A
        # guard that turns a wrong catalogue into a failed job is still the right
        # trade -- but on a fused group, "fails loudly" is not a cheap outcome,
        # and that is an argument for never reaching this state rather than for
        # relaxing the guard.
        #
        # WHY THE TILE_LOCAL WARNING ABOVE DOES NOT ALREADY COVER IT, which is
        # the whole reason this needs its own note: tile_local() sits in the
        # params.pre of tile_vignets, tile_ngmix AND tile_make_cat, so an edit
        # there drags tile_vignets back in and the store is rebuilt.
        # range_hash is on tile_ngmix ALONE. Same family, opposite shape -- it
        # is the ONE-rule fingerprint change that is dangerous, not the
        # three-rule kind that warning describes, and reading it as covering
        # both is the mistake to avoid. job_head.sh's runbook sweep for the
        # guard state does not help either: it skips any tile with a final_cat,
        # which is exactly the set this breaks.
        #
        # INVERTED BY clean_tiles, so read the default carefully. With
        # reclamation ON the tombstoned tile has lost tile_detect.json, and the
        # structural "Input files updated by another job" propagation puts
        # tile_vignets back in the group (clean_tile below measures this) --
        # the store is rebuilt and nothing trips. The SHIPPED DEFAULT
        # clean_tiles: false is therefore the dangerous configuration, which is
        # the opposite of how reclamation reads everywhere else in this file.
        #
        # So: land this hash, and every later edit to ngmix_range.py, at a
        # campaign boundary on a fresh root. Same rule as tile_local()'s, same
        # escape hatch and same caveat if a resume is unavoidable --
        # `--rerun-triggers mtime code software-env`.
        range_hash = NGMIX_RANGE_HASH
    # ONE core, not four. `-b {threads}` is shapepipe_run's SMP BATCH SIZE
    # (pipeline/args.py) -- joblib Parallel(n_jobs=batch_size) over
    # filehd.process_list, i.e. parallelism ACROSS INPUT FILE SETS. An ngmix
    # chunk is one catalogue, so process_list has exactly one entry: every real
    # log says "Batch size: 4 / Total number of processes: 1". There is no
    # internal parallelism either (no multiprocessing/Pool/joblib/threading in
    # ngmix_package/ngmix.py), and OMP/BLAS are pinned to 1 by both the prologue
    # and the profile. The reserved cores 2-4 never had anything to run.
    #
    # Worth being precise about what this saves. nibi bills
    # TRESBillingWeights=CPU=1000,Mem=250G under PriorityFlags=MAX_TRES, i.e.
    # max(cores, mem_GB/4). At 14 GB the memory term alone is 3.5 core-
    # equivalents, so 4 -> 1 core moves the reservation from 4.0 to 3.5, a 12%
    # saving -- NOT 75%. The real core-hour win is the elapsed-time collapse
    # from killing the NFS random reads, not this.
    threads: 1
    retries: 2
    benchmark:
        f"{TILE_DIR}/manifests/tile_ngmix_{{chunk}}.benchmark.tsv"
    resources:
        # 5000, down from 14000, and this is the campaign's largest single cost
        # saving — but read what the number means before moving it again.
        #
        # sacct's MaxRSS here is NOT process memory. nibi runs
        # JobAcctGatherType=jobacct_gather/cgroup, so it reports the cgroup's
        # memory.current, which under cgroup v2 CHARGES PAGE CACHE to the job.
        # Cache is reclaimable — the kernel evicts it before it kills anything —
        # so the cgroup high-water is an upper bound on what the job NEEDS, not
        # a floor. What a tight reservation can still do is squeeze the cache
        # that keeps the node-local store resident, which is part of why the
        # fused tile is fast.
        #
        # MEASURED ACROSS 31 TILES (campaign smk-g4, 2026-08-30): cgroup
        # high-water 27.40 GiB max, 22.14 GiB mean, against the 109.4 GiB the
        # group was reserving. The eight chunks are SIBLINGS in the group's
        # toposort so snakemake SUMS their mem_mb; at 5000 the group asks
        # 40000 MiB = 39.1 GiB, a 1.43x margin over the worst tile observed.
        #
        # The ANONYMOUS half of that, which is the part that can actually OOM,
        # is the number to argue about, and this repo carries two estimates that
        # disagree by 1.8x. Snakemake's psutil benchmark says ~1.25 GiB per
        # chunk (~10 GiB for eight); config_tile_Ng_template.ini's own
        # SAVE_BATCH note says ~2.3 GB per worker from an A/B test (job
        # 17607877), i.e. ~18.4 GiB. Prefer the larger: snakemake samples RSS on
        # a 30-second grid (BENCHMARK_INTERVAL), and a SAVE_BATCH = 250 flush
        # cycle is exactly the sawtooth such a grid misses. So read the margin
        # as ~2.1x over anonymous memory with ~20 GiB left for cache against an
        # 8.1 GiB store — comfortable, but not the 4x a psutil-only reading
        # would suggest. An OOM is also self-healing: mem_mb scales with
        # attempt, so a retry asks 80000.
        #
        # WHY THIS IS THE BIG ONE. nibi bills max(cores, mem_GB/4), so at
        # 112 GB the fused group billed 28 core-equivalents for 8 real cores —
        # 3.4x, measured live at 905 billed against 266 allocated. At 40 GB it
        # bills 10. That is not a saving so much as a schedule: the account's
        # fairshare target is ~250 CE, which buys 9 tiles in flight at 28 and 25
        # at 10, and DR6's wall clock is (tiles / tiles-in-flight) x elapsed.
        #
        # 4000 is the floor and is NOT recommended yet: it would take the group
        # to 32000 MiB, where tile_vignets' own 32000 becomes the binding term
        # and the group finally bills its 8 real cores — but that is a 1.14x
        # margin over the worst tile measured, and the first thing to give would
        # be the page cache holding the store. Take it only with a measurement
        # of cache behaviour under pressure, not on the arithmetic alone.
        mem_mb = lambda wc, attempt: 5000 * attempt,
        # 120 on the FIRST attempt, and ATTEMPT-SCALED after it. MEASURED
        # two ways, and the second is why the margin is thinner than it looks:
        #   * alone (job 20795277, tile 198.305 chunk 1): ~76 min for 3547
        #     objects = 1.29 s/object, against a 7h34m median on NFS -- the
        #     same ~50 min of TotalCPU either way, so the collapse is pure I/O.
        #   * EIGHT-WIDE on one node (job 20799387, tile 186.307, 4412
        #     objects/chunk): 1.54 s/object steady-state, i.e. concurrency
        #     costs ~19%, and the chunk lands at ~113 min. The campaign's
        #     largest tile (198.306, 4678 objects/chunk) projects to ~120 --
        #     exactly this number, with nothing left over.
        #
        # Inside a group SLURM enforces only the GROUP's wall (165 min), never
        # a member's, so 120 is a budgeting term rather than a kill line and
        # the worst tile still lands ~127 min inside 165. What is NOT safe is a
        # flat retry: a group that TIMEOUTs re-queues against the identical
        # wall and fails identically, burning three 165-minute allocations to
        # learn nothing. Scaling with `attempt` keeps the happy path in
        # cpubase_bycore_b1 (20+120+10+15 = 165 <= 180) and gives a retry real
        # headroom (285 min, which is b2) instead of a rerun of the same
        # failure. Attempt 1 is unchanged, so this does not perturb the
        # benchmark -- it only makes the failure branch mean something.
        runtime = lambda wc, attempt: 120 * attempt,
        slurm_extra = TILE_SLURM_EXTRA
    shell:
        sp_shell("tile_ngmix", "config_tile_Ng_template.ini")


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
        sp_shell("tile_merge_cats", "config_merge_sep_cats.ini")

# The run's science product. make_cat also reads the vignette store's
# psfex_interp output, so it — not ngmix — is the store's last reader.
#
# No protected(): the full default rerun-triggers govern, and protected() only
# ever forced people through a `--forcerun` detour.
rule tile_make_cat:
    group: TILE_GROUP
    input:
        # No store input: it is node-local, written by tile_vignets in this same
        # group job. make_cat reads its psfex_interp output through
        # $NGMIX_VIGNET_DIR (config_tile_Mc.ini).
        ms    = rules.tile_merge_cats.output.manifest,
    output:
        manifest  = f"{TILE_DIR}/manifests/tile_make_cat.json",
        final_cat = f"{PROD_TILE_DIR}/final_cat-{{tile}}.fits",
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
                      '  cp -f "$(ls -1 "$SP_RUN"/output/run_sp_Mc/make_cat_runner'
                      '/output/final_cat*.fits | head -1)" {output.final_cat}\n'
                      "fi\n")


# --- tile reclamation (D5) --------------------------------------------------
# The tile-side counterpart of clean_exposure, and the differences are the whole
# design. Read that rule's commentary in exposure.smk first; this one only says
# where the tile case departs from it.
#
# WHY IT EXISTS. clean_exposure reclaims exposure stores; nothing reclaimed
# tiles, so a FINISHED tile held 1.19 GiB across 137 inodes on scratch forever
# (measured on 186.307; the per-directory breakdown is in clean_tile.py). At
# DR6's 23,114 tiles that is 26.9 TiB against a 1 TiB quota and 3.17M inodes
# against a 1M one. BOTH BOUNDS BIND, and the byte one binds first: without this
# rule no batch may exceed ~859 tiles. The inode arithmetic is why the whole
# directory goes rather than just `output/` — the 62 non-output inodes per tile
# are still 1.43M at DR6, i.e. over quota on their own.
#
# ELIGIBILITY IS TRIVIAL, AND THAT IS THE POINT. An exposure's consumer set
# closes over the CAMPAIGN (~7-10 tiles read each one), which is why
# clean_targets() has to test every consumer and why clean_exposure carries the
# set in params so a late append makes the tombstone stale. A tile's store has
# NO consumer outside that tile: tiles read exposures, and nothing reads another
# tile's store. So there is no consumer set, no eligibility test, no staleness
# to detect — only an ordering edge on the tile's own final_cat.
#
# THE INPUT IS final_cat, ON THE PERSISTENT ROOT, AND NOT ancient(). It is
# already the campaign's designated tile-finished marker (see final_cat()'s
# docstring and tile_finished() above), and an ordinary input is what orders the
# clean after the tile — which is what makes reclamation ROLLING: within one
# invocation a tile is reclaimed as soon as its own chain lands, so the scratch
# high-water tracks the tiles in flight instead of the tiles in the batch. That
# is the whole reason the rule is worth having; ancient() would keep the
# dependency but drop the ordering, and the batch would peak at its full size.
#
# SCOPE: TILES_READY only, never every tile in the index. clean_targets() has to
# reason about eligibility; this one has to reason about the converse, and for a
# different reason. Requesting an OUT-OF-SCOPE tile's tombstone would put that
# tile's final_cat — hence its tile_make_cat job, hence the whole fused
# tile_shape group — into this DAG, where the `params` rerun-trigger can find it
# out of date (see the TILE_LOCAL warning above) and reschedule a finished tile
# against exposure stores clean_exposure has already reclaimed. In scope, that
# costs nothing: `rule all` requests the same final_cat anyway, so the clean adds
# an edge and no jobs. Deferral, never loss: the tile list accumulates, so every
# later invocation has the tile in scope and picks it up — which is also how
# flipping `clean_tiles` on later reclaims retroactively (the tiles already
# cleaned are exactly the ones with a tombstone).
#
# FOUR SURVIVORS, not one, and clean_tile.py argues each. Two of them are the
# corrections this rule's first design needed: `manifests/tile_find_exposures.json`
# (rule prepare_all_tiles declares it for EVERY tile in the list, and the list
# accumulates — delete it and the prepare group reruns per cleaned tile per
# invocation, writing uncompress's 382 MB back into the store) and the Fe
# exposure list under `output/` (build_index.exp_list_path re-reads it at EVERY
# compute parse; delete it and the default SP_MISSING_THRESHOLD of 0.0 makes the
# first `sp run` after the first reclamation a fatal parse). Ten inodes per tile
# survive in total, ~231k at DR6.
#
# THE SURVIVING tile_vignets.json DOES NOT TRIP TILE_VIGNET_REQUIRED, AND THE
# REASON IS WHAT GETS DELETED, NOT WHAT GETS KEPT. Read that guard above first.
# A cleaned tile is, on its face, exactly the state it exists to catch: a valid
# tile_vignets manifest over a node-local store that is not there (here because
# it never outlived its job, rather than because a fused group half-ran).
# smk-g4's job_head.sh carries a runbook sweep for that state, and it skips any
# tile with a final_cat — i.e. every tile this rule ever touches — so a cleaned
# tile is invisible to it either way.
#
# Force a cleaned tile back anyway (delete its final_cat) and tile_vignets IS
# rescheduled, so the store is rebuilt and no chunk trips the guard. Measured on
# the fixture: 19 jobs, the whole chain from tile_get_images, tile_vignets inside
# the tile_shape group. Snakemake names the mechanism itself --
#   reason: Input files updated by another job: .../exp_forest,
#           .../manifests/tile_find_exposures.json, .../manifests/tile_detect.json
# -- and it is NOT the `mtime` trigger. It is the structural propagation that
# makes every dependent of a rerunning job rerun, the same edge the cascade
# commentary above says ancient() cannot suppress. Worth the distinction: that
# propagation is not in profiles/nibi's rerun-triggers list and so cannot be
# switched off there, where an mtime argument could be.
#
# So the guard stays silent only because tile_detect.json and exp_forest/ GO.
# Counterfactual on the same fixture, upstream manifests restored and mtimes
# controlled so nothing reruns on timestamp -- i.e. the "delete only output/"
# design -- schedules 8 x tile_ngmix + tile_merge_cats + tile_make_cat with
# tile_vignets ABSENT, and every chunk trips the guard.
# ADDING A MANIFEST TO THE SURVIVOR LIST IS THEREFORE NOT FREE: tile_vignets is
# safe because nothing upstream of it survives, and tile_detect.json in
# particular must never join the list.
#
# A LOCALRULE (declared in the Snakefile), same as clean_exposure: it is an
# rmtree, not science, and one sbatch per tile is ~23k scheduler submissions at
# DR6 scale for work shorter than the scheduling latency. Unlike exp_star_cat it
# is a DAG LEAF, so being local can never make it both a dependency and a
# dependent of a group and there is no cycle to avoid — no `group:` label here,
# and none possible.
#
# WHAT THIS SHARPENS ELSEWHERE: the TILE_LOCAL warning above says an edit to
# tile_local() mid-campaign reruns finished tiles unsatisfiably because their
# EXPOSURE stores are reclaimed. With this rule on, the tile's own store is gone
# too, so the rerun has even less to stand on. The recommendation is unchanged —
# land those edits between campaigns — but the margin for being wrong is smaller.
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
