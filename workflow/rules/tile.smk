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
rows to it before the multi-epoch post-processing (MATCH_CATALOGUE, set
through SP_MATCH_CATALOGUE). Each row pairs with its mutual nearest catalogue
object within 1 px and takes its NUMBER, from which make_cat builds
``TILE_UNIQUE_ID``, the key shared with the photometry and photo-z
catalogues; unpaired rows leave the catalogue. The catalogue is SExtractor run
by MegaPipe on the same DR6 image with the same configuration, so the join is
exact; on other pixels (a DR5 image) it fails loudly.

``blend_handling: uberseg``, the default, gives ngmix UberSeg's hard mask
instead of its own noise-fill default. Both read one coadd segmentation stamp
per object from the sexcat's SEG_VIGNET column: uberseg masks by it, and
noisefill takes as neighbours only the -1e30 VIGNET markers on another
object's footprint in it. tile_detect cuts it from SExtractor's SEGMENTATION
check image on the grid of each object's VIGNET, before the join, which
relabels every stamp to the catalogue's NUMBER (-1 on the footprints of rows
that leave, 0 on sky). The switch is one prologue variable and nothing else
(blend_env, below).

There is no `tile_mask` rule, and there will not be one (PR #847). ShapePipe
generates no masks: tiles have no instrument flag image of their own, so
tile_detect runs SExtractor with FLAG_IMAGE = False against
default_noimaflags.param (config_tile_Sx.ini — what used to be the "sx_nomask"
variant, now the only one). Sky-fixed masks reach the tile as CATALOGUE columns
instead: ``tile_make_cat``'s make_cat queries the configured healsparse maps at
every object's (RA, Dec) and writes one ``MASK_<band>`` column per band, which
is what downstream selections cut on.
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
if [ ! -d "$NGMIX_VIGNET_DIR/vignetmaker_runner/output" ]; then
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

# --- blend handling -------------------------------------------------------

# The blend-handling exports into a stage's prologue. tile_detect writes
# SEG_VIGNET (config_tile_Sx.ini, SEG_VIGNET = ${SP_SEG_VIGNET:-False}) under
# both handlings: noisefill keeps only the -1e30 markers on another object's
# footprint in it (split_tile_markers), uberseg masks by it. uberseg (the
# default) also sets BLEND_HANDLING (config_tile_Ng_template.ini,
# ${SP_BLEND_HANDLING:-noisefill}), which reruns only tile_ngmix on a flip.
BLEND_ENV = {
    "tile_detect": {"SP_SEG_VIGNET": "True"},
    **({"tile_ngmix": {"SP_BLEND_HANDLING": "uberseg"}}
       if BLEND_HANDLING == "uberseg" else {}),
}


def blend_env(stage):
    return BLEND_ENV.get(stage, {})


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
    """tile_detect's prologue exports: the catalogue to join, if any,
    SP_SEG_VIGNET, and for image simulations the epoch-membership CCD bounds.

    SP_MATCH_CATALOGUE is empty under tile_detection: sextractor (image
    simulations), so that no value left exported in the submitting shell
    reaches the job. SP_CCD_SIZE is exported only for image simulations (their
    trimmed CCDs); data runs take config_tile_Sx.ini's default.
    """
    env = {"SP_MATCH_CATALOGUE": "", **blend_env("tile_detect")}
    if TILE_DETECTION == "unions_catalogue":
        gic = f"{tile_dir(tile)}/output/run_sp_tile_Gic/get_images_runner/output"
        env["SP_MATCH_CATALOGUE"] = f"{gic}/CFIS_cat{unit_num(tile)}.cat"
    elif INPUT_TYPE == "image_sims":
        env["SP_CCD_SIZE"] = config["ccd_size"]
    return env


# tile_detect's extra memory for SEG_VIGNET: add_seg_vignet reads the 400 MB
# SEGMENTATION check image and adds an int32 SEG_VIGNET column the size of
# VIGNET, which the join then holds twice. MEASURED on DR6 tile 202.301
# (36,022 detections; sextractor_runner as tile_detect runs it, without the
# post-processing): peak RSS 1.53 GiB without SEG_VIGNET, 2.93 GiB with it.
# Scaled to the 1.84 GiB worst tile above, it needs ~3.5 GiB; 7000 MB is 2x
# that.
DETECT_SEG_MEM_MB = 3000


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
    # One core: the container's SExtractor is built without threads
    # (NTHREADS 4 warns and runs in the same 36 s), default_tile.sex sets
    # NTHREADS 1, the join and the post-processing are serial Python, and
    # `-b {threads}` only sets the SMP batch over input sets, of which a tile
    # is one.
    #
    # MEASURED on eight DR6 tiles chosen to span conditions (high latitude,
    # b = 18 deg, the A2199 cluster, Alioth's halo, two survey-edge tiles at
    # 96% zero weight; 1-10 exposures), the step as shapepipe_run runs it on
    # candide: wall time 23-160 s (SExtractor + join <= 93 s, post-processing
    # <= 97 s), peak RSS 1.84 GiB, set by the join holding the ~430 MB
    # catalogue twice. Image simulations run the same step on tiles of the same
    # size (10000 x 10000 px); SExtractor alone on one takes 34 s and 0.74 GiB.
    # 4000 MB is 2.1x the worst tile, and an OOM retries at 8000. The runtime
    # is ~7x the slowest tile, for /scratch I/O on nibi (a 400 MB image and
    # weight in, the catalogue and check images out). The SEG_VIGNET cut, the
    # join and the post-processing each rewrite the whole catalogue, so
    # config_tile_Sx.ini's WORK_DIR keeps them on node-local disk: MEASURED on
    # DR6 191.307 under uberseg (31,736 rows, 8 epochs), the step took 44 min
    # in the loaded smk-g13-us allocation with the rewrites on /scratch, 112 s
    # on idle /scratch, and 36 s through WORK_DIR. nibi bills
    # max(cores, mem_GB/4), so the job bills 1 core-equivalent; plus
    # DETECT_SEG_MEM_MB above.
    threads: 1
    retries: 1
    resources:
        mem_mb = lambda wc, attempt: (4000 + DETECT_SEG_MEM_MB) * attempt,
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
        # tile_vignets now PRODUCES the chunk partition (TILE_NGMIX_RANGES), so
        # it carries the splitter's fingerprint too. Same hash, same constant —
        # the Snakefile's NGMIX_RANGE_HASH argues what it guards, and
        # tile_ngmix's copy below carries the mid-campaign-edit warning.
        range_hash = NGMIX_RANGE_HASH,
        # tile_ngmix exports them, so they are a params change there; carried
        # here too so that changing either reruns the whole group from
        # tile_vignets, never the chunks alone (see the mid-campaign deletion
        # note on tile_ngmix's range_hash).
        defect_fill = DEFECT_FILL,
        defect_weighting = DEFECT_WEIGHTING
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

# Each ngmix chunk's extra memory for SEG_VIGNET: Tile_cat holding the tile's
# SEG_VIGNET column, a view into the loaded table like VIGNET. Measured on DR6
# 202.301 (36065 rows, 51x51): Tile_cat peak and resident RSS 955 MiB with the
# column, 595 without (+360 MiB, the column's size). 500 covers 186.307's
# ~37.9k rows (~375 MiB) with a third to spare.
NGMIX_SEG_MEM_MB = 500

# ngmix shape measurement — N chunks per tile (D4). Each chunk LOOKS UP its own
# CLOSED catalogue-row range in the file tile_vignets materialised at the top of this
# group job (TILE_NGMIX_RANGES); the ranges are knowable only at EXECUTION time,
# from this tile's own sexcat, which is why a params function cannot supply them.
# Closed, not open-ended: `ID_OBJ_MAX = -1` on the last chunk was the 13-hour
# straggler's root cause (ngmix treats id_obj_max <= 0 as unbounded).
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
            env={"SP_NGMIX_CHUNK": wc.chunk, "NGMIX_N_CHUNKS": NGMIX_CHUNKS,
                 "SP_DEFECT_FILL": DEFECT_FILL,
                 "SP_DEFECT_WEIGHTING": DEFECT_WEIGHTING,
                 **blend_env("tile_ngmix")},
            # Two steps, not `eval "$(...)"`: a command substitution inside eval
            # discards the script's exit status, so a missing sexcat would fall
            # through to shapepipe_run with an unset range and fail as something
            # else. Capture, check, then eval — the range script fails as itself.
            # tile_local FIRST, because it is what exports $SP_LOCAL — and the
            # ranges file the lookup reads lives there, written once by
            # tile_vignets (TILE_NGMIX_RANGES). This chunk only looks its row
            # up; it never recomputes, and the script refuses to.
            #
            # Two steps, not `eval "$(...)"`: a command substitution inside eval
            # discards the script's exit status, so a missing ranges file would
            # fall through to shapepipe_run with an unset range and fail as
            # something else. Capture, check, then eval — the lookup fails as
            # itself.
            pre_run=[tile_local(wc.tile), TILE_VIGNET_REQUIRED,
                     f'ngmix_range_out=$(python {SCRIPTS}/ngmix_range.py '
                     f'--read "{NGMIX_RANGES}" --chunk {wc.chunk}) || exit 1',
                     'eval "$ngmix_range_out"']),
        script_hash = SCRIPT_HASH,
        # `pre_run` already puts the INVOCATION in params, but the invocation is
        # invariant to the script's body, and what that body decides is a
        # PARTITION (the Snakefile's NGMIX_RANGE_HASH argues the corruption).
        # Carried on tile_vignets too, since the split moved there.
        #
        # WHY THE PARTITION IS MATERIALISED ONCE, kept as history because the
        # design only reads as over-careful until you have seen this. The eight
        # chunks used to each run ngmix_range.py in their own shell and trust the
        # others to have landed on the same boundaries. A fused group holds them
        # open for hours (smk-g4: 6,236-7,762 s elapsed per chunk), and each chunk
        # read the script only when its own shell started -- against the LIVE
        # checkout, as every job did before the launch code snapshot (bin/sp) --
        # so an edit landed mid-flight was read by some and not others.
        # Seen once, live: an invocation four seconds after a rewrite returned
        # chunk 5 of 186.307 as 17649..22060 where the other seven had been split
        # 17129..21229 — 520 objects orphaned, 831 measured twice, and
        # merge_sep_cats concatenates whatever it is handed, so the tile would
        # have completed green. (Twelve sequential and thirty-six concurrent runs
        # against a stable checkout gave the identical correct partition; the
        # race was purely the edit.) That is now STRUCTURALLY closed: tile_vignets
        # writes the ranges once and the chunks only look their row up, so there
        # is nothing left for them to disagree about. The hash below no longer
        # buys sibling agreement — it buys REPRODUCIBILITY across a resume.
        #
        # WHAT A ONE-RULE FINGERPRINT DOES MID-CAMPAIGN: IT DELETES SCIENCE
        # PRODUCTS. Observed when range_hash lived on tile_ngmix ALONE, which is
        # no longer the case -- tile_vignets carries it too now that it is the
        # rule that runs the splitter, so an edit drags tile_vignets back into
        # the group, the store and the ranges are rebuilt, and nothing trips.
        # Kept because the mechanism is general and the next single-rule param is
        # one edit away. A new param entry replans the fused group with "Params have changed since last
        # execution", scheduling the eight chunks, tile_merge_cats and
        # tile_make_cat but NOT tile_vignets or tile_detect (their manifests
        # exist, so missing_output never queues them). That is exactly the state
        # TILE_VIGNET_REQUIRED catches: every chunk trips the guard, the group
        # fails, and GroupJob.postprocess(error=True) removes every member's
        # EXISTING outputs -- tile_make_cat's final_cat on the persistent root
        # among them.
        #
        # OBSERVED end to end, not derived: SLURM job 20818649 against smk-g5,
        # 2026-08-30, 32 seconds. The dry run scheduled exactly that member set
        # (the group's SLURM label came back
        # tile_shape_tile_make_cat_tile_merge_cats_tile_ngmix, no tile_vignets);
        # every chunk tripped the guard; the group failed; and
        # final_cat-186.307.fits was GONE from /project afterwards, the tile's
        # directory on the persistent root empty. It had never been seen before
        # because "Group jobs: inactive (local execution)" puts it out of reach
        # of any login-node fixture. (Restored from a copy taken first; design
        # and transcript in sp-products/smk-g5/EXPERIMENT_postprocess_deletion.md.)
        #
        # WORTH SITTING WITH: the guard is what makes this loud rather than
        # silent, and the loud failure is precisely what triggers the deletion.
        # Still the right trade -- but it argues for never reaching this state,
        # not for relaxing the guard.
        #
        # THE RULE THE INCIDENT LEFT BEHIND: a fingerprint that changes on the
        # chunks but NOT on tile_vignets is the dangerous shape. tile_local()
        # was always safe for exactly this reason -- it sits in three rules'
        # params.pre, so an edit there pulls tile_vignets in and the store is
        # rebuilt. range_hash is now the same shape by construction. job_head.sh's
        # runbook sweep would not have helped either: it skips any tile with a
        # final_cat, exactly the set this breaks.
        #
        # INVERTED BY clean_tiles, so read the default carefully. With
        # reclamation ON the tombstoned tile has lost tile_detect.json and the
        # structural "Input files updated by another job" propagation puts
        # tile_vignets back in the group (measured at clean_tile below), so the
        # store is rebuilt and nothing trips. A campaign override setting
        # clean_tiles: false leaves the dangerous configuration; the shipped
        # default is true. Epoch-cut diagnostics no longer need that override.
        #
        # So: land this hash, and every later edit to ngmix_range.py, at a
        # campaign boundary on a fresh root. The same rule and direct command as
        # tile_local()'s -- `--rerun-triggers mtime code software-env`.
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
        # 4000 is the current lower bound and is NOT recommended yet: it would take the group
        # to 32000 MiB, where tile_vignets' own 32000 becomes the binding term
        # and the group finally bills its 8 real cores — but that is a 1.14x
        # margin over the worst tile measured, and the first thing to give would
        # be the page cache holding the store. Take it only with a measurement
        # of cache behaviour under pressure, not on the arithmetic alone.
        #
        # Each chunk also holds the tile's SEG_VIGNET column,
        # int32 and the size of VIGNET: NGMIX_SEG_MEM_MB below.
        mem_mb = lambda wc, attempt: (5000 + NGMIX_SEG_MEM_MB) * attempt,
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
        stats = f"{PROD_TILE_DIR}/tile_make_cat.json",
    log:
        f"{TILE_DIR}/logs/tile_make_cat.json"
    params:
        pre = lambda wc: unit_pre("tile_make_cat", wc.tile,
                                  env={"NGMIX_N_CHUNKS": NGMIX_CHUNKS},
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
                      '  cp -f {output.manifest} {output.stats}.tmp\n'
                      '  mv -f {output.stats}.tmp {output.stats}\n'
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
        cat = lambda wc: final_cat(wc.tile),
        stats = rules.tile_make_cat.output.stats,
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
        " --tombstone {output.tombstone} --final-cat {input.cat}"


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
