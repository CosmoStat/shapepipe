"""Exposure chain — per exposure, keyed by exp base id (dedup is structural).

    exp_get_images -> exp_split -> exp_psf -> exp_persist

Each in the exposure's own sharded work dir, chained by manifests; every config
reads fixed ``$SP_RUN/output/run_sp_exp_*`` INPUT_DIRs, so nothing resolves a
run log. There is no `prepare_exposures` aggregation target: these chains hang
off the compute DAG (`all` <- final_cat <- tile chain <- exposure manifests).

NO MASK RULE, and that is the design (PR #847). ShapePipe generates no masks.
The only mask that reaches pixels is the instrument flag image delivered with
the exposure, which ``exp_split`` splits per CCD alongside image and weight and
SExtractor reads directly. Sky-fixed masks are healsparse maps, queried once per
object: ``mask_query`` (inside exp_psf's config chain) writes ``FLAG_EXT`` onto
each CCD's SExtractor catalogue for setools' star cut, and ``make_cat`` writes
the per-band ``MASK_<band>`` columns on the tile side. Neither needs a rule, a
star catalogue, or a network fetch — hence no ``star_catalogue`` / ``exp_star_cat``
here, and no ``exp_mask``.

``exp_persist`` is the one rule here that writes to the PERSISTENT root: it
packs the PSF products named by `persist_exp:` into one tar per exposure off
/scratch before the purge (or clean_exposure) can take them. It is a separate
rule from exp_psf precisely so that editing that list costs a re-pack and not a
four-hour refit; the full
argument is in workflow/scripts/persist_exp.py.

NO temp() anywhere in this file, ever (D5). Exposures overlap tiles by
construction (~7-10 tiles each), so their consumer set closes over the CAMPAIGN,
not over one invocation — reclamation here is clean_exposure's job (S5), driven
by the accumulating index. A temp() here would delete an exposure the moment
this invocation's readers finished and cascade destructive reruns across spatial
neighbours the next time a tile is appended.

NO GROUPING. The ``exp_short`` group existed to fuse exp_split and exp_mask —
two rules whose medians were 1:28 and 1:54, both well under the 15-minute floor
Alliance policy asks us to bundle away — into one sbatch per exposure. With
exp_mask gone there is nothing to fuse: a group of one rule submits exactly the
job the ungrouped rule submits, and the label would only obscure that. The
composition rules, should a second short rule ever appear here, are in
prepare.smk's docstring. exp_get_images stays separate for the same reason it
always did (a download, retried on its own), and exp_psf is heavy (16 GB, 4 h)
and never fuses with a short rule.

NUMBER_LIST ($SP_UNIT_NUM, see unit_num in the Snakefile) is set only for
exp_split, whose numbering scheme IS the exposure id; never for get_images /
exp_psf, whose per-CCD or download numbering would turn tolerated per-CCD
attrition into a whole-exposure hard failure. It is a property of the committed
configs (config_exp_Sp.ini alone carries the entry).
"""

rule exp_get_images:
    output:
        manifest = f"{EXP_DIR}/manifests/exp_get_images.json"
    log:
        f"{EXP_DIR}/logs/exp_get_images.json"
    params:
        pre = lambda wc: unit_pre("exp_get_images", wc.exp,
                                  exp_name=exp_name(wc.exp)),
        script_hash = SCRIPT_HASH
    threads: 1
    retries: 2
    resources:
        mem_mb = lambda wc, attempt: 4000 * attempt,
        runtime = 60
    shell:
        sp_shell("exp_get_images", "config_exp_Gie.ini")

# Split the multi-HDU exposure into single-CCD files (+ headers-*.npy, which the
# tiles' merge_headers reads).
rule exp_split:
    input:
        rules.exp_get_images.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_split.json"
    log:
        f"{EXP_DIR}/logs/exp_split.json"
    params:
        pre = lambda wc: unit_pre("exp_split", wc.exp),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 120
    shell:
        sp_shell("exp_split", "config_exp_Sp.ini")

# SExtractor -> mask_query (FLAG_EXT) -> setools star selection -> PSFEx model
# -> psfex_interp, per CCD.
# setools may reject a sparse CCD (~0.2% attrition) — tolerated by the floor's
# :warn on psfex_interp_runner.
rule exp_psf:
    input:
        rules.exp_split.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_psf.json"
    log:
        f"{EXP_DIR}/logs/exp_psf.json"
    params:
        pre = lambda wc: unit_pre("exp_psf", wc.exp),
        script_hash = SCRIPT_HASH
    threads: 8
    retries: 2
    benchmark:
        # BESIDE manifests/, not inside it: clean_exposure deletes manifests/
        # wholesale, and this tsv is the measured-memory feed for mem_mb sizing
        # (D4). Inside manifests/ it died with the first reclamation and took
        # the campaign's only record of exp_psf's real footprint with it.
        f"{EXP_DIR}/exp_psf.benchmark.tsv"
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        runtime = 240
    shell:
        sp_shell("exp_psf", f"config_exp_{PSF_MODEL}.ini")


# --- persistence (D5) -------------------------------------------------------
# The counterpart of reclamation, and it must come first in the DAG: this packs
# the exposure's keepable PSF products into one tar on the persistent root, and
# clean_exposure below takes its manifest as an input so the store is never
# reclaimed before the keepers have left /scratch. The purge would take them
# anyway — that, not clean_exposure, is what this rule exists for
# (persist_exp.py's docstring argues both halves, and config.yaml's
# `persist_exp:` block carries the keep list and its candidates).
#
# A LOCALRULE (declared in the Snakefile), by exactly the arithmetic that made
# clean_exposure one: the body is a `tar` of a few MB from one shared filesystem
# to another, seconds of work, and one sbatch per exposure would be ~20k
# submissions at DR6 scale for jobs shorter than the scheduling latency. The
# grouping constraint that binds mid-chain localrules (this file's docstring)
# does not bite here: exp_persist's only neighbours are exp_psf, which is too
# heavy to ever fuse, and clean_exposure, which is local itself.
#
# ONE DECLARED OUTPUT, AND IT IS A MANIFEST, NOT THE TAR OR A directory(). The
# tar is not declared: a directory output would attest that a directory exists,
# where what we want written down is WHICH files were packed and how big each was —
# the provenance a rho-statistics run months from now needs in order to know
# what it is reading. The manifest is byte-stable, so a no-op rerun does not
# move its mtime and does not make clean_exposure look out of date.
#
# THE KEEP LIST RIDES ON params. That is the entire reason this is not three
# lines of tar appended to exp_psf's shell: `params` is a rerun trigger, so
# adding a pattern reruns the packing and leaves the PSF chain alone.
rule exp_persist:
    input:
        rules.exp_psf.output.manifest
    output:
        manifest = f"{PROD_EXP_DIR}/manifests/exp_persist.json"
    # No `log:`: the script's only failure modes are "nothing matched" and a
    # name collision, both of which it reports on stderr and neither of which
    # has a per-CCD verdict worth a completeness record.
    params:
        patterns    = " ".join(f"--pattern '{p}'" for p in PERSIST_EXP),
        exp_dir     = lambda wc: exp_dir(wc.exp),
        dest        = lambda wc: f"{prod_exp_dir(wc.exp)}/psf",
        script_hash = PERSIST_HASH
    threads: 1
    retries: 2
    resources:
        mem_mb = 2000,
        runtime = 10
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/persist_exp.py"
        " --exp-dir '{params.exp_dir}' --exp {wildcards.exp}"
        " --dest '{params.dest}' --manifest {output.manifest}"
        " {params.patterns}"


# --- reclamation (D5) -------------------------------------------------------
# The one exception to "no reclamation in this file": clean_exposure OWNS
# exposure-level deletion, and it is a real job, not temp() bookkeeping, because
# an exposure's consumer set closes over the CAMPAIGN. The index supplies that
# set (EXP_TILES, accumulated across invocations); the input is every consuming
# tile's tile_vignets manifest — vignets is the last stage that reads exposure
# products, everything after it reads tile-level files.
#
# What the job deletes, and why a late append still behaves, is argued in
# clean_exposure.py's docstring; params.consumers is what makes a grown consumer
# set stale (same file).
#
# The tile side reads the exposure manifests through ancient() and cuts the
# reclaimed edges of finished tiles (see tile.smk), which is what keeps this
# deletion from rebuilding every neighbouring tile. This rule's OWN inputs are
# deliberately not ancient: a tile that really did rebuild its vignets must
# reschedule the cleans of the exposures it read.
#
# A localrule (declared in the Snakefile). Local execution serialises the cleans
# under local-cores, which costs nothing at rmtree speed and never blocks the
# compute chains (this rule is in none of them).
rule clean_exposure:
    input:
        # ONLY the consumers this invocation may actually build. A consumer that
        # is out of scope had its vignets manifest checked for existence at parse
        # time (clean_targets' eligibility test) — declaring it here as well would
        # pull that finished tile's whole chain into the DAG, where a rebuilt
        # shared exposure then reruns it. That is how one damaged tile reached its
        # spatial neighbours. In-scope consumers keep their edge: they may run in
        # this DAG, so the clean must be ordered after them.
        lambda wc: [tile_manifest(t, "tile_vignets")
                    for t in clean_consumers(wc.exp) if t in READY_SET],
        # The keepers must be off /scratch before the store goes. Unlike the
        # consumer edges above, this edge does not depend on scope: it is the
        # same exposure's own rule, so it drags nothing into the DAG that this
        # exposure's chain did not already put there. It is conditional only on
        # there being a keep list at all — with `persist_exp:` empty, "keep
        # nothing" is a coherent instruction and must not become a dependency on
        # a rule that would fail for having nothing to copy.
        lambda wc: ([prod_exp_manifest(wc.exp, "exp_persist")]
                    if PERSIST_EXP else [])
    output:
        tombstone = f"{EXP_DIR}/cleaned.json"
    params:
        consumers   = lambda wc: ",".join(clean_consumers(wc.exp)),
        script_hash = CLEAN_HASH
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    shell:
        f"python {SCRIPTS}/clean_exposure.py"
        " --exp-dir $(dirname {output.tombstone}) --exp {wildcards.exp}"
        " --tombstone {output.tombstone} --consumers '{params.consumers}'"


# --- the campaign's star catalogue ------------------------------------------
# ONE job per campaign: every exposure's every CCD's `validation_psf-<exp>-<ccd>.fits`,
# stacked into `<products_dir>/full_starcat-0000000.fits`. That file is the
# rho/tau statistics input and sp_validation reads it at exactly that path,
# doing no merging of its own; the old bash chain built it with
# `combine_runs.bash psf` + a `merge_starcat_runner` pass, and the workflow
# emitted neither. The stacking itself is `MergeStarCatPSFEX` — the same class
# the old runner called, reused rather than restated, so a column added to the
# module is a column added here (merge_star_cat.py argues the reuse and the
# tar-member reading).
#
# THE INPUT IS star_cat_inputs() (Snakefile): every exposure of TILES_READY whose
# PSF products are on the persistent root — the live ones through the exp_persist
# manifest edge `rule all` already requests, the RECLAIMED ones through their TAR,
# which no rule declares and which therefore requires nothing to be built. That
# asymmetry is not a flourish; requesting a reclaimed exposure's manifest
# rebuilds its whole chain from VOS, and ancient() does not prevent it (measured
# — the Snakefile carries the numbers). Nothing new enters the DAG either way. It
# is read through an INPUT FUNCTION rather than at module level so that only a
# parse which actually builds this job pays for the walk.
#
# THE PATHS DO NOT REACH THE SHELL, and that is not a style choice: ~20k manifest
# paths is an order of magnitude over Linux's 128 KiB MAX_ARG_STRLEN for a single
# argv entry, so `{input}` here would be a job that dies on exec at DR6 scale.
# The job is handed the two small files the Snakefile itself started from — the
# tile list and the index — and derives THE SAME SET from them; `params.inputs`
# carries that set's FINGERPRINT, which is the rerun trigger. The equality is
# the point: a job that stacked anything the fingerprint did not see would be
# rows no rerun trigger could notice, which is what a glob over products_dir
# would have given on a root shared with an earlier, larger tile list.
# Byte-stable output otherwise (tmp-then-cmp-then-mv), so a no-op rerun does not
# move its mtime.
#
# NOT A LOCALRULE. exp_persist is local because it is 20k jobs of seconds; this
# is one job that holds a campaign's stars in memory (~800k catalogues at DR6
# scale). mem_mb is a guess scaled by attempt, not a measurement — the campaigns
# run so far are 127 exposures, three orders of magnitude short of the case this
# sizing is for, and the first DR6-scale run should replace this number with a
# benchmark.
#
# NO JOB AT ALL when `persist_exp:` keeps no validation catalogue, or when every
# exposure in scope is tombstoned: star_cat_targets() (Snakefile) simply does not
# request the output, and the parse says so rather than a node failing later.
rule star_cat_merge:
    input:
        lambda wc: star_cat_inputs()
    output:
        star_cat = full_starcat()
    params:
        products_dir = str(PRODUCTS_DIR),
        tile_list    = str(config["tile_list"]),
        index_db     = str(INDEX_DB),
        inputs       = unit_fingerprint(star_cat_exposures()),
        script_hash  = MERGE_STAR_HASH
    threads: 1
    resources:
        # Sized on the campaign's own member bytes, slope and intercept
        # measured (the Snakefile's sizing block carries both points, and the
        # ceiling this rule runs into at DR6 scale). Still * attempt, because a
        # measured slope on synthetic tars is not a guarantee about real ones.
        mem_mb = lambda wc, attempt: attempt * (
            STAR_MEM_BASE_MB + STAR_MEM_FACTOR * star_cat_bytes() // 1_000_000),
        # ~2 min per GB of members on the measurement above, doubled, over a
        # floor that covers the fixed cost of opening ~40 members per exposure.
        runtime = lambda wc, attempt: attempt * (
            30 + 4 * star_cat_bytes() // 1_000_000_000)
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/merge_star_cat.py"
        " --products-dir '{params.products_dir}'"
        " --tile-list '{params.tile_list}' --index-db '{params.index_db}'"
        " --output {output.star_cat}"
        f" --psf-model {PSF_MODEL}"
