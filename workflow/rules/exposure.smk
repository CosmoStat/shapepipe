"""Exposure chain — per exposure, keyed by exp base id (dedup is structural).

    exp_get_images -> exp_split ---> exp_mask -> exp_psf
                   -> exp_star_cat -/

Each in the exposure's own sharded work dir, chained by manifests; every config
reads fixed ``$SP_RUN/output/run_sp_exp_*`` INPUT_DIRs, so nothing resolves a
run log. There is no `prepare_exposures` aggregation target: these chains hang
off the compute DAG (`all` <- final_cat <- tile chain <- exposure manifests).

NO temp() anywhere in this file, ever (D5). Exposures overlap tiles by
construction (~7-10 tiles each), so their consumer set closes over the CAMPAIGN,
not over one invocation — reclamation here is clean_exposure's job (S5), driven
by the accumulating index. A temp() here would delete an exposure the moment
this invocation's readers finished and cascade destructive reruns across spatial
neighbours the next time a tile is appended.

NUMBER_LIST is set only for exp_split (its numbering scheme IS the exposure id);
never for get_images / exp_mask / exp_psf, whose per-CCD or download numbering
would make the #746 startup validation turn tolerated per-CCD attrition into a
whole-exposure hard failure. That is now a property of the committed configs
(config_exp_Sp.ini has NUMBER_LIST = $SP_UNIT_NUM; Gie/Ma/psfex have none).
"""

rule exp_get_images:
    output:
        manifest = f"{EXP_DIR}/manifests/exp_get_images.json"
    log:
        f"{EXP_DIR}/logs/exp_get_images.json"
    params:
        pre = lambda wc: unit_pre("exp_get_images", "exp", wc.exp,
                                  exp_name=EXP[wc.exp]),
        script_hash = SCRIPT_HASH
    threads: 1
    retries: 2
    resources:
        mem_mb = lambda wc, attempt: 4000 * attempt,
        runtime = 60
    shell:
        sp_shell("exp_get_images", "config_exp_Gie.ini")

# --- mask star catalogues ---------------------------------------------------
# One GSC 2.3 (Vizier) cone query per exposure, covering the whole MegaCam focal
# plane, plus the 40 per-CCD symlinks the mask module's numbering scheme needs.
#
# It reads the exp_get_images output dir directly: that is a symlink farm of
# image-<exp>.fitsfz, exactly what create_star_cat.py's `-k exp` mode consumes —
# one multi-extension header read, one query, ~6 s measured.
#
# A LOCALRULE (declared in the Snakefile), so it runs in the head process — an
# sbatch job on a compute node, and nibi compute nodes have internet (verified by
# curl to vizier.cds.unistra.fr from a compute job). Local execution does NOT
# serialise the queries — an earlier version of this comment claimed it did, and
# that was false. In cluster mode snakemake runs local jobs in a worker pool
# sized by `--local-cores`, and this rule takes one thread each, so as many CDS
# queries (and apptainer execs) run at once as that pool is wide. The nibi
# profile leaves --local-cores unset, so it is the head process's CPU count —
# 5 in the measured allocation.
#
# The localrule stays anyway. 5-way concurrency against Vizier is ordinary
# politeness, not abuse, and the alternative is one sbatch per exposure: ~20k
# submissions for six seconds of work each, which is exactly the sub-15-minute
# scheduler churn cluster policy asks us to bundle away. If CDS ever objects,
# --local-cores is the knob — it caps these queries directly.
#
# The CATALOGUE is cached run-independently under config `star_cats`, keyed by
# exposure. create_star_cat.py skips a cat that exists, so retries, reruns and
# later campaigns re-link an existing file and issue no query.
#
# The per-unit farm is a REAL directory holding exactly this exposure's 40
# numbers, and that is load-bearing: config_exp_Ma.ini reads it as an INPUT_DIR
# and the file handler INTERSECTS the numbers found across INPUT_DIRs, so a
# symlink to a shared whole-store pool contributes every other exposure's numbers
# and the intersection is empty ("numbers ... do not intersect", live).
#
# TWO declared outputs, and the second one is the point.
#
# The manifest keeps the "one rule, one manifest" currency of every other rule:
# written last, unique to this rule, a record of what the farm points at, and
# deleted by clean_exposure so a reclaimed exposure rebuilds its farm from the
# cache at no query cost.
#
# But a manifest attests FOREVER, and the two things it attests to both live
# outside the unit's manifests/ dir: the cache FITS on /scratch (60-day purge)
# and the farm itself. Either can vanish under a manifest that still says
# "complete", and then exp_mask runs against nothing. So the ccd-0 farm link is
# declared too — one link stands for all 40, they are created by the same loop
# in the same instant, and declaring 40 buys nothing. Snakemake's existence test
# is os.path.exists, which FOLLOWS symlinks and is therefore False for a link
# whose target the purge removed. A purged cache or a deleted farm makes the
# rule out of date, it reruns, and it re-queries or re-links as needed — the
# recovery config.yaml's star_cats comment promises, now actually wired up.

# The container's certifi bundle. The host leaks SSL_CERT_FILE / CURL_CA_BUNDLE
# pointing at a path that does not exist inside the image, so requests is pointed
# at the bundle explicitly (proven in the p3-batch1 bash precedent).
STAR_CAT_CA = "/app/.venv/lib/python3.12/site-packages/certifi/cacert.pem"


def star_cat_cmd(exp):
    """The whole rule body, as bash — carried as a params value, never inlined
    in ``shell:``: it contains literal ``{}`` (the manifest JSON) and snakemake
    formats a shell string once, which would consume those braces."""
    cache = f"{STAR_CATS}/exp"
    cat = f"{cache}/star_cat-{exp}.fits"
    work = exp_dir(exp)
    farm = f"{work}/star_cat_exp"
    images = f"{work}/output/run_sp_exp_Gie/get_images_runner/output"
    manifest = exp_manifest(exp, "exp_star_cat")
    body = json.dumps({
        "stage": "exp_star_cat", "level": "exp", "unit": exp,
        "status": "complete", "cache": cat, "link_dir": farm, "n_links": 40,
    }, indent=2, sort_keys=True)
    return "\n".join([
        "set -euo pipefail",
        # LEGACY-SYMLINK HAZARD. Unit dirs built before this rule existed carry
        # star_cat_exp as a SYMLINK into the old shared star-cat pool. `mkdir -p`
        # is a no-op on an existing symlink-to-directory, so the 40-link loop
        # below followed it and wrote this exposure's links INTO THE SHARED POOL
        # (520 stray links found live). Replace the link — never `rm -rf` it,
        # which would recurse into the pool, and never touch a real directory:
        # a real farm is this rule's own output and `ln -sfn` refreshes it.
        f"[ -L '{farm}' ] && rm -f '{farm}' || true",
        f"mkdir -p '{cache}' '{farm}' '{work}/manifests'",
        # Same binds and isolation as the profile's apptainer-args; the SDM does
        # not wrap this rule (container: None below), because the farm loop and
        # the manifest write belong to the host side of the job.
        f"apptainer exec --cleanenv --home {Path.home()} --bind /project --bind /scratch"
        f" --env PYTHONPATH={REPO_DIR}/src,REQUESTS_CA_BUNDLE={STAR_CAT_CA}"
        f",SSL_CERT_FILE={STAR_CAT_CA},CURL_CA_BUNDLE={STAR_CAT_CA}"
        f" '{config['container']}'"
        f" python {REPO_DIR}/scripts/python/create_star_cat.py"
        f" -i '{images}' -o '{cache}' -k exp",
        f"test -s '{cat}'",
        # The fan-out the file handler's NUMBERING_SCHEME wants: 40 links to the
        # one focal-plane catalogue (pattern from the p3-batch1 precedent).
        f"for ccd in $(seq 0 39); do ln -sfn '{cat}' "
        f"'{farm}/star_cat-{exp}-'\"$ccd\"'.fits'; done",
        # Byte-stable, and written only after the links exist: an unconditional
        # write would move the mtime, which is a rerun-trigger.
        f"tmp='{manifest}.tmp'",
        "cat > \"$tmp\" <<'SP_STAR_CAT_JSON'",
        body,
        "SP_STAR_CAT_JSON",
        f"cmp -s \"$tmp\" '{manifest}' && rm -f \"$tmp\" || mv -f \"$tmp\" '{manifest}'",
    ])


rule exp_star_cat:
    input:
        rules.exp_get_images.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_star_cat.json",
        # The sentinel: ccd-0 of the 40-link farm (see above).
        link     = f"{EXP_DIR}/star_cat_exp/star_cat-{{exp}}-0.fits"
    # No `log:`. Every other rule's log carries a completeness VERDICT, and this
    # rule computes none: it is not a shapepipe_run, it has no count floor, and
    # under `set -euo pipefail` it either completes or aborts at the failing
    # step. Snakemake's own captured job stderr is the evidence for that.
    params:
        cmd = lambda wc: star_cat_cmd(wc.exp),
        # create_star_cat.py is external to the shell string, so the `code`
        # rerun-trigger does not see it — same reason SCRIPT_HASH exists.
        script_hash = STAR_CAT_HASH
    # No SDM wrapping: the rule calls apptainer itself, because only the query
    # runs in the container (bin/sp has loaded the apptainer module).
    container:
        None
    threads: 1
    retries: 2
    resources:
        mem_mb = 2000,
        runtime = 10
    shell:
        "{params.cmd}"

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
        pre = lambda wc: unit_pre("exp_split", "exp", wc.exp),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 120
    shell:
        sp_shell("exp_split", "config_exp_Sp.ini")

rule exp_mask:
    input:
        # Both inputs are real INPUT_DIRs of config_exp_Ma.ini: the split CCDs
        # and this exposure's own star_cat_exp farm.
        rules.exp_split.output.manifest,
        rules.exp_star_cat.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_mask.json"
    log:
        f"{EXP_DIR}/logs/exp_mask.json"
    params:
        pre = lambda wc: unit_pre("exp_mask", "exp", wc.exp),
        script_hash = SCRIPT_HASH
    threads: 4
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 120
    shell:
        sp_shell("exp_mask", "config_exp_Ma.ini")

# SExtractor -> setools star selection -> PSFEx model -> psfex_interp, per CCD.
# setools may reject a sparse CCD (~0.2% attrition) — tolerated by the floor's
# :warn on psfex_interp_runner.
rule exp_psf:
    input:
        rules.exp_mask.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_psf.json"
    log:
        f"{EXP_DIR}/logs/exp_psf.json"
    params:
        pre = lambda wc: unit_pre("exp_psf", "exp", wc.exp),
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
        sp_shell("exp_psf", "config_exp_psfex.ini")


# --- reclamation (D5) -------------------------------------------------------
# The one exception to "no reclamation in this file": clean_exposure OWNS
# exposure-level deletion, and it is a real job, not temp() bookkeeping, because
# an exposure's consumer set closes over the CAMPAIGN. The index supplies that
# set (EXP_TILES, accumulated across invocations); the input is every consuming
# tile's tile_vignets manifest — vignets is the last stage that reads exposure
# products, everything after it reads tile-level files.
#
# Three properties make the late append behave (see clean_exposure.py):
#   * the job deletes the exposure's manifests too, so a tile appended after the
#     clean sees an unbuilt chain and regenerates it instead of running against
#     an empty store. The tombstone deliberately does NOT stand in for those
#     manifests — it is not an input to anything but itself.
#   * a finished tile is not disturbed: Snakemake demands a missing intermediate
#     only when something downstream of it must run.
#   * params.consumers carries the consumer set, so growing it makes the
#     tombstone stale under the default `params` rerun-trigger; the clean job
#     reruns after the new tile's vignets, against the enlarged set.
#
# The tile side reads the exposure manifests through ancient() (see tile.smk),
# which is what keeps this deletion from rebuilding every neighbouring tile.
# This rule's OWN inputs are deliberately not ancient: a tile that really did
# rebuild its vignets must reschedule the cleans of the exposures it read.
#
# A localrule (declared in the Snakefile): it is an rmtree, not science, and one
# sbatch per exposure would be ~20k scheduler submissions at DR6 scale. Local
# execution serialises them under local-cores, which costs nothing at rmtree
# speed and never blocks the compute chains (this rule is in none of them).
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
                    for t in clean_consumers(wc.exp) if t in READY_SET]
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
