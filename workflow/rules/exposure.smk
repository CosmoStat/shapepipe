"""Exposure chain — per exposure, keyed by exp base id (dedup is structural).

    exp_get_images -> exp_split -> exp_mask -> exp_psf

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

# Split the multi-HDU exposure into single-CCD files (+ headers-*.npy, which the
# tiles' merge_headers reads).
rule exp_split:
    input:
        rules.exp_get_images.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_split.json"
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
        rules.exp_split.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_mask.json"
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
