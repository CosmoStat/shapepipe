"""Exposure chain — per exposure, keyed by exp base id (dedup is structural).

    exp_get_images -> exp_split -> exp_psf -> exp_persist -> exp_maps

Each in the exposure's own sharded work dir, chained by manifests; every config
reads fixed ``$SP_RUN/output/run_sp_exp_*`` INPUT_DIRs, so nothing resolves a
run log. There is no `prepare_exposures` aggregation target: these chains hang
off the compute DAG (`all` <- final_cat <- tile chain <- exposure manifests).

There is no mask-generation rule. ``exp_split`` splits the delivered instrument
flag image per CCD alongside image and weight; SExtractor reads it directly.
Sky-fixed masks are queried within the PSF and tile config chains, not fetched
or rasterized by a separate rule. See ``config/cfis/star_selection.setools``
for star-selection mask cuts and ``config/cfis/config_tile_Mc.ini`` for the
catalogue mask columns.

``exp_persist`` and ``exp_maps`` write to the persistent root. Retention policy
and tar layout are documented in workflow/scripts/persist_exp.py; the exposure
map format is documented in workflow/scripts/exp_maps.py.

NO temp() on exposure products: exposures overlap ~7-10 tiles each, so their
consumer set spans the campaign, not one invocation. ``clean_exposure`` uses
the accumulating index to reclaim them. A temp() would delete an exposure as
soon as this invocation's readers finished, causing destructive reruns of
neighbouring tiles when the tile list grows.

The exposure rules are ungrouped: exp_get_images is a separately retried
download, exp_split has no adjacent short compute rule to fuse with, and
exp_psf is a heavy job (16 GB, 4 h). Group resource composition is documented
in prepare.smk.

NUMBER_LIST ($SP_UNIT_NUM, see unit_num in the Snakefile) is set only for
exp_split, whose numbering scheme is the exposure id. The implementing config
is config_exp_Sp.ini; get_images and exp_psf use different numbering schemes.
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

# PSFEx: SExtractor -> mask_query -> setools -> PSFEx -> psfex_interp, per CCD.
# For the fake PSF, this stage runs SExtractor only. Failure policy lives in
# completeness.py's COMPLETENESS table.
#
# The MCCD branch reserves 2 cores: a measured ~85 min focal-plane fit for
# ~2500 stars uses one thread (8.3 GB), with only 1.75x speedup on 8 BLAS threads.
# The per-CCD stages run 2 wide without reserving 8 cores for the serial fit.
# MCCD campaigns are refused by Snakefile::refuse_unpersistable_psf.
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
    threads: 2 if PSF_MODEL == "mccd" else 8
    retries: 2
    benchmark:
        # Keep the memory-sizing benchmark outside manifests/, which
        # clean_exposure deletes wholesale.
        f"{EXP_DIR}/exp_psf.benchmark.tsv"
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        runtime = 240
    shell:
        sp_shell("exp_psf", f"config_exp_{PSF_MODEL}.ini")


# --- persistence -----------------------------------------------------------
# Packs PSF products into a tar on the persistent root; see persist_exp.py for
# retention, provenance and byte-stable manifest semantics. This is a localrule
# (Snakefile): seconds of packing do not warrant ~20k SLURM submissions at DR6.
# The keep list is a params rerun trigger, so an edit re-packs without refitting.
rule exp_persist:
    input:
        rules.exp_psf.output.manifest
    output:
        manifest = f"{PROD_EXP_DIR}/manifests/exp_persist.json"
    # No `log:`: the script's only failure modes are "nothing matched" and a
    # name collision, both of which it reports on stderr and neither of which
    # has a per-CCD verdict worth a completeness record.
    params:
        # Optional products only; persist_exp.py always packs psf_validation.
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


# --- the exposure's footprint and defects -----------------------------------
# One HealSparse fragment per exposure on the persistent root: the sky pixels
# its valid-PSF CCDs cover, and how many flagged CCD pixels fall in each
# (workflow/scripts/exp_maps.py). It reads the split images' headers and flag
# splits from the scratch store, so clean_exposure waits for it exactly as for
# exp_persist; its one declared input is exp_persist's manifest, which names
# the valid-PSF CCDs. ~8 s per exposure.
# @sc [decision:masking.defect_map_from_flags]
# @sc [decision:masking.nexp_map_valid_psf_ccds]
rule exp_maps:
    input:
        persist = lambda wc: prod_exp_manifest(wc.exp, "exp_persist")
    output:
        manifest = f"{PROD_EXP_DIR}/manifests/exp_maps.json"
    params:
        exp_dir     = lambda wc: exp_dir(wc.exp),
        fragment    = lambda wc: prod_exp_maps(wc.exp),
        script_hash = EXP_MAPS_HASH
    threads: 1
    retries: 2
    resources:
        # 0.15 GB typical; 1.5 GB for a fully flagged CCD (9.4M pixels)
        mem_mb = 3000,
        runtime = 20
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/exp_maps.py"
        " --exp-dir '{params.exp_dir}' --exp {wildcards.exp}"
        " --persist-manifest {input.persist}"
        " --fragment '{params.fragment}' --manifest {output.manifest}"


# --- reclamation -----------------------------------------------------------
# Reclaims exposure stores after their in-scope consumers finish tile_vignets,
# the last tile stage that reads exposure products. Deletion and consumer-set
# staleness are documented in clean_exposure.py; tile.smk owns reclaimed-edge
# handling. Inputs here are not ancient: rebuilt vignets must reschedule clean.
# This is a localrule (Snakefile), serialized under local-cores.
rule clean_exposure:
    input:
        # In-scope consumers only; clean_targets checks out-of-scope vignets
        # at parse time. Declaring them here would pull finished tiles back
        # into the DAG. See workflow/CONTRACTS:
        # clean-exposure-waits-on-persist-iff-psf.
        lambda wc: [tile_manifest(t, "tile_vignets")
                    for t in clean_consumers(wc.exp) if t in READY_SET],
        # Preserve-before-delete ordering, gated by PERSISTS_PSF. durable_edge
        # (Snakefile) uses the manifest for a live store and the durable product
        # for a reclaimed store; the custody contract above governs this edge.
        lambda wc: (durable_edge(wc.exp, "exp_persist", prod_exp_tar(wc.exp))
                    if PERSISTS_PSF else []),
        # The exposure maps' fragment, by the same rule.
        lambda wc: (durable_edge(wc.exp, "exp_maps", prod_exp_maps(wc.exp))
                    if PERSISTS_PSF else [])
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
# One compute job per campaign writes <products_dir>/full_starcat_<run>.hdf5,
# the rho/tau statistics input. merge_star_cat.py owns the reconciled format
# and tar reading; star_cat_inputs/targets (Snakefile) select durable inputs
# and suppress the job when none remain. sp_validation reads this format once
# CosmoStat/sp_validation#340 lands; until then it opens the flat FITS catalogue.
#
# Input paths are DAG edges, not shell arguments: ~20k paths exceed Linux's
# 128 KiB limit for one argv entry. The script derives the same exposure set
# from the tile list and index; params.inputs fingerprints that set for reruns.
# See star_cat_inputs (Snakefile) for the matching-set constraint.
rule star_cat_merge:
    input:
        lambda wc: star_cat_inputs()
    output:
        star_cat = full_starcat()
    params:
        products_dir = str(PRODUCTS_DIR),
        tile_list    = str(config["tile_list"]),
        index_db     = str(INDEX_DB),
        campaign     = CAMPAIGN,
        snapshot     = str(SNAPSHOT_JSON),
        inputs       = unit_fingerprint(star_cat_exposures()),
        script_hash  = MERGE_STAR_HASH
    threads: 1
    resources:
        # Memory scales with the largest exposure, held one at a time;
        # Snakefile's sizing block owns the measurements. Retry scaling allows
        # for real catalogues exceeding the measured synthetic footprint.
        mem_mb = lambda wc, attempt: capped_mem(attempt * (
            STAR_MEM_BASE_MB
            + STAR_MEM_FACTOR * star_cat_max_bytes() // 1_000_000),
            "star_cat_merge"),
        # Runtime scales with total member bytes: ~2 min/GB measured, doubled,
        # plus a floor for opening the per-CCD members.
        runtime = lambda wc, attempt: attempt * (
            30 + 4 * star_cat_bytes() // 1_000_000_000)
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/merge_star_cat.py"
        " --products-dir '{params.products_dir}'"
        " --tile-list '{params.tile_list}' --index-db '{params.index_db}'"
        " --output {output.star_cat}"
        " --campaign '{params.campaign}'"
        " --snapshot-json '{params.snapshot}'"


# --- the campaign's exposure maps -------------------------------------------
# ONE job per campaign: every fragment of the campaign's exposures, summed into
# <products_dir>/nexp_<run>.hsp (exposures with a valid PSF model per sky pixel)
# and nflagged_<run>.hsp (their flagged CCD pixels per sky pixel). Rebuilt
# whole when the exposure set or a fragment changes; inputs, fingerprint and
# the job's own rediscovery of the set follow star_cat_merge. Memory is the two
# maps, 3 MiB per nside-128 coverage pixel the campaign touches: ~13 per
# exposure, capped at the ~23k of the UNIONS footprint (~70 GB at DR6).
# @sc [decision:masking.defect_map_from_flags]
# @sc [decision:masking.nexp_map_valid_psf_ccds]
rule exposure_maps:
    input:
        lambda wc: exposure_maps_inputs()
    output:
        nexp    = nexp_map(),
        nflagged = nflagged_map()
    params:
        products_dir = str(PRODUCTS_DIR),
        tile_list    = str(config["tile_list"]),
        index_db     = str(INDEX_DB),
        inputs       = unit_fingerprint(exposure_maps_exposures()),
        script_hash  = MERGE_MAPS_HASH
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: capped_mem(attempt * (
            1000 + 3.2 * min(13 * len(exposure_maps_exposures()), 23_000)),
            "exposure_maps"),
        runtime = lambda wc, attempt: attempt * (
            30 + len(exposure_maps_exposures()) // 30)
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/merge_exposure_maps.py"
        " --products-dir '{params.products_dir}'"
        " --tile-list '{params.tile_list}' --index-db '{params.index_db}'"
        " --nexp {output.nexp} --nflagged {output.nflagged}"
