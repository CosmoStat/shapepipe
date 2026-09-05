"""Exposure chain — per exposure, keyed by exp base id (dedup is structural).

    exp_get_images -> exp_split ---> exp_mask -> exp_psf -> exp_persist -> exp_footprint
                   -> exp_star_cat -/
    star_catalogue -------------/

``star_catalogue`` is campaign-level, not per-exposure: one fetch of the whole
footprint's stars, which every exposure's ``exp_star_cat`` then cuts locally.

Each in the exposure's own sharded work dir, chained by manifests; every config
reads fixed ``$SP_RUN/output/run_sp_exp_*`` INPUT_DIRs, so nothing resolves a
run log. There is no `prepare_exposures` aggregation target: these chains hang
off the compute DAG (`all` <- final_cat <- tile chain <- exposure manifests).

The last two rules write to the PERSISTENT root, and both exist because the
exposure store is on /scratch and the purge takes it. ``exp_persist`` packs the
PSF products named by `persist_exp:` into one tar per exposure; it is a separate
rule from exp_psf precisely so that editing that list costs a `tar` and not a
four-hour refit (workflow/scripts/persist_exp.py). ``exp_footprint`` then records
where on the sky the CCDs that got a PSF model are, which is the coverage mask's
raw material (workflow/scripts/exp_footprint.py, rules/coverage.smk).

NO temp() anywhere in this file, ever (D5). Exposures overlap tiles by
construction (~7-10 tiles each), so their consumer set closes over the CAMPAIGN,
not over one invocation — reclamation here is clean_exposure's job (S5), driven
by the accumulating index. A temp() here would delete an exposure the moment
this invocation's readers finished and cascade destructive reruns across spatial
neighbours the next time a tile is appended.

GROUPING (``group: "exp_short"``) covers exp_split and exp_mask, and only them —
one sbatch per exposure for two jobs whose medians are 1:28 and 1:54, well under
the 15-minute floor Alliance policy asks us to bundle away. The composition
rules are in prepare.smk's docstring; this chain is linear too, so the group
asks max(mem_mb) = 8000*attempt, max(threads) = 8, sum(runtime) = 240 min.

The two rules NOT in it are structural, not taste:
  * exp_psf is heavy (16 GB, 4 h) and never fuses with a short rule;
  * exp_get_images cannot join, because ``exp_star_cat`` — a LOCALRULE, and so
    ungroupable — sits between it and exp_mask. Pulling get_images in would make
    the group both a dependency and a dependent of exp_star_cat, i.e. a cycle.
    Starting the group at exp_split leaves star_cat's inputs entirely upstream
    of it, so the group has one clean external edge.
Different exposures share no DAG edge, so this is one group job per exposure.

NUMBER_LIST ($SP_UNIT_NUM, see unit_num in the Snakefile) is set only for
exp_split, whose numbering scheme IS the exposure id; never for get_images /
exp_mask / exp_psf, whose per-CCD or download numbering would turn tolerated
per-CCD attrition into a whole-exposure hard failure. It is a property of the
committed configs (config_exp_Sp.ini alone carries the entry).
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

# --- mask star catalogues ---------------------------------------------------
# Two rules, and the split between them is the design: the NETWORK is a function
# of the campaign's sky area, the per-exposure catalogue is a local cut.
#
# `star_catalogue` fetches the footprint's GSC 2.3 stars once, one Vizier query
# per HEALPix chunk, into a run-independent chunk store under config `star_cats`.
# `exp_star_cat` then reads the chunks covering an exposure's focal plane and
# cuts them to it — no network at all. workflow/scripts/star_cats.py holds both
# halves, the geometry they must agree on, and the arithmetic that motivates the
# split; its docstring is the reference for chunking, padding and query counts.

# The container's certifi bundle. The host leaks SSL_CERT_FILE / CURL_CA_BUNDLE
# pointing at a path that does not exist inside the image, so requests is pointed
# at the bundle explicitly (proven in the p3-batch1 bash precedent).
STAR_CAT_CA = "/app/.venv/lib/python3.12/site-packages/certifi/cacert.pem"

# The rules run star_cats.py inside the container (healpy, astroquery, astropy)
# but call apptainer THEMSELVES rather than letting the SDM wrap them
# (`container: None` on both): the CA bundle above and the exposure rule's
# host-side farm loop both need the explicit exec. bin/sp has loaded the
# apptainer module.
#
# WHICH image and WHICH arguments are not this file's to decide, and hand-rolling
# them here was a real divergence: these were the only two rules that ignored a
# user's dev sandbox, because they read `config['container']` — the shared /project
# fallback — instead of the image the Snakefile resolved for everything else.
# `_image` is that resolution (sandbox -> cached SIF -> config), and the profile's
# own apptainer-args are the same string the SDM splices onto every other rule,
# PYTHONPATH pin for shapepipe.utilities.{vizier,cfis} included.
# container.profile_apptainer_args() exists precisely so this file can read them
# rather than restate them. Only the CA bundle is added on top, and only when the
# command touches the network.
def in_container(cmd, *, network=False):
    args = list(_container.profile_apptainer_args())
    if not args:
        # Silently falling back would run these two rules with no --cleanenv and
        # no PYTHONPATH pin, i.e. against a different src/ than every other rule.
        raise WorkflowError(
            f"Could not read apptainer-args from {_container.PROFILE_FILE}; "
            f"star_catalogue and exp_star_cat build their apptainer line from it.")
    if network:
        # The container's certifi bundle, one --env per variable (the profile's
        # own PYTHONPATH entry uses the same one-assignment-per-flag form).
        args += [a for k in ("REQUESTS_CA_BUNDLE", "SSL_CERT_FILE",
                             "CURL_CA_BUNDLE")
                 for a in ("--env", f"{k}={STAR_CAT_CA}")]
    return f"apptainer exec {' '.join(args)} '{_image}' {cmd}"


# The campaign's star catalogue: a first-class durable science product, keyed by
# sky rather than by run. Chunk-need is recomputed from the tile list on every
# run and only the missing chunks are fetched, so appending tiles costs exactly
# the chunks they add.
#
# A LOCALRULE (declared in the Snakefile): it is one job of network I/O, and the
# fetch loop is a 4-wide thread pool inside it — the same modest concurrency the
# per-exposure rule reached by accident through --local-cores, now an explicit
# number that does not scale with the head node's CPU count.
#
# `tile_list_hash` is what makes the incremental behaviour visible to the DAG.
# The tile list is parse-time config, not a rule input (and the profile drops the
# `input` rerun-trigger anyway), so appending tiles would otherwise leave this
# rule up to date against a footprint that has grown. Hashing the list into a
# param reruns it, and the rerun fetches only what is new.
STAR_CAT_MANIFEST = f"{RUN_DIR}/manifests/star_catalogue.json"


rule star_catalogue:
    output:
        manifest = STAR_CAT_MANIFEST
    # No `log:` — see write_manifest() in star_cats.py.
    # `cmd` is a params value, so placeholders in it are not formatted (see
    # unit_pre in the Snakefile). Hence the explicit manifest path.
    params:
        cmd = in_container(
            f"python {SCRIPTS}/star_cats.py fetch"
            f" --tile-list '{config['tile_list']}' --store '{STAR_CATS}'"
            f" --manifest '{STAR_CAT_MANIFEST}'", network=True),
        tile_list_hash = hashlib.md5(
            Path(config["tile_list"]).read_bytes()).hexdigest()[:12],
        script_hash = STAR_CAT_HASH
    container:
        None
    threads: 4
    retries: 2
    resources:
        mem_mb = 4000,
        runtime = 720
    shell:
        "set -euo pipefail\n{params.cmd}"


# The per-exposure catalogue and the 40 per-CCD symlinks the mask module's
# numbering scheme needs. Local: one header read for the focal-plane footprint,
# a load of the chunks covering it, a radial cut.
#
# A LOCALRULE, for the reason the Snakefile's localrules line gives.
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
# chunk store at no network cost.
#
# But a manifest attests FOREVER, and the two things it attests to both live
# outside the unit's manifests/ dir: the cut catalogue on /scratch (60-day purge)
# and the farm itself. Either can vanish under a manifest that still says
# "complete", and then exp_mask runs against nothing. So the ccd-0 farm link is
# declared too — one link stands for all 40, they are created by the same loop
# in the same instant, and declaring 40 buys nothing. Snakemake's existence test
# is os.path.exists, which FOLLOWS symlinks and is therefore False for a link
# whose target the purge removed. A purged cut or a deleted farm makes the rule
# out of date, it reruns, and it re-cuts or re-links as needed.
def star_cat_cmd(exp):
    """The whole rule body, as bash — carried as a params value because it
    contains literal ``{}`` (the manifest JSON); see unit_pre in the Snakefile."""
    cut_dir = f"{STAR_CATS}/exp"
    cat = f"{cut_dir}/star_cat-{exp}.fits"
    work = exp_dir(exp)
    farm = f"{work}/star_cat_exp"
    images = f"{work}/output/run_sp_exp_Gie/get_images_runner/output"
    manifest = exp_manifest(exp, "exp_star_cat")
    body = json.dumps({
        "stage": "exp_star_cat", "level": "exp", "unit": exp,
        "status": "complete", "cat": cat, "link_dir": farm, "n_links": 40,
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
        f"mkdir -p '{cut_dir}' '{farm}' '{work}/manifests'",
        in_container(f"python {SCRIPTS}/star_cats.py cut"
                     f" --images '{images}' --store '{STAR_CATS}'"
                     f" --out '{cat}'"),
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
        rules.exp_get_images.output.manifest,
        # The chunks this cut reads. star_cats.py fails loudly on a chunk that is
        # missing anyway, but the edge is what makes the fetch happen first.
        rules.star_catalogue.output.manifest
    output:
        manifest = f"{EXP_DIR}/manifests/exp_star_cat.json",
        # The sentinel: ccd-0 of the 40-link farm (see above).
        link     = f"{EXP_DIR}/star_cat_exp/star_cat-{{exp}}-0.fits"
    # No `log:`, for the same reason as star_catalogue above.
    params:
        cmd = lambda wc: star_cat_cmd(wc.exp),
        # star_cats.py is external to the shell string, so the `code`
        # rerun-trigger does not see it — same reason SCRIPT_HASH exists.
        script_hash = STAR_CAT_HASH
    container:
        None
    threads: 1
    retries: 2
    resources:
        mem_mb = 4000,
        runtime = 10
    shell:
        "{params.cmd}"

# Split the multi-HDU exposure into single-CCD files (+ headers-*.npy, which the
# tiles' merge_headers reads).
rule exp_split:
    group: "exp_short"
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

rule exp_mask:
    group: "exp_short"
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
        pre = lambda wc: unit_pre("exp_mask", wc.exp),
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
        sp_shell("exp_psf", "config_exp_psfex.ini")


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
# exp_star_cat one: the body is a `tar` of a few MB from one shared filesystem
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


# --- per-CCD sky footprints -------------------------------------------------
# One JSON per exposure recording, for each CCD that HAS a PSF model, the four
# sky corners of that CCD. It is the raw material of the campaign's coverage mask
# (rules/coverage.smk), and it replaces the v1.x chain's VOSpace header download
# + summary-scrape subtraction with a local read of two things this workflow
# already wrote. workflow/scripts/exp_footprint.py argues both inputs and the
# CCD-index invariant tests/unit/test_exp_footprint.py pins.
#
# A LOCALRULE (declared in the Snakefile), by exp_persist's arithmetic and then
# some: one pickle load and ~160 pixel_to_world calls, milliseconds, against a
# scheduling latency of seconds and ~20k exposures at DR6 scale.
#
# ONE DECLARED INPUT, AND IT IS THE PERSIST MANIFEST — deliberately NOT
# exp_psf's as well, even though the WCS array this rule reads is written by
# exp_split and lives in the same scratch store. exp_persist already orders this
# rule after the whole PSF chain, so the second edge would buy no ordering; what
# it WOULD buy is a scratch manifest (the one clean_exposure deletes) in the
# input list of a rule whose output is durable. A persistent-root manifest
# outliving a purged scratch store is a real state, and in it that edge would
# schedule a four-hour VOS rebuild of the exposure to satisfy a few KB of
# provenance. Declaring only the durable input makes the same state a loud
# one-line failure from the script instead.
#
# WRITES TO THE PERSISTENT ROOT, beside exp_persist's manifest and for the same
# reason: the record must outlive both reclamation and the purge — a coverage
# map is campaign-cumulative, so a record written today is read by every map
# built after it.
rule exp_footprint:
    input:
        lambda wc: prod_exp_manifest(wc.exp, "exp_persist")
    output:
        manifest = f"{PROD_EXP_DIR}/manifests/exp_footprint.json"
    # No `log:`, for exp_persist's reason: the failure modes are "no WCS array"
    # and "the two stages disagree about the focal plane", both of which the
    # script reports on stderr and neither of which has a per-CCD verdict.
    params:
        exp_dir     = lambda wc: exp_dir(wc.exp),
        persist     = lambda wc: prod_exp_manifest(wc.exp, "exp_persist"),
        script_hash = FOOTPRINT_HASH
    threads: 1
    retries: 2
    resources:
        mem_mb = 2000,
        runtime = 10
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/exp_footprint.py"
        " --exp-dir '{params.exp_dir}' --exp {wildcards.exp}"
        " --persist-manifest '{params.persist}'"
        " --manifest {output.manifest}"


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
                    if PERSIST_EXP else []),
        # And the footprint, for the same ordering reason one layer further out:
        # it is derived from headers-<exp>.npy, which lives in the store this job
        # deletes. Reclamation must not overtake the read, and unlike the purge
        # this deletion is ours to order. Conditional on the keep list naming the
        # PSF products (Snakefile's PERSIST_HAS_PSF) because that is exactly when
        # the footprint rule exists at all.
        lambda wc: ([prod_exp_manifest(wc.exp, "exp_footprint")]
                    if PERSIST_HAS_PSF else [])
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
