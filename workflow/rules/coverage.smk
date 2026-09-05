"""Coverage — the campaign's HealSparse nexp mask, from the exposure footprints.

    exp_footprint (per exposure, exposure.smk) -> coverage_map (per campaign)

One rule, and it is the only campaign-level product besides the report. The
precedent it follows is ``star_catalogue`` (exposure.smk): a rule keyed by the
campaign rather than by a unit, whose inputs are the units' own records.

CAMPAIGN-CUMULATIVE. The declared inputs are the IN-SCOPE, non-tombstoned
footprint manifests — ordering and rerun semantics for free, without dragging
out-of-scope tiles into the DAG through a new target. The SCRIPT then reads every
footprint record on the persistent root, tombstoned exposures included: their
records outlive their scratch stores and are still valid sky. So appending tiles
grows the map instead of replacing it, which is what a survey coverage mask
should do, and rebuilding is one job rather than a campaign.

NOT A LOCALRULE, unlike every other rule that writes to the persistent root: at
DR6 scale this stamps ~1M polygons at nside=131072 in a Python loop
(coverage_map_builder.build_map). The resources below are a first sizing from
that count and not a measurement — the polygon loop is unmeasured above a few
thousand CCDs.

PLOTS STAY OUT OF THE DAG. `plot_coverage_map -i <products_dir>/coverage/
coverage.hsp ...`, by hand, with the windows in config.yaml's `coverage.plot`
block — the same argument that keeps run_report.py a standalone script: it is a
human act on a durable product.
"""

# `coverage:` is OFF by default: the map is a campaign-end product, and a
# half-finished campaign's map is a picture of how far it has got rather than of
# the survey. flag() because `--config` delivers booleans as strings.
COVERAGE = config.get("coverage") or {}
COVERAGE_ENABLED = flag(COVERAGE.get("enabled", False))
NSIDE_COVERAGE = int(COVERAGE.get("nside_coverage", 128))
NSIDE = int(COVERAGE.get("nside", 131072))

COVERAGE_DIR = f"{PRODUCTS_DIR}/coverage"
COVERAGE_HSP = f"{COVERAGE_DIR}/coverage.hsp"
COVERAGE_MANIFEST = f"{COVERAGE_DIR}/manifests/coverage_map.json"
COVERAGE_HASH = script_hash("coverage_map.py")

# PARSE-TIME GUARDS, and both of them fail before submission rather than after.
#
# The first is the precondition the whole chain rests on: exp_footprint reads its
# valid-PSF CCD set off exp_persist.json's members, so with `validation_psf-*`
# absent from `persist_exp:` no footprint job exists, no record exists, and
# `coverage:` on would otherwise mean a green run producing nothing — or, worse,
# a map built from whatever stale records happened to be on the products root.
if COVERAGE_ENABLED and not PERSIST_HAS_PSF:
    raise WorkflowError(
        f"coverage: is enabled but persist_exp: {PERSIST_EXP} packs no "
        f"{PSF_PRODUCT_PROBE}-like files. The coverage map is built from the "
        f"valid-PSF CCD set, which is read off exp_persist's tar members — "
        f"without that pattern there is no such set. Add it, or turn coverage "
        f"off.")

# The second: healsparse requires powers of two, and it says so only once a job
# has reached a node.
for _key, _n in (("nside_coverage", NSIDE_COVERAGE), ("nside", NSIDE)):
    if COVERAGE_ENABLED and (_n <= 0 or _n & (_n - 1)):
        raise WorkflowError(f"coverage.{_key} must be a power of 2, got {_n}")


def coverage_targets():
    """The campaign map, when `coverage:` asks for one.

    HEAD PROCESS ONLY, for the reason clean_targets() gives: this feeds `rule
    all` at module level, so it is evaluated on every per-job re-parse too, none
    of which can schedule `all`.
    """
    if not COVERAGE_ENABLED or not workflow.is_main_process:
        return []
    return [COVERAGE_HSP, COVERAGE_MANIFEST]


rule coverage_map:
    input:
        footprint_targets()
    output:
        hsp      = COVERAGE_HSP,
        manifest = COVERAGE_MANIFEST
    # No `log:`: this is one job with one verdict, and its stderr is the job's.
    params:
        products       = PRODUCTS_DIR,
        nside_coverage = NSIDE_COVERAGE,
        nside          = NSIDE,
        script_hash    = COVERAGE_HASH
    threads: 1
    resources:
        mem_mb = 32000,
        runtime = 240
    shell:
        "set -euo pipefail\n"
        f"python {SCRIPTS}/coverage_map.py"
        " --products-dir '{params.products}'"
        " --out {output.hsp} --manifest {output.manifest}"
        " --nside-coverage {params.nside_coverage} --nside {params.nside}"
