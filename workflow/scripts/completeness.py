#!/usr/bin/env python3
"""The count-based completeness table — the single failure policy.

Each runner must meet its nominal ``expect`` count unless marked ``warn``;
counts above nominal also pass. On a nibi run of 127 exposures, 64 tiles and
512 ngmix chunks, setools produced 80/80 in every exposure, supporting a
mandatory PSFEx-side ``setools_runner`` count. ``split_exp`` is structurally
all-or-nothing because it raises on an HDU-count mismatch.

A runner below ``expect`` fails its unit unless it has ``warn=True``.
A warning-only shortfall gives the unit status ``warn``; any mandatory
shortfall makes it ``failed``. Scraped error signatures explain shortfalls,
but do not determine whether they pass.

This file is also the ``check`` CLI, called after ShapePipe in each checked
rule's shell line. The rules capture ShapePipe's return code rather than ``&&``-ing
onto it, so the check runs — and the verdict is recorded — even when
``shapepipe_run`` failed::

    rc=0
    shapepipe_run -c $SP_CONFIG/config_exp_Sp.ini -b {threads} || rc=$?
    completeness.py check exp_split {output} --log {log} --job-rc "$rc" || rc=1
    exit $rc

It counts the unit's products under ``$SP_RUN`` and exits nonzero iff a mandatory
runner is below ``expect`` or ``--job-rc`` is nonzero. The verdict combines
the counts and shapepipe_run's own exit status, because a runner can raise after
the counted ones have written their files.

The verdict is written to two files with different jobs:

  * The log (``--log``, the rule's Snakemake ``log:``) records every verdict:
    counts, per-runner detail, scraped reasons and nonzero ``job_rc``.
    It survives failed-job output cleanup for ``run_report.py`` to read.
  * The manifest (the rule's declared ``output:``) receives that verdict only
    on success, including warning-only success. Snakemake removes declared
    outputs of failed jobs; this script does not unlink them.

A surviving ``<stage>.json`` represents success after Snakemake handles failed
outputs. A head-process SIGKILL between failure and output cleanup can leave a
stale success manifest. The profile's ``rerun-incomplete`` handles scheduling;
``run_report.py`` takes the worst status across log and manifest for reporting.

Neither file carries wall-clock time. ``write_if_changed`` preserves unchanged
content and mtime so unrelated ``--forcerun`` invocations do not trigger
downstream reruns.

Per-runner fields:
    expect   nominal file count for a fully complete unit; below it fails
    warn     if True a shortfall warns instead of failing the unit (e.g.
             exposure-side ``psfex_interp_runner``; its tile-side count is
             mandatory)
    subpath  count entries in ``<runner>/output/<subpath>/`` instead of
             ``<runner>/output/`` (e.g. setools split catalogues)

Counts include directory entries as well as files and live symlinks.
Dangling symlinks are excluded.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

# How the tile's galaxy sample is defined: SExtractor on the tile image alone,
# or those detections joined to the UNIONS per-tile catalogue, whose NUMBER
# they take. The run config's `tile_detection:` picks one.
TILE_DETECTIONS = ("sextractor", "unions_catalogue")

# stage -> {runner_subdir: {expect, [warn], [subpath]}}
# exp_psf and tile_vignets are selected by $SP_PSF at check time.
# @sc [decision:per_unit_completeness]
COMPLETENESS = {
    # --- tile prepare (phase A) ---
    # The nibi symlink configs produce one file per INPUT_FILE_PATTERN entry:
    # image+weight for tiles, image+weight+flag for exposures. Tile counts were
    # verified on a nibi run (100 files across 50 tiles).
    "tile_get_images":     {"get_images_runner":      dict(expect=2)},
    "tile_uncompress":     {"uncompress_fits_runner": dict(expect=1)},
    "tile_find_exposures": {"find_exposures_runner":  dict(expect=1)},

    # --- exposure chain ---
    "exp_get_images": {"get_images_runner": dict(expect=3)},
    "exp_split":      {"split_exp_runner":  dict(expect=121)},
    # SExtractor writes three files per CCD (sexcat, background, background_rms),
    # measured at 120 per exposure on nibi. mask_query writes one sexcat_ext
    # per CCD by querying the healsparse map at the detections.
    "exp_psf": {
        "psfex": {
            "sextractor_runner":   dict(expect=120),
            "mask_query_runner":   dict(expect=40),
            "setools_runner":      dict(expect=80, subpath="rand_split"),
            "psfex_runner":         dict(expect=80),
            "psfex_interp_runner":  dict(expect=40, warn=True),
        },
        # MCCD shares the chain up to setools with PSFEx, then fits one
        # focal-plane model per exposure. Preprocessing is a serial runner that
        # merges the 40 CCDs' split catalogues into one training and one test
        # catalogue; fit_val writes the model (fitted_model-<exp>.npy, what the
        # tiles interpolate) and its validation catalogue; merge_starcat and
        # mccd_plots add per-exposure diagnostics. Every runner here is :warn.
        "mccd": {
            "sextractor_runner":          dict(expect=120, warn=True),
            "mask_query_runner":          dict(expect=40, warn=True),
            "setools_runner":             dict(expect=80, warn=True,
                                                subpath="rand_split"),
            # mccd_preprocessing merges the exposure's per-CCD star catalogues
            # into one train and one test catalogue, not one output per CCD.
            "mccd_preprocessing_runner":  dict(expect=2, warn=True),
            # Fit/validation is exposure-wide: one model and one validation
            # catalogue.
            "mccd_fit_val_runner":        dict(expect=2, warn=True),
            "merge_starcat_runner":       dict(expect=1, warn=True),
            # config_exp_mccd enables the ten meanshape and six histogram plots.
            "mccd_plots_runner":          dict(expect=16, warn=True),
        },
        # Image simulations with the true PSF (psf_model: fake): no PSF fit on
        # the exposures, only the SExtractor pass whose background/background_rms
        # checkimages the tile vignets read (config_exp_fake.ini in
        # config/cfis_image_sims). Same 40 CCDs x (sexcat, background, rms).
        "fake": {
            "sextractor_runner":   dict(expect=120),
        },
    },

    # --- tile post ---
    "tile_merge_headers": {"merge_headers_runner": dict(expect=1)},
    # The fetched UNIONS catalogue.
    "tile_get_catalogue": {"get_images_runner":     dict(expect=1)},
    "tile_detect":        {"sextractor_runner":     dict(expect=2)},
    "tile_vignets": {
        "psfex": {
            "psfex_interp_runner":     dict(expect=1),
            "vignetmaker_runner_run_1": dict(expect=1),
            # Five stores per tile: image, weight, flag, background and rms.
            # Each feeds ngmix, so all five are required.
            "vignetmaker_runner_run_2": dict(expect=5),
        },
        # mccd_interp writes galaxy_psf from exposure-wide focal-plane models.
        "mccd": {
            "mccd_interp_runner":        dict(expect=1),
            "vignetmaker_runner_run_1":   dict(expect=1),
            "vignetmaker_runner_run_2":   dict(expect=5),
        },
        # Image simulations: fake_interp_runner writes the same galaxy_psf
        # sqlite psfex_interp_runner writes, from the simulation's PSF dictionary.
        "fake": {
            "fake_interp_runner":       dict(expect=1),
            "vignetmaker_runner_run_1": dict(expect=1),
            "vignetmaker_runner_run_2": dict(expect=5),
        },
    },
    # One check runs inside run_sp_tile_ngmix_Ng${SP_NGMIX_CHUNK}u per chunk,
    # so expect=1 is the correct per-chunk count.
    "tile_ngmix":     {"ngmix_runner":          dict(expect=1)},
    "tile_merge_cats": {"merge_sep_cats_runner": dict(expect=1)},
    "tile_make_cat":  {"make_cat_runner":       dict(expect=1)},
}


def count_products(run_dir, runner, spec):
    """Count entries in ``run_dir/<runner>/output[/<subpath>]/``, excluding dead links.

    scandir, not iterdir: the dirent already says whether an entry is a symlink,
    so only the symlinks need the follow-stat that drops dead links. A plain
    ``p.exists()`` per entry stats every one of them, and at DR6 scale this runs
    once per runner per job over directories of tens to hundreds of files on a
    network filesystem.
    """
    out = run_dir / runner / "output"
    if "subpath" in spec:
        out = out / spec["subpath"]
    if not out.is_dir():
        return 0
    n = 0
    try:
        with os.scandir(out) as entries:
            for e in entries:
                # Only symlinks need a follow-stat to exclude missing targets.
                if not e.is_symlink() or os.path.exists(e.path):
                    n += 1
    except OSError:
        return 0
    return n


def check_counts(stage, run_dir):
    """Return (ok, details); mandatory shortfalls below expect make ok false.

    ``details`` is a list of (runner, n_found, expect, warn) tuples.
    """
    table = COMPLETENESS[stage]
    if stage in ("exp_psf", "tile_vignets"):
        psf_model = os.environ.get("SP_PSF", "psfex")
        try:
            table = table[psf_model]
        except KeyError as exc:
            raise ValueError(
                f"Invalid SP_PSF={psf_model!r}; expected one of {sorted(table)}."
            ) from exc
    details, ok = [], True
    for runner, spec in table.items():
        n = count_products(run_dir, runner, spec)
        warn = spec.get("warn", False)
        details.append((runner, n, spec["expect"], warn))
        if not warn and n < spec["expect"]:
            ok = False
    return ok, details


# --- where a stage writes -------------------------------------------------
#
# stage -> (level, run_sp_<prefix> dir under $SP_RUN/output/). These are the
# committed configs' RUN_NAMEs (RUN_DATETIME=False makes them fixed), so
# the check never resolves a run-log. The ngmix entry interpolates the same env
# var its config does, so chunk K's check looks at chunk K's dir.
#
# Entries correspond to workflow rules. Mask queries run inside exp_psf and
# tile_make_cat, so there is no separate mask stage.
STAGE_DIR = {
    "tile_get_images":     ("tile", "run_sp_tile_Git"),
    "tile_uncompress":     ("tile", "run_sp_tile_Uz"),
    "tile_find_exposures": ("tile", "run_sp_tile_Fe"),
    "exp_get_images":      ("exp",  "run_sp_exp_Gie"),
    "exp_split":           ("exp",  "run_sp_exp_Sp"),
    "exp_psf":             ("exp",  "run_sp_exp_SxSePsf"),
    "tile_merge_headers":  ("tile", "run_sp_tile_Mh_exp"),
    "tile_get_catalogue":  ("tile", "run_sp_tile_Gic"),
    "tile_detect":         ("tile", "run_sp_tile_Sx"),
    "tile_vignets":        ("tile", "run_sp_tile_PiViVi"),
    "tile_ngmix":          ("tile", "run_sp_tile_ngmix_Ng${SP_NGMIX_CHUNK}u"),
    "tile_merge_cats":     ("tile", "run_sp_tile_Ms"),
    "tile_make_cat":       ("tile", "run_sp_tile_Mc"),
}


# --- failure reasons ------------------------------------------------------

# Lines worth showing a human who asks "why is this runner short?". Deliberately
# crude: the point is a pointer into the logs, not a taxonomy (there is no error
# whitelist in this design — the count policy is the policy).
_ERROR_RE = re.compile(
    r"traceback|exception|\berror\b|\bfailed\b|no such file|not found|"
    r"killed|out of memory|oom|segmentation fault|bad chi2",
    re.IGNORECASE)
_TS_RE = re.compile(r"^\d{2}/\d{2}/\d{4} \d{2}:\d{2}:\d{2}\s*")
_NOISE_RE = re.compile(r"A total of 0 errors were recorded")

MAX_LOG_FILES = 40      # logs are per-CCD; a handful is enough to characterise
MAX_TAIL_LINES = 120    # per file
MAX_REASONS = 3         # per runner


def _normalise(line: str) -> str:
    """Collapse a log line to its shape, so 40 per-CCD copies dedupe to one."""
    line = _TS_RE.sub("", line.strip())
    line = re.sub(r"/\S+", "<path>", line)      # paths differ per CCD
    line = re.sub(r"\d+", "N", line)
    return line[:200]


def scrape_reasons(stage_dir, runner):
    """Best-effort, bounded: distinct error-looking lines from a runner's logs.

    Two sources, in order of usefulness: the runner's per-process worker logs
    (``<runner>/logs/process-*.log`` — where the module's own exception lands),
    and the stage's ``logs/log_sp.log`` (where ShapePipe records its error
    tally). Sorted, truncated, deduped by shape — a manifest must stay
    byte-stable for a given tree.
    """
    seen, reasons = {}, []
    candidates = []
    for d in (stage_dir / runner / "logs", stage_dir / "logs"):
        if d.is_dir():
            candidates += sorted(p for p in d.iterdir() if p.is_file())
    for path in candidates[:MAX_LOG_FILES]:
        try:
            lines = path.read_text(errors="replace").splitlines()[-MAX_TAIL_LINES:]
        except OSError:
            continue
        for raw in lines:
            if not _ERROR_RE.search(raw) or _NOISE_RE.search(raw):
                continue
            shape = _normalise(raw)
            if shape in seen:
                seen[shape] += 1
                continue
            seen[shape] = 1
            reasons.append([path.name, _TS_RE.sub("", raw.strip())[:300], shape])
    out = []
    for name, text, shape in reasons[:MAX_REASONS]:
        n = seen[shape]
        out.append(f"{name}: {text}" + (f"  [x{n}]" if n > 1 else ""))
    return out


# --- manifest -------------------------------------------------------------

def build_manifest(stage, run_dir, unit, stage_subdir=None):
    """Count, classify and (on shortfall) scrape. Returns (manifest, ok).

    Stages absent from the table fall back to a zero-output check: any product
    anywhere under the stage dir passes, nothing at all fails.
    """
    level, subdir = STAGE_DIR.get(stage, (None, None))
    subdir = stage_subdir or (os.path.expandvars(subdir) if subdir else None)
    stage_dir = run_dir / "output" / subdir if subdir else run_dir
    manifest = {
        "stage": stage,
        "level": level,
        "unit": unit,
        "run_dir": str(run_dir),
        "stage_dir": str(stage_dir),
        "runners": {},
        "failures": [],
    }

    if stage not in COMPLETENESS:
        produced = list(stage_dir.glob("**/output/*")) if stage_dir.is_dir() else []
        ok = bool(produced)
        manifest["status"] = "complete" if ok else "failed"
        manifest["n_products"] = len(produced)
        if not ok:
            manifest["failures"].append(
                {"runner": None, "found": 0, "expect": 1, "warn": False,
                 "status": "failed", "reasons": [f"zero output under {stage_dir}"]})
        return manifest, ok

    ok, details = check_counts(stage, stage_dir)
    short = False
    for runner, n, expect, warn in details:
        below = n < expect
        if below:
            short = True
        status = "complete" if not below else "warn" if warn else "failed"
        manifest["runners"][runner] = {
            "found": n, "expect": expect, "warn": warn, "status": status,
        }
        if below:
            manifest["failures"].append({
                "runner": runner, "found": n, "expect": expect,
                "warn": warn, "status": status,
                "reasons": scrape_reasons(stage_dir, runner),
            })
    manifest["status"] = "failed" if not ok else ("warn" if short else "complete")
    return manifest, ok


def write_if_changed(path: Path, text: str) -> None:
    """Write only when the bytes differ (see the module docstring on mtime).

    @sc [label:operations] completeness-stable-verdict-publication
    Leave identical content untouched to preserve mtime; publish changes by
    atomic replacement so readers never see a truncated verdict.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        tmp = path.with_name(f".{path.name}.tmp{os.getpid()}")
        try:
            tmp.write_text(text)
            os.replace(tmp, path)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise


def _unit_from_run_dir(run_dir):
    """The human unit ID: the basename of ``$SP_RUN`` (``210.282``, ``2605805``).

    Use the store basename (``210.282``) rather than ShapePipe's dashed
    ``SP_UNIT_NUM`` (``-210-282``).
    """
    return Path(str(run_dir)).name or "unknown"


def main(argv=None) -> int:
    """Run the CLI and persist the per-unit verdict.

    @sc [label:policy] exact-counts-fail-the-unit
    ``--job-rc`` can fail a stage even when product counts pass; the log records
    every verdict, while the manifest is emitted only for success.
    """
    p = argparse.ArgumentParser(description="ShapePipe per-unit completeness check")
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="count products, write the manifest")
    c.add_argument("stage")
    c.add_argument("manifest", type=Path)
    c.add_argument("--log", type=Path, required=True,
                   help="the rule's log: path — the verdict is written here every "
                        "run, success or failure")
    c.add_argument("--run-dir", type=Path, default=None,
                   help="the unit's $SP_RUN (default: the env var)")
    c.add_argument("--unit", default=None,
                   help="override the unit ID (default: basename of $SP_RUN)")
    c.add_argument("--stage-dir", default=None,
                   help="override the run_sp_* subdir (default: the stage table)")
    c.add_argument("--job-rc", type=int, default=0,
                   help="shapepipe_run's exit status, composed into the verdict")
    args = p.parse_args(argv)

    run_dir = args.run_dir or Path(os.environ.get("SP_RUN", ""))
    if not str(run_dir):
        print("[completeness] FATAL: $SP_RUN unset and --run-dir not given",
              file=sys.stderr)
        return 2
    unit = args.unit or _unit_from_run_dir(run_dir)

    manifest, ok = build_manifest(args.stage, Path(run_dir), unit, args.stage_dir)

    # Compose job status with counts (see this function's verdict contract).
    # Include job_rc only on failure; successful verdicts need no extra field.
    if args.job_rc != 0:
        ok = False
        manifest["status"] = "failed"
        manifest["job_rc"] = args.job_rc
        manifest["failures"].append({
            "runner": "shapepipe_run", "found": 0, "expect": 1,
            "warn": False, "status": "failed",
            "reasons": [f"shapepipe_run exited {args.job_rc} "
                        f"(counts met expect)"],
        })

    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"

    # See the function contract and module docstring for log/manifest custody.
    write_if_changed(args.log, text)
    if ok:
        write_if_changed(args.manifest, text)

    for runner, r in manifest["runners"].items():
        tag = {"complete": "OK", "warn": "warn", "failed": "<-- BELOW expect"}
        print(f"[completeness]   {runner}: {r['found']}/{r['expect']} "
              f"{tag[r['status']]}", file=sys.stderr)
    print(f"[completeness] {args.stage} {unit}: {manifest['status']} "
          f"-> {args.log}" + (f" + {args.manifest}" if ok else ""), file=sys.stderr)
    for f in manifest["failures"]:
        for reason in f["reasons"]:
            print(f"[completeness]   {f['runner']}: {reason}", file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
