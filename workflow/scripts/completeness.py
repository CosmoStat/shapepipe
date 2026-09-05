#!/usr/bin/env python3
"""The count-floor completeness table — the single failure policy.

This is the ported ``complete_check`` count table from the v2.0 bash layer
(``run_job_sp_canfar_v2.0.bash`` job dispatch, survey §4). It is the *only*
failure policy in the design: there is no 3-class taxonomy and no error-signature
whitelist. A stage is a real failure iff a mandatory runner produced fewer than
its ``floor`` files; per-CCD attrition (a sparse CCD setools rejects, ~0.2%)
sits between ``floor`` and ``expect`` and is tolerated.

This file is also the ``check`` CLI, the second half of every rule's shell line
(PRD D2/D3). The rules capture ShapePipe's return code rather than ``&&``-ing
onto it, so the check runs — and the verdict is recorded — even when
``shapepipe_run`` failed::

    rc=0
    shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini -b {threads} || rc=$?
    completeness.py check exp_mask {output} --log {log} --job-rc "$rc" || rc=1
    exit $rc

It counts the unit's products under ``$SP_RUN`` and exits nonzero iff a mandatory
runner is below its floor OR ``--job-rc`` is nonzero — the verdict is COMPOSED of
the counts and shapepipe_run's own exit status, because a runner can raise after
the counted ones have written their files.

The verdict is written to two files with different
different jobs:

  * the LOG (``--log``, the rule's snakemake ``log:``) gets the full verdict on
    EVERY run, success or failure — counts against floors, per-runner detail,
    scraped failure reasons, ``job_rc`` when nonzero. Snakemake never deletes a
    log file, so it survives the failed job that wrote it and is the post-mortem
    evidence ``run_report.py`` reads.
  * the MANIFEST (the rule's declared ``output:``) gets that same verdict ONLY
    when it is a success. Snakemake deletes a failed job's declared output
    natively, so nothing here has to unlink anything.

``<stage>.json`` therefore means "this stage succeeded": a resume cannot schedule
a downstream stage on top of a failed one. The removal of a PREVIOUS success's
manifest is snakemake's to do, not this script's, and the one gap that leaves is
a head process SIGKILLed between the job's failure and that deletion — a
success-named manifest then outlives the failure it no longer describes. The
profile's ``rerun-incomplete`` covers the DAG side, and ``run_report.py`` takes
the WORST status across a stage's log and manifest, so the report is right even
in that window.

Neither file carries wall-clock, and each is rewritten ONLY when its content
changes — identical on-disk state must leave a byte-identical manifest with an
UNMOVED mtime, or the mtime rerun-trigger churns the cone on every unrelated
``--forcerun``.

Per-runner fields:
    expect   nominal file count for a fully complete unit (report yardstick)
    floor    the fail-loud minimum (below this the job exits nonzero)
    warn     if True the runner never fails the unit at all (bash ``:warn`` —
             e.g. psfex_interp on tiles missing some epochs)
    subpath  count files in ``<runner>/output/<subpath>/`` instead of
             ``<runner>/output/`` (bash ``:rand_split`` — setools split cats)

Counts are file counts in the runner's output dir, matching the bash
``ls <out_dir>/ | wc -l`` semantics (broken symlinks excluded by the caller).
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

# stage -> {runner_subdir: {expect, floor, [warn], [subpath]}}
# exp_psf and tile_vignets are selected by $SP_PSF at check time.
COMPLETENESS = {
    # --- tile prepare (phase A) ---
    # get_images counts are CONFIG-FLAVOR-DEPENDENT: the v2.0 bash table said 4/6
    # for the canfar vos flavor; the nibi symlink configs produce one file per
    # INPUT_FILE_PATTERN entry (tile: image+weight=2; exp: image+weight+flag=3),
    # verified against the p3-batch1 baseline tree (100 files / 50 tiles).
    "tile_get_images":     {"get_images_runner":      dict(expect=2, floor=2)},
    "tile_uncompress":     {"uncompress_fits_runner": dict(expect=1, floor=1)},
    "tile_find_exposures": {"find_exposures_runner":  dict(expect=1, floor=1)},

    # --- exposure chain ---
    "exp_get_images": {"get_images_runner": dict(expect=3, floor=3)},
    "exp_split":      {"split_exp_runner":  dict(expect=121, floor=41)},
    "exp_mask":       {"mask_runner":       dict(expect=40, floor=1)},
    # sextractor expect is nibi-flavor: 3 files/CCD (sexcat + background +
    # background_rms; v2.0's 80 assumed 2/CCD), verified against the P0 tree
    # AND the bash baseline (both 120/exposure).
    "exp_psf": {
        "psfex": {
            "sextractor_runner":   dict(expect=120, floor=2),
            "setools_runner":      dict(expect=80, floor=2, subpath="rand_split"),
            "psfex_runner":         dict(expect=80, floor=2),
            "psfex_interp_runner":  dict(expect=40, floor=0, warn=True),
        },
        # MCCD counts are derived from config_exp_mccd.ini and its per-CCD
        # runners, but this chain has not been exercised through this workflow.
        # Keep the expected counts visible while making the unverified branch
        # warning-only until a real campaign validates its floors.
        "mccd": {
            "sextractor_runner":          dict(expect=120, floor=0, warn=True),
            "setools_runner":             dict(expect=80, floor=0, warn=True,
                                                subpath="rand_split"),
            "mccd_preprocessing_runner":  dict(expect=80, floor=0, warn=True),
            # Fit/validation is exposure-wide: one model and one validation
            # catalogue, unlike the per-CCD preprocessing outputs.
            "mccd_fit_val_runner":        dict(expect=2, floor=0, warn=True),
            "merge_starcat_runner":       dict(expect=1, floor=0, warn=True),
            # config_exp_mccd enables the ten meanshape and six histogram plots.
            "mccd_plots_runner":          dict(expect=16, floor=0, warn=True),
        },
    },

    # --- tile post ---
    "tile_merge_headers": {"merge_headers_runner": dict(expect=1, floor=1)},
    "tile_detect":        {"sextractor_runner":    dict(expect=2, floor=2)},
    "tile_vignets": {
        "psfex": {
            "psfex_interp_runner":     dict(expect=1, floor=1),
            "vignetmaker_runner_run_1": dict(expect=1, floor=1),
            # 5 sqlites/tile on nibi (image/weight/flag/background/background_rms);
            # v2.0's 4 was the canfar flavor. floor follows the tile-post pattern
            # (expect=floor: all-or-nothing, every vignette feeds ngmix).
            "vignetmaker_runner_run_2": dict(expect=5, floor=5),
        },
        # MCCD is wired but unvalidated here; retain the expected runner names
        # and counts as warnings until a workflow campaign exercises them.
        "mccd": {
            "mccd_interp_runner":        dict(expect=1, floor=0, warn=True),
            "vignetmaker_runner_run_1":   dict(expect=1, floor=0, warn=True),
            "vignetmaker_runner_run_2":   dict(expect=5, floor=0, warn=True),
        },
    },
    "tile_ngmix":     {"ngmix_runner":          dict(expect=1, floor=1)},
    "tile_merge_cats": {"merge_sep_cats_runner": dict(expect=1, floor=1)},
    "tile_make_cat":  {"make_cat_runner":       dict(expect=1, floor=1)},
}


def count_products(run_dir, runner, spec):
    """Count files in ``run_dir/<runner>/output[/<subpath>]/`` (live links only).

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
                # A dead symlink must not count (the bash
                # `ls | wc -l` semantics this ports counted live files only), and
                # a symlink is the only entry that can be dead — so it is the
                # only one worth a follow-stat.
                if not e.is_symlink() or os.path.exists(e.path):
                    n += 1
    except OSError:
        return 0
    return n


def check_floor(stage, run_dir):
    """Return (ok, details). ok is False iff a mandatory runner is below floor.

    ``details`` is a list of (runner, n_found, floor, expect, warn) tuples.
    """
    table = COMPLETENESS[stage]
    if stage in ("exp_psf", "tile_vignets"):
        psf_model = os.environ.get("SP_PSF", "psfex")
        try:
            table = table[psf_model]
        except KeyError as exc:
            raise ValueError(
                f"Invalid SP_PSF={psf_model!r}; expected one of psfex, mccd."
            ) from exc
    details, ok = [], True
    for runner, spec in table.items():
        n = count_products(run_dir, runner, spec)
        details.append((runner, n, spec["floor"], spec["expect"],
                        spec.get("warn", False)))
        if not spec.get("warn", False) and n < spec["floor"]:
            ok = False
    return ok, details


# --- where a stage writes -------------------------------------------------
#
# stage -> (level, run_sp_<prefix> dir under $SP_RUN/output/). These are the
# committed configs' RUN_NAMEs (RUN_DATETIME=False makes them fixed, PRD D2), so
# the check never resolves a run-log. The ngmix entry interpolates the same env
# var its config does, so chunk K's check looks at chunk K's dir.
#
# EVERY ENTRY HERE (and in COMPLETENESS above) HAS A RULE. The table used to
# carry two stages that did not: `tile_mask` (run_sp_tile_Ma, mask_runner 1/1)
# and `tile_detect_uc` (run_sp_tile_Uc). The committed config chain is the
# "sx_nomask" tile_detect variant and no tile-mask config was committed, so both
# were unreachable — tile.smk's docstring is where the masked variant is argued,
# and it is a config plus a rule plus these two rows, added back together.
STAGE_DIR = {
    "tile_get_images":     ("tile", "run_sp_tile_Git"),
    "tile_uncompress":     ("tile", "run_sp_tile_Uz"),
    "tile_find_exposures": ("tile", "run_sp_tile_Fe"),
    "exp_get_images":      ("exp",  "run_sp_exp_Gie"),
    "exp_split":           ("exp",  "run_sp_exp_Sp"),
    "exp_mask":            ("exp",  "run_sp_exp_Ma"),
    "exp_psf":             ("exp",  "run_sp_exp_SxSePsf"),
    "tile_merge_headers":  ("tile", "run_sp_tile_Mh_exp"),
    "tile_detect":         ("tile", "run_sp_tile_Sx"),
    "tile_vignets":        ("tile", "run_sp_tile_PiViVi"),
    "tile_ngmix":          ("tile", "run_sp_tile_ngmix_Ng${SP_NGMIX_CHUNK}u"),
    "tile_merge_cats":     ("tile", "run_sp_tile_Ms"),
    "tile_make_cat":       ("tile", "run_sp_tile_Mc"),
}


# --- failure reasons ------------------------------------------------------

# Lines worth showing a human who asks "why is this runner short?". Deliberately
# crude: the point is a pointer into the logs, not a taxonomy (there is no error
# whitelist in this design — the count floor is the policy).
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

    Stages absent from the table fall back to the zero-output floor: any product
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
                {"runner": None, "found": 0, "floor": 1,
                 "reasons": [f"zero output under {stage_dir}"]})
        return manifest, ok

    ok, details = check_floor(stage, stage_dir)
    short = False
    for runner, n, floor, expect, warn in details:
        below = n < floor
        if n < expect:
            short = True
        manifest["runners"][runner] = {
            "found": n, "expect": expect, "floor": floor, "warn": warn,
            "status": ("complete" if n >= expect else
                       "warn" if (warn or not below) else "below_floor"),
        }
        if below:
            manifest["failures"].append({
                "runner": runner, "found": n, "floor": floor, "expect": expect,
                "warn": warn, "reasons": scrape_reasons(stage_dir, runner),
            })
    manifest["status"] = "failed" if not ok else ("warn" if short else "complete")
    return manifest, ok


def write_if_changed(path: Path, text: str) -> None:
    """Write only when the bytes differ — see the module docstring on mtime."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def _unit_from_run_dir(run_dir):
    """The human unit ID: the basename of ``$SP_RUN`` (``210.282``, ``2605805``).

    NOT ``SP_UNIT_NUM``, which carries ShapePipe's dashed numbering form
    (``-210-282``) and would put ``210-282`` in the manifest — a key that joins
    to nothing. ``run_report`` keys units on the store directory name, which is
    exactly this basename, so the two now agree.
    """
    return Path(str(run_dir)).name or "unknown"


def main(argv=None) -> int:
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

    # The verdict is COMPOSED of two independent statements: the count floors
    # (above) and shapepipe_run's own exit status (here). Counts alone are not
    # enough — a runner can raise AFTER the counted runners have written their
    # files, so the floors are met while the job died. Without this, such a job
    # publishes a SUCCESS manifest that snakemake then deletes as a failed job's
    # output — and the log would claim "complete" for a stage with no manifest,
    # which reads as a bookkeeping bug rather than as the failure it is.
    #
    # Recorded only when nonzero, which keeps every existing success manifest
    # byte-identical (the mtime rerun-trigger reads those bytes).
    if args.job_rc != 0:
        ok = False
        manifest["status"] = "failed"
        manifest["job_rc"] = args.job_rc
        manifest["failures"].append({
            "runner": "shapepipe_run", "found": 0, "floor": 1, "expect": 1,
            "warn": False,
            "reasons": [f"shapepipe_run exited {args.job_rc} "
                        f"(counts above floor)"],
        })

    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"

    # The log ALWAYS gets the verdict; the manifest gets it only on success. No
    # unlink of anything: snakemake deletes a failed job's declared output, so
    # the presence of <stage>.json IS the statement "this stage succeeded" and
    # the DAG can never build on top of a failure. A failure->success transition
    # publishes the manifest; a success->failure has the manifest removed for us,
    # and the log is overwritten with the new verdict either way.
    write_if_changed(args.log, text)
    if ok:
        write_if_changed(args.manifest, text)

    for runner, r in manifest["runners"].items():
        tag = {"complete": "OK", "warn": "warn", "below_floor": "<-- BELOW floor"}
        print(f"[completeness]   {runner}: {r['found']}/{r['expect']} "
              f"(floor {r['floor']}) {tag[r['status']]}", file=sys.stderr)
    print(f"[completeness] {args.stage} {unit}: {manifest['status']} "
          f"-> {args.log}" + (f" + {args.manifest}" if ok else ""), file=sys.stderr)
    for f in manifest["failures"]:
        for reason in f["reasons"]:
            print(f"[completeness]   {f['runner']}: {reason}", file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
