#!/usr/bin/env python3
"""Thin per-unit ShapePipe rule wrapper — the single entrypoint every rule calls.

Under the one-allocation local executor + apptainer software-deployment,
Snakemake runs each rule's shell *inside* the container on the allocation's node.
So this script already runs in-container: it calls ``shapepipe_run`` directly
(never apptainer itself), which is why the workflow never hand-invokes apptainer.

What it owns — only what is intrinsically ShapePipe's, nothing Snakemake already
does:

1. **Per-unit work dir + v2.0 isolation furniture.** Each tile/exposure gets its
   own ``$SP_RUN`` (``tiles/<ID>/`` or ``exp/<base>/``). Isolation is by work-dir
   *content*, the proven v2.0 mechanism — NOT ``-e/--exclusive``:
     - a ``cfis`` symlink to the canonical config dir (so ``$SP_CONFIG`` / ``.mask``
       refs resolve);
     - ``star_cat_exp`` / ``star_cat_tiles`` symlinks to the shared cat pool;
     - for a tile: ``tile_numbers.txt`` = this tile ID (dot format, real data —
       what get_images reads);
     - for an exposure: a fabricated
       ``output/run_sp_tile_Fe/find_exposures_runner/output/exp_numbers-000-000.txt``
       = this exposure ID (v2.0 ``init_exp_work_dir`` pseudo-Fe), so the Gie
       config's ``last:find_exposures_runner`` resolver finds exactly one unit.

2. **Per-unit config copy.** Forces ``RUN_DATETIME = False`` (deterministic
   ``output/run_sp_<prefix>/`` the rule declares as its ``directory()`` output),
   sets ``SMP_BATCH_SIZE = {threads}`` (fork width == cpus_per_task, so SLURM
   packing density and true per-job process count are one number), and — only
   when ``isolate`` (tile-scheme stages + exp split; NEVER get_images / exp_mask /
   exp_psf, whose per-CCD/download numbering would make #746 turn tolerated
   per-CCD attrition into a whole-exposure hard failure) — ``NUMBER_LIST = -<ID>``.

3. **Log sync at START and end.** ``update_runs_log_file.py`` regenerates
   ``log_run_sp.txt`` before ``shapepipe_run`` (after the pre-delete, so
   resolution never sees a just-deleted dir) and again after. Skipped for ngmix
   chunks: the log is already complete from the prior tile stages, and the N
   concurrent chunk jobs share the tile work dir — a full log *rewrite* would
   race, whereas ShapePipe's own single-line ``O_APPEND`` of each chunk's run is
   atomic. ngmix chunk configs additionally carry *literal* INPUT_DIRs
   (deterministic under RUN_DATETIME=False), so they need no log read at all.

4. **Count-floor check** (``completeness.py`` — the ported bash ``complete_check``
   table, the single failure policy: no 3-class taxonomy, no whitelist). After
   the run, count products per mandatory runner and exit nonzero below the floor.
   ``keep-going`` isolates that unit's cone; per-CCD attrition between floor and
   ``expect`` is tolerated. Stages absent from the table fall back to a
   zero-output floor.

Stage-specific verbs (grounded in ``runs/p3-batch1/job_p3_batch1.sh``):

- ``--ngmix-chunk K --ngmix-nchunks N`` — read this tile's SExtractor object
  count from its ``run_sp_tile_Sx`` output, split ``[1, N_obj]`` into N disjoint
  closed ID ranges, and generate chunk K's config from
  ``config_tile_Ng_template_batch.ini``.
- ``--n-split-max N`` — generate the merge config from
  ``config_merge_sep_cats_template.ini`` with ``N_SPLIT_MAX = N``.
- ``--final-cat PATH`` — after ``make_cat``, copy the produced ``final_cat`` FITS
  to PATH (the rule's declared, protected science-product output).

Run by hand to reproduce exactly what a rule does (the escape hatch — it
pre-deletes the run dir, just like Snakemake's delete-before-rerun)::

    python sp_rule.py --stage exp_split --config config_exp_Sp.ini \\
        --unit 2243881 --level exp --run-dir $SP_RUN --config-src $SP_CONFIG \\
        --star-cats $STAR_CATS --threads 8
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from completeness import COMPLETENESS, check_floor  # noqa: E402


def unit_work_dir(run_dir: Path, level: str, unit: str) -> Path:
    """Per-unit work dir in the forest: tiles/<ID> or exp/<base>.

    Flat (no IDra/prefix sharding): the shard existed only to keep the bash
    layer's ``ls exp/`` bearable, but Snakemake never globs the forest — the
    run index drives every path — so flattening keeps single-wildcard rules.
    ``run`` level (aggregation stages) work directly in ``run_dir``.
    """
    if level == "tile":
        return run_dir / "tiles" / unit
    if level == "exp":
        base = unit[:-1] if unit[-1].isalpha() else unit
        return run_dir / "exp" / base
    return run_dir


def materialise_unit(work: Path, level: str, unit: str, config_src: Path,
                     star_cats: Path, exp_name: str | None = None) -> None:
    """Create the v2.0 unit-isolation furniture in the work dir."""
    (work / "output").mkdir(parents=True, exist_ok=True)

    cfis = work / "cfis"  # $SP_CONFIG -> canonical config dir (.mask refs resolve)
    if not cfis.exists():
        cfis.symlink_to(config_src)

    for name, sub in (("star_cat_exp", "exp"), ("star_cat_tiles", "tiles")):
        link = work / name
        if not link.exists() and (star_cats / sub).exists():
            link.symlink_to(star_cats / sub)

    if level == "tile":
        (work / "tile_numbers.txt").write_text(unit + "\n")  # dot format, real data
    elif level == "exp":
        fe = work / "output" / "run_sp_tile_Fe" / "find_exposures_runner" / "output"
        fe.mkdir(parents=True, exist_ok=True)
        exp_numbers = fe / "exp_numbers-000-000.txt"
        # Written UNCONDITIONALLY (cheap, deterministic): an exists-guard once
        # pinned a stale pre-fix file with the bare base id. Content is the
        # ORIGINAL exposure name (with its 'p'-style suffix, from the index) —
        # the on-disk store is <name>.fits.fz; the bare base id matches nothing.
        # v2.0 copies the name verbatim from the tile's Fe output; the index
        # carries it for us.
        exp_numbers.write_text((exp_name or unit) + "\n")


def normalize_config_text(text: str, unit: str, isolate: bool, threads: int) -> str:
    """Force RUN_DATETIME=False, SMP_BATCH_SIZE=threads, and (if isolating)
    NUMBER_LIST=-<ID dashed>.

    ``NUMBER_LIST`` uses ShapePipe's image-number convention (dots -> dashes),
    the ``set_config_number_list`` mechanism that replaced the retired
    ``-e/--exclusive`` flag (#746).
    """
    lines = text.splitlines()
    out, in_file, saw_datetime, saw_numberlist, saw_smp = [], False, False, False, False
    number = "-" + unit.replace(".", "-")
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_file = stripped.upper() == "[FILE]"
        if stripped.upper().startswith("RUN_DATETIME"):
            out.append("RUN_DATETIME = False")
            saw_datetime = True
            continue
        if stripped.upper().startswith("SMP_BATCH_SIZE"):
            out.append(f"SMP_BATCH_SIZE = {threads}")
            saw_smp = True
            continue
        if isolate and stripped.upper().startswith("NUMBER_LIST"):
            out.append(f"NUMBER_LIST = {number}")
            saw_numberlist = True
            continue
        out.append(line)
        if in_file and stripped.upper() == "[FILE]" and isolate and not saw_numberlist:
            out.append(f"NUMBER_LIST = {number}")
            saw_numberlist = True
    if not saw_datetime:
        # ShapePipe reads RUN_DATETIME from [DEFAULT]; inject after its header.
        patched, done = [], False
        for line in out:
            patched.append(line)
            if not done and line.strip().upper() == "[DEFAULT]":
                patched.append("RUN_DATETIME = False")
                done = True
        out = patched if done else ["[DEFAULT]", "RUN_DATETIME = False"] + out
    if not saw_smp:
        # SMP_BATCH_SIZE lives in [JOB]; inject after its header.
        patched, done = [], False
        for line in out:
            patched.append(line)
            if not done and line.strip().upper() == "[JOB]":
                patched.append(f"SMP_BATCH_SIZE = {threads}")
                done = True
        out = patched if done else out + ["[JOB]", f"SMP_BATCH_SIZE = {threads}"]
    return "\n".join(out) + "\n"


def run_prefix_from_text(text: str, config_label: str) -> str:
    """The RUN_NAME (== run_sp_<prefix>) the module will write under."""
    for line in text.splitlines():
        s = line.strip()
        if s.upper().startswith("RUN_NAME"):
            return s.split("=", 1)[1].strip()
    raise SystemExit(f"no RUN_NAME in {config_label}")


# --- ngmix chunking --------------------------------------------------------

# Literal, deterministic INPUT_DIRs for the ngmix chunk (RUN_DATETIME=False), so
# concurrent chunks need no log read at all (finding 18).
NGMIX_LITERAL_INPUT = (
    "./output/run_sp_tile_Sx/sextractor_runner/output, "
    "./output/run_sp_tile_PiViVi/psfex_interp_runner/output, "
    "./output/run_sp_tile_PiViVi/vignetmaker_runner_run_2/output, "
    "./output/run_sp_tile_Mh_exp/merge_headers_runner/output")


def ngmix_id_ranges(n_obj: int, n_chunks: int) -> list[tuple[int, int]]:
    """Split ``[1, n_obj]`` into ``n_chunks`` disjoint *closed* ID ranges.

    SExtractor's ``NUMBER`` column (what ngmix reads as ``obj_id``) runs 1..N
    contiguous, so covering ``[1, N]`` processes every object exactly once. The
    first n_chunks-1 chunks get ``n_obj // n_chunks`` each; the last gets the
    remainder — a *closed* upper bound at ``n_obj`` (never ``ID_OBJ_MAX = -1``:
    ngmix treats ``id_obj_max <= 0`` as unbounded, which would double-count).

    Deviation from the nibi monolith (deliberate): the monolith sizes every
    tile's chunks from the *average* object count across all tiles and leaves the
    last chunk open (``-1``); here each tile is an independent job, so we size
    from *this* tile's own count and close the last chunk. The merged catalogue
    is identical either way (merge_sep_cats recombines all chunks, ngmix fits
    each object independently) — this just removes the open-ended overflow.
    """
    base, rem = divmod(n_obj, n_chunks)
    ranges, lo = [], 1
    for k in range(1, n_chunks + 1):
        size = base + (rem if k == n_chunks else 0)
        hi = lo + size - 1
        ranges.append((lo, hi))
        lo = hi + 1
    return ranges


def tile_object_count(work: Path) -> int:
    """Object count for this tile from its run_sp_tile_Sx SExtractor catalogue.

    Matches ``get_number_objects.py``: the ``NAXIS2`` of the catalogue's last
    HDU. Per-unit isolation means exactly one ``sexcat*.fits`` lives here.
    """
    from astropy.io import fits

    sexcats = sorted((work / "output" / "run_sp_tile_Sx").glob(
        "sextractor_runner/output/sexcat*.fits"))
    if not sexcats:
        raise SystemExit(
            f"[sp_rule] FATAL ngmix: no sexcat under {work}/output/run_sp_tile_Sx")
    with fits.open(sexcats[0]) as hdul:
        return int(hdul[-1].header["NAXIS2"])


def apply_ngmix_subs(text: str, chunk: int, lo: int, hi: int) -> str:
    """Fill the ngmix batch template for one chunk (the monolith's sed recipe)
    plus literal INPUT_DIRs so the chunk needs no shared-log read."""
    text = text.replace("NgXu", f"Ng{chunk}u").replace("X_interp", "psfex_interp")
    out = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("ID_OBJ_MIN"):
            out.append(f"ID_OBJ_MIN = {lo}")
        elif s.startswith("ID_OBJ_MAX"):
            out.append(f"ID_OBJ_MAX = {hi}")
        elif s.startswith("INPUT_DIR") and "sextractor_runner" in s:
            out.append(f"INPUT_DIR = {NGMIX_LITERAL_INPUT}")
        else:
            out.append(line)
    return "\n".join(out)


def apply_merge_subs(text: str, n_split_max: int) -> str:
    """Fill N_SPLIT_MAX in the merge template (the monolith's sed recipe)."""
    out = []
    for line in text.splitlines():
        if line.strip().startswith("N_SPLIT_MAX"):
            out.append(f"N_SPLIT_MAX = {n_split_max}")
        else:
            out.append(line)
    return "\n".join(out)


# --- in-job count floor ----------------------------------------------------

def count_floor(run_out: Path, stage: str) -> None:
    """Fail loudly if any mandatory runner is below its ``floor`` (completeness.py).

    Stages absent from the table skip this floor (only the zero-output floor
    applies to them).
    """
    if stage not in COMPLETENESS:
        return
    ok, details = check_floor(stage, run_out)
    for runner, n, floor, expect, warn in details:
        tag = "warn" if warn else ("OK" if n >= floor else "  <-- BELOW floor")
        print(f"[sp_rule]   {runner}: {n}/{expect} (floor {floor}) {tag}",
              file=sys.stderr)
    if not ok:
        print(f"[sp_rule] FATAL {stage}: count floor breached in {run_out}",
              file=sys.stderr)
        sys.exit(1)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", required=True)
    p.add_argument("--config", required=True, help="config filename in --config-src")
    p.add_argument("--unit", required=True)
    p.add_argument("--level", required=True, choices=["tile", "exp", "run"])
    p.add_argument("--run-dir", required=True, type=Path, help="$SP_RUN root")
    p.add_argument("--config-src", required=True, type=Path,
                   help="dir holding the canonical configs")
    p.add_argument("--star-cats", type=Path, default=None,
                   help="shared star-cat pool (tiles/ + exp/ subdirs)")
    p.add_argument("--threads", type=int, default=1,
                   help="SMP_BATCH_SIZE == cpus_per_task (fork width)")
    p.add_argument("--exp-forest", type=Path, default=None,
                   help="$SP_EXP: per-tile exposure forest (tile post-stages)")
    p.add_argument("--exp-name", default=None,
                   help="original exposure name incl. suffix (e.g. 2605805p); "
                        "written into the fabricated exp_numbers list so "
                        "get_images finds <name>.fits.fz in the store")
    p.add_argument("--no-isolate", action="store_true",
                   help="skip NUMBER_LIST (get_images / exp_mask / exp_psf)")
    p.add_argument("--ngmix-chunk", type=int, default=None,
                   help="1-indexed ngmix chunk K (with --ngmix-nchunks)")
    p.add_argument("--ngmix-nchunks", type=int, default=None,
                   help="total ngmix chunks N (with --ngmix-chunk)")
    p.add_argument("--n-split-max", type=int, default=None,
                   help="N_SPLIT_MAX for the merge_sep_cats template")
    p.add_argument("--final-cat", type=Path, default=None,
                   help="copy the produced final_cat FITS here (make_cat)")
    args = p.parse_args()

    ngmix_mode = args.ngmix_chunk is not None
    if ngmix_mode and args.ngmix_nchunks is None:
        raise SystemExit("--ngmix-chunk requires --ngmix-nchunks")

    work = unit_work_dir(args.run_dir, args.level, args.unit)
    materialise_unit(work, args.level, args.unit, args.config_src,
                     args.star_cats or args.run_dir, exp_name=args.exp_name)

    # Build the config copy. ngmix chunks and merge fill template placeholders
    # first; every copy is then normalized (RUN_DATETIME=False, SMP_BATCH_SIZE,
    # [+ NUMBER_LIST]).
    src_text = (args.config_src / args.config).read_text()
    isolate = not args.no_isolate
    if ngmix_mode:
        n_obj = tile_object_count(work)
        lo, hi = ngmix_id_ranges(n_obj, args.ngmix_nchunks)[args.ngmix_chunk - 1]
        print(f"[sp_rule] ngmix chunk {args.ngmix_chunk}/{args.ngmix_nchunks}: "
              f"n_obj={n_obj} -> ID_OBJ [{lo}, {hi}]", file=sys.stderr)
        src_text = apply_ngmix_subs(src_text, args.ngmix_chunk, lo, hi)
        cfg_name = f"config_tile_Ng{args.ngmix_chunk}u.ini"  # distinct per chunk
    elif args.n_split_max is not None:
        src_text = apply_merge_subs(src_text, args.n_split_max)
        cfg_name = args.config
    else:
        cfg_name = args.config

    cfg_mod_dir = work / "cfis_mod"
    cfg_mod = cfg_mod_dir / cfg_name
    cfg_text = normalize_config_text(src_text, args.unit, isolate, args.threads)
    cfg_mod_dir.mkdir(parents=True, exist_ok=True)
    cfg_mod.write_text(cfg_text)

    prefix = run_prefix_from_text(cfg_text, str(cfg_mod))
    run_out = work / "output" / prefix

    # Pre-delete the deterministic run dir (delete-before-rerun for the hand-run
    # escape hatch; ShapePipe's FileHandler.mkdir raises on an existing run dir).
    if run_out.exists():
        shutil.rmtree(run_out)

    env = dict(os.environ)
    env["SP_RUN"] = str(work)
    env["SP_CONFIG"] = str(work / "cfis")  # canonical dir (.mask refs), not cfis_mod
    if args.exp_forest:
        env["SP_EXP"] = str(args.exp_forest)
    # Thread caps also come via apptainer-args; keep here for the hand-run path.
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    env.setdefault("MALLOC_ARENA_MAX", "2")
    env.setdefault("MALLOC_TRIM_THRESHOLD_", "0")

    # Log sync before the run (after pre-delete). Skipped for ngmix chunks: N
    # concurrent chunks would race a full rewrite (literal INPUT_DIRs make it moot).
    if not ngmix_mode:
        subprocess.call(["update_runs_log_file.py"], cwd=str(work), env=env)

    print(f"[sp_rule] stage={args.stage} unit={args.unit} work={work}",
          file=sys.stderr)
    rc = subprocess.call(["shapepipe_run", "-c", str(cfg_mod)],
                         cwd=str(work), env=env)
    if rc != 0:
        print(f"[sp_rule] WARN shapepipe_run exit {rc} — checking count floor",
              file=sys.stderr)

    if not ngmix_mode:
        subprocess.call(["update_runs_log_file.py"], cwd=str(work), env=env)

    # Zero-output floor for stages without a count table; count floor otherwise.
    if args.stage in COMPLETENESS:
        count_floor(run_out, args.stage)
    else:
        produced = list(run_out.glob("**/output/*")) if run_out.exists() else []
        if not produced:
            print(f"[sp_rule] FATAL {args.stage} {args.unit}: zero output in {run_out}",
                  file=sys.stderr)
            sys.exit(1)

    # make_cat: publish the science product to its declared, protected path.
    if args.final_cat is not None:
        finals = sorted(run_out.glob("**/output/final_cat*.fits"))
        if not finals:
            print(f"[sp_rule] FATAL {args.stage} {args.unit}: no final_cat under {run_out}",
                  file=sys.stderr)
            sys.exit(1)
        args.final_cat.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(finals[0], args.final_cat)
        print(f"[sp_rule] final_cat -> {args.final_cat}", file=sys.stderr)

    print(f"[sp_rule] OK {args.stage} {args.unit}", file=sys.stderr)


if __name__ == "__main__":
    main()
