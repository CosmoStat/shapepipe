#!/usr/bin/env python3
"""The count-floor completeness table — the single failure policy.

This is the ported ``complete_check`` count table from the v2.0 bash layer
(``run_job_sp_canfar_v2.0.bash`` job dispatch, survey §4). It is the *only*
failure policy in the design: there is no 3-class taxonomy and no error-signature
whitelist. A stage is a real failure iff a mandatory runner produced fewer than
its ``floor`` files; per-CCD attrition (a sparse CCD setools rejects, ~0.2%)
sits between ``floor`` and ``expect`` and is tolerated.

Two consumers share this table:
  * ``sp_rule.py`` — after ``shapepipe_run``, counts products per runner and
    exits nonzero if any mandatory runner is below its floor (``keep-going``
    then isolates that unit's cone; ``retries`` handle transient failures).
  * ``run_report.py`` — disk-scans finished trees against ``expect`` to
    enumerate shortfalls for the human, mid-run or after.

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

# stage -> {runner_subdir: {expect, floor, [warn], [subpath]}}
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
    "exp_psf": {
        "sextractor_runner":   dict(expect=80, floor=2),
        "setools_runner":      dict(expect=80, floor=2, subpath="rand_split"),
        "psfex_runner":        dict(expect=80, floor=2),
        "psfex_interp_runner": dict(expect=40, floor=0, warn=True),
    },

    # --- tile post ---
    "tile_merge_headers": {"merge_headers_runner": dict(expect=1, floor=1)},
    "tile_mask":          {"mask_runner":          dict(expect=1, floor=1)},
    "tile_detect":        {"sextractor_runner":    dict(expect=2, floor=2)},
    "tile_vignets": {
        "psfex_interp_runner":     dict(expect=1, floor=1),
        "vignetmaker_runner_run_1": dict(expect=1, floor=1),
        "vignetmaker_runner_run_2": dict(expect=4, floor=4),
    },
    "tile_ngmix":     {"ngmix_runner":          dict(expect=1, floor=1)},
    "tile_merge_cats": {"merge_sep_cats_runner": dict(expect=1, floor=1)},
    "tile_make_cat":  {"make_cat_runner":       dict(expect=1, floor=1)},
}


def count_products(run_dir, runner, spec):
    """Count files in ``run_dir/<runner>/output[/<subpath>]/`` (live links only)."""
    out = run_dir / runner / "output"
    if "subpath" in spec:
        out = out / spec["subpath"]
    if not out.is_dir():
        return 0
    return sum(1 for p in out.iterdir() if p.exists())  # p.exists() drops dead links


def check_floor(stage, run_dir):
    """Return (ok, details). ok is False iff a mandatory runner is below floor.

    ``details`` is a list of (runner, n_found, floor, expect, warn) tuples.
    """
    table = COMPLETENESS[stage]
    details, ok = [], True
    for runner, spec in table.items():
        n = count_products(run_dir, runner, spec)
        details.append((runner, n, spec["floor"], spec["expect"],
                        spec.get("warn", False)))
        if not spec.get("warn", False) and n < spec["floor"]:
            ok = False
    return ok, details
