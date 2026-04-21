#!/usr/bin/env python
"""Test and compare tile object-detection modes (sx vs uc).

Two stages, each independently switchable:

  run     -- execute ShapePipe jobs 128 and 256 (both modes) for each tile
  compare -- read sexcat FITS-LDAC catalogues, produce per-quantity plots for
             each flag filter, and compile a LaTeX/PDF summary document

By default both stages run.  Use -R/--no_run or -C/--no_compare to skip one.
All comparison options (quantities, flag columns, ranges) are read from a
config file (default: test_tile_det.cfg in the current directory).

Usage examples
--------------
    # Run jobs AND compare for one tile:
    test_tile_det.py -e 301.279

    # Multiple tiles, compare only, custom config:
    test_tile_det.py -e "301.279 301.001" -R -c my.cfg

    # Dry run to preview commands:
    test_tile_det.py -e 301.279 -n

    # Force re-run of jobs even if output exists:
    test_tile_det.py -e 301.279 -f
"""

import configparser
import glob
import os
import subprocess
import sys

from cs_util import args as cs_args
from cs_util import logging as cs_logging

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


# ---------------------------------------------------------------------------
# Default parameters (cs_util.args style)
# ---------------------------------------------------------------------------

def params_default():
    """Return default parameter dicts for cs_util.args.parse_options."""

    _params = {
        "exclusive":  None,
        "base_dir":   os.path.expanduser("~/v2.0"),
        "psf":        "psfex",
        "config":     "test_tile_det.cfg",
        "output_dir": "test_tile_det_output",
        "force":      False,
        "dry_run":    False,
        "no_run":     False,
        "no_compare": False,
    }

    _short_options = {
        "exclusive":  "-e",
        "base_dir":   "-d",
        "psf":        "-p",
        "config":     "-c",
        "output_dir": "-o",
        "force":      "-f",
        "dry_run":    "-n",
        "no_run":     "-R",
        "no_compare": "-C",
    }

    _types = {
        "force":      "bool",
        "dry_run":    "bool",
        "no_run":     "bool",
        "no_compare": "bool",
    }

    _help_strings = {
        "exclusive":  "tile ID(s), space-separated, e.g. '301.279 301.001'; default={}",
        "base_dir":   "v2.0 base directory (contains tiles/); default={}",
        "psf":        "PSF model, psfex|mccd; default={}",
        "config":     "comparison config file; default={}",
        "output_dir": "output directory for plots and summary; default={}",
        "force":      "force re-run even if output already exists; default={}",
        "dry_run":    "print commands without executing; default={}",
        "no_run":     "skip the ShapePipe job stage; default={}",
        "no_compare": "skip the catalogue comparison/plotting stage; default={}",
    }

    return _params, _short_options, _types, _help_strings


def check_params(params):
    """Validate parameters. Return False on error."""
    if params["exclusive"] is None:
        print("Error: tile ID(s) required, use -e")
        return False
    if params["psf"] not in ("psfex", "mccd"):
        print(f"Error: psf must be 'psfex' or 'mccd', got '{params['psf']}'")
        return False
    return True


# ---------------------------------------------------------------------------
# Stage 1 — run ShapePipe jobs
# ---------------------------------------------------------------------------

def run_job(cmd, dry_run=False):
    """Print and optionally execute *cmd*. Returns the return code."""
    print(f"  + {cmd}")
    if dry_run:
        return 0
    result = subprocess.run(cmd, shell=True)
    return result.returncode


def run_tile_jobs(tile_id, base_dir, psf, force, dry_run, cfg):
    """Run prerequisite jobs then the detection job for each mode.

    Job numbers and modes are read from *cfg* [jobs]:
      prereq_jobs  -- bit-coded job sum to run first (ensures all inputs exist)
      det_job      -- the detection job number (e.g. 256)
      det_modes    -- space-separated tile_det modes (e.g. 'sx uc')
    """
    script = "run_job_sp_canfar_v2.0.bash"
    force_flag = "--force" if force else ""
    common = f"{script} -e {tile_id} -d {base_dir} -p {psf} {force_flag}".strip()

    prereq_jobs = cfg.getint("jobs", "prereq_jobs")
    det_job     = cfg.getint("jobs", "det_job")
    det_modes   = cfg.get("jobs", "det_modes").split()

    print(f"\n--- Tile {tile_id}: prerequisite jobs (bit-coded {prereq_jobs}) ---")
    rc = run_job(f"{common} -j {prereq_jobs}", dry_run)
    if rc != 0:
        print(
            f"WARNING: prerequisite jobs ({prereq_jobs}) returned {rc} "
            f"for tile {tile_id}",
            file=sys.stderr,
        )

    for mode in det_modes:
        print(f"\n--- Tile {tile_id}: job {det_job} --tile_det {mode} ---")
        rc = run_job(f"{common} -j {det_job} --tile_det {mode}", dry_run)
        if rc != 0:
            print(
                f"WARNING: job {det_job} --tile_det {mode} returned {rc} "
                f"for tile {tile_id}",
                file=sys.stderr,
            )


# ---------------------------------------------------------------------------
# Stage 2 — catalogue I/O helpers
# ---------------------------------------------------------------------------

def find_sexcat(tile_dir, mode):
    """Return the most recent sexcat FITS-LDAC file for *mode* ('sx' or 'uc').

    Raises FileNotFoundError if no match is found.
    """
    if mode == "sx":
        pattern = os.path.join(
            tile_dir, "output", "run_sp_tile_Sx*",
            "sextractor_runner", "output", "sexcat-*.fits",
        )
    else:
        pattern = os.path.join(
            tile_dir, "output", "run_sp_tile_Uc*",
            "read_ext_cat_runner", "output", "sexcat-*.fits",
        )
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No sexcat for mode '{mode}' in {tile_dir}\n"
            f"  (searched: {pattern})"
        )
    return matches[-1]


def read_ldac_objects(path):
    """Return the LDAC_OBJECTS table from *path* as an astropy FITS rec."""
    with fits.open(path) as hdul:
        return hdul["LDAC_OBJECTS"].data


def apply_flag_filter(data, flag_col, condition):
    """Return *data* masked according to *condition* on *flag_col*.

    Parameters
    ----------
    data : FITS rec
    flag_col : str
        Column name carrying integer FLAGS.
    condition : str
        ``'all'``, ``'eq0'`` (FLAGS==0), or ``'le2'`` (FLAGS<=2).
    """
    if condition == "all":
        return data
    if flag_col not in data.names:
        print(
            f"  WARNING: flag column '{flag_col}' not found; "
            "skipping flag filter",
            file=sys.stderr,
        )
        return data
    flags = np.asarray(data[flag_col], dtype=int)
    if condition == "eq0":
        mask = flags == 0
    elif condition == "le2":
        mask = flags <= 2
    else:
        return data
    return data[mask]


# ---------------------------------------------------------------------------
# Stage 2 — comparison plots
# ---------------------------------------------------------------------------

# Human-readable filter labels (plain text version for plot titles)
FILTER_LABELS = {
    "all": "All objects",
    "eq0": "FLAGS = 0",
    "le2": "FLAGS ≤ 2",
}

# LaTeX versions for the PDF summary
FILTER_LABELS_TEX = {
    "all": "All objects",
    "eq0": r"FLAGS $= 0$",
    "le2": r"FLAGS $\leq 2$",
}


def compare_quantity(
    tile_dirs, qty, cfg, output_dir, flag_col_sx, flag_col_uc, verbose=False
):
    """Produce one PNG for a single quantity with all filters side by side.

    Layout: n_tiles rows × n_filters columns.

    Parameters
    ----------
    tile_dirs : list of str
    qty : str
        Section name in *cfg* (e.g. ``'magnitude'``).
    cfg : configparser.ConfigParser
    output_dir : str
    flag_col_sx, flag_col_uc : str
    verbose : bool

    Returns
    -------
    filename of the generated PNG.
    """
    label     = cfg.get(qty, "label")
    col_sx    = cfg.get(qty, "col_sx")
    col_uc    = cfg.get(qty, "col_uc")
    range_min = cfg.getfloat(qty, "range_min")
    range_max = cfg.getfloat(qty, "range_max")
    bin_width = cfg.getfloat(qty, "bin_width")

    filters   = cfg.get("comparison", "filters").split()
    bins      = np.arange(range_min, range_max + bin_width, bin_width)
    n_tiles   = len(tile_dirs)
    n_filters = len(filters)

    fig, axes = plt.subplots(
        n_tiles, n_filters,
        figsize=(3.5 * n_filters, 3.2 * n_tiles),
        squeeze=False,
    )

    for row, tile_dir in enumerate(tile_dirs):
        tile_id = os.path.basename(tile_dir)

        for col, filt in enumerate(filters):
            ax = axes[row][col]

            for mode, mode_col, flag_col, mode_label, color in [
                ("sx", col_sx, flag_col_sx, "SExtractor", "steelblue"),
                ("uc", col_uc, flag_col_uc, "UNIONS",     "tomato"),
            ]:
                try:
                    path = find_sexcat(tile_dir, mode)
                    data = read_ldac_objects(path)
                    data = apply_flag_filter(data, flag_col, filt)

                    if mode_col not in data.names:
                        print(
                            f"  WARNING: column '{mode_col}' not in {path}",
                            file=sys.stderr,
                        )
                        continue

                    values = np.asarray(data[mode_col], dtype=float)
                    values = values[
                        (values > range_min) & (values < range_max)
                    ]
                    ax.hist(
                        values, bins=bins,
                        histtype="step", color=color, linewidth=1.5,
                        label=f"{mode_label}  N={len(values):,}",
                    )
                    if verbose:
                        print(
                            f"  {tile_id} {mode} {filt}: "
                            f"{len(values):,} objects  ({path})"
                        )

                except (FileNotFoundError, KeyError) as exc:
                    print(f"  WARNING: {exc}", file=sys.stderr)

            filt_label = FILTER_LABELS.get(filt, filt)
            title = f"Tile {tile_id} — {filt_label}" if n_tiles > 1 else filt_label
            ax.set_xlabel(label, fontsize=8)
            ax.set_ylabel(f"N / {bin_width}", fontsize=8)
            ax.set_title(title, fontsize=8)
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=7)

    fig.suptitle(f"{label} — sx vs uc", fontsize=10)
    plt.tight_layout()

    fname = f"{qty}.png"
    fpath = os.path.join(output_dir, fname)
    plt.savefig(fpath, dpi=150)
    plt.close()
    print(f"  Saved {fpath}")
    return fname


def compare_tiles(tile_dirs, cfg, output_dir, verbose=False):
    """Run all comparisons defined in *cfg* and collect generated plot paths.

    Returns
    -------
    list of (qty, label, filename) tuples.
    """
    os.makedirs(output_dir, exist_ok=True)

    quantities  = cfg.get("comparison", "quantities").split()
    flag_col_sx = cfg.get("flags", "col_sx")
    flag_col_uc = cfg.get("flags", "col_uc")

    all_plots = []
    for qty in quantities:
        label = cfg.get(qty, "label")
        print(f"\n  Quantity: {label}")
        fname = compare_quantity(
            tile_dirs, qty, cfg, output_dir,
            flag_col_sx, flag_col_uc, verbose=verbose,
        )
        all_plots.append((qty, label, fname))

    return all_plots


# ---------------------------------------------------------------------------
# Stage 2 — LaTeX / PDF summary
# ---------------------------------------------------------------------------

def _tex_escape(s):
    """Escape special LaTeX characters in plain-text strings."""
    return s.replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")


def generate_summary(tile_ids, all_plots, output_dir):
    """Write summary.tex and compile to summary.pdf via pdflatex.

    Parameters
    ----------
    tile_ids : list of str
    all_plots : list of (qty, label, filename) tuples
    output_dir : str
    """
    tex_path = os.path.join(output_dir, "summary.tex")

    tile_list = ", ".join(
        f"\\texttt{{{_tex_escape(tid)}}}" for tid in tile_ids
    )

    lines = [
        r"\documentclass[a4paper,11pt]{article}",
        r"\usepackage{graphicx}",
        r"\usepackage[margin=2cm]{geometry}",
        r"\usepackage{hyperref}",
        r"\begin{document}",
        r"\title{Tile detection comparison: \texttt{sx} vs \texttt{uc}}",
        r"\author{}",
        r"\date{\today}",
        r"\maketitle",
        "",
        r"\noindent Tiles: " + tile_list + r"\\",
        r"\noindent Modes: \texttt{sx} (SExtractor detection), "
        r"\texttt{uc} (UNIONS catalogue).\\",
        "",
    ]

    for qty, label, fname in all_plots:
        lines += [
            f"\\section{{{_tex_escape(label)}}}",
            "",
            r"\begin{figure}[h!]",
            r"  \centering",
            f"  \\includegraphics[width=\\linewidth]{{{fname}}}",
            r"\end{figure}",
            r"\clearpage",
            "",
        ]

    lines.append(r"\end{document}")

    with open(tex_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"  Written {tex_path}")

    # Compile (run twice so hyperref bookmarks are correct)
    try:
        for _ in range(2):
            result = subprocess.run(
                [
                    "pdflatex",
                    "-interaction=nonstopmode",
                    "-output-directory", output_dir,
                    tex_path,
                ],
                capture_output=True,
                text=True,
                cwd=output_dir,
            )
        if result.returncode == 0:
            pdf_path = os.path.join(output_dir, "summary.pdf")
            print(f"  Compiled {pdf_path}")
        else:
            print(
                f"  WARNING: pdflatex returned {result.returncode}; "
                "check summary.log",
                file=sys.stderr,
            )
    except FileNotFoundError:
        print(
            "  WARNING: pdflatex not found — LaTeX written but not compiled",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    _params, _short_options, _types, _help_strings = params_default()

    options = cs_args.parse_options(
        _params, _short_options, _types, _help_strings,
    )
    for key in options:
        _params[key] = options[key]

    cs_logging.log_command(argv)

    if not check_params(_params):
        return 1

    tile_ids  = _params["exclusive"].split()
    base_dir  = os.path.abspath(_params["base_dir"])
    tiles_dir = os.path.join(base_dir, "tiles")
    tile_dirs = [
        os.path.join(tiles_dir, tid.split(".")[0], tid)
        for tid in tile_ids
    ]
    output_dir = os.path.abspath(_params["output_dir"])

    # Load config (needed by both stages)
    cfg_path = _params["config"]
    if not os.path.exists(cfg_path):
        print(f"Error: config file '{cfg_path}' not found", file=sys.stderr)
        return 1
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

    # --- Stage 1: run jobs ---
    if not _params["no_run"]:
        print("=" * 60)
        print("Stage 1: running ShapePipe jobs")
        print("=" * 60)
        for tid in tile_ids:
            run_tile_jobs(
                tile_id=tid,
                base_dir=base_dir,
                psf=_params["psf"],
                force=_params["force"],
                dry_run=_params["dry_run"],
                cfg=cfg,
            )

    # --- Stage 2: compare ---
    if not _params["no_compare"]:
        print("=" * 60)
        print("Stage 2: comparing catalogues")
        print("=" * 60)

        all_plots = compare_tiles(
            tile_dirs=tile_dirs,
            cfg=cfg,
            output_dir=output_dir,
            verbose=_params.get("verbose", False),
        )

        print("\n  Generating summary document …")
        generate_summary(tile_ids, all_plots, output_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
