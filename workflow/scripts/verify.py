#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["pyyaml", "h5py", "numpy", "astropy", "matplotlib"]
# ///
"""Run one ShapePipe stage through workflow/bin/sp on a small fixture and inspect the output.

Usage:
    uv run workflow/scripts/verify.py --stage psf
    uv run workflow/scripts/verify.py --stage all --baseline /path/to/baseline

Environment defaults (overridden by CLI flags):
    SP_VERIFY_INPUTS     Staged inputs root (tiles/, exposures/, catalogues/)
    SP_VERIFY_MASKS      External sky masks directory
    SP_VERIFY_BASELINE   Known-good run directory to compare against
    SP_CONTAINER         ShapePipe container image (.sif)
    SP_SNAKEMAKE_ENV     Snakemake virtualenv (passed through to workflow/bin/sp)
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SIGMA_TO_FWHM = 2.3548
DEFAULT_TILES = ["592.232"]
STAGES = ("prepare", "psf", "detect", "vignets", "all")

# Fir defaults — used only when the corresponding env var and CLI flag are unset
# and the path exists on this machine.
FIR_DEFAULTS = {
    "inputs": "/project/rrg-kyi/cdaley/unions/analyses/shapepipe-patch-testbed/inputs",
    "masks": "/project/def-waerbeke/cdaley/fir-assets/masks/dr6-2026-09-29",
    "container": "/project/def-waerbeke/cdaley/fir-assets/shapepipe-develop-240b37e4.sif",
    "baseline": "/project/rrg-kyi/cdaley/unions/analyses/shapepipe-verify-skill/baseline_develop",
    "snakemake_env": str(Path.home() / "sp/snakemake-env"),
}


def default_path(env_key: str, fir_key: str) -> str | None:
    val = os.environ.get(env_key)
    if val:
        return val
    cand = FIR_DEFAULTS.get(fir_key)
    return cand if cand and Path(cand).exists() else None


def git_info(repo: Path) -> dict:
    def run(*args: str) -> str:
        p = subprocess.run(
            ["git", "-C", str(repo), *args],
            capture_output=True, text=True, check=False,
        )
        return p.stdout.strip() if p.returncode == 0 else ""

    return {
        "repo": str(repo),
        "commit": run("rev-parse", "HEAD"),
        "branch": run("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(run("status", "--porcelain")),
    }


def deep_merge(base: dict, extra: dict) -> dict:
    out = dict(base)
    for k, v in extra.items():
        out[k] = deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def build_config(
    out_dir: Path,
    work_dir: Path,
    inputs_dir: Path,
    masks_dir: Path,
    container: str,
    ngmix_chunks: int,
    overrides_path: Path | None,
    profile: str,
) -> Path:
    cfg = {
        "run": "verify",
        "machine": profile,
        "input_type": "data",
        "retrieve": "symlink",
        "container": container,
        "tile_list": str(out_dir / "tiles.txt"),
        "inputs": {
            "tiles": str(inputs_dir / "tiles"),
            "exposures": str(inputs_dir / "exposures"),
            "catalogues": str(inputs_dir / "catalogues"),
            "masks": str(masks_dir),
        },
        "outputs": {
            "run_dir": str(work_dir),
            "products_dir": str(out_dir / "products"),
            "index_db": str(out_dir / "products/index/run_index.sqlite"),
        },
        "tile_store_root": str(work_dir / "tile-store"),
        "ngmix_chunks": ngmix_chunks,
        "persist_exp": [
            "star_selection", "star_train", "star_test",
            "star_stats", "psf_model", "psfex_cat",
        ],
        "clean": False,
        "clean_tiles": False,
    }
    if overrides_path:
        cfg = deep_merge(cfg, yaml.safe_load(overrides_path.read_text()) or {})
    cfg_path = out_dir / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return cfg_path


def run_launcher(
    stage: str,
    tiles: list[str],
    repo: Path,
    out_dir: Path,
    work_dir: Path,
    inputs_dir: Path,
    masks_dir: Path,
    container: str,
    cfg_path: Path,
    cores: int,
    profile: str,
) -> dict:
    sp_bin = repo / "workflow/bin/sp"
    if not sp_bin.exists():
        sys.exit(f"verify: launcher not found at {sp_bin}")

    state_dir = work_dir / "state"
    for sub in ("tmp", "cache", "state", "tile-store", "exp-store"):
        (work_dir / sub).mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)
    (out_dir / "products").mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["TMPDIR"] = str(work_dir / "tmp")
    env["APPTAINER_TMPDIR"] = str(work_dir / "tmp")
    env["APPTAINER_CACHEDIR"] = str(work_dir / "cache")
    env["SP_STATE_DIR"] = str(state_dir)
    env["SP_MISSING_THRESHOLD"] = "0"
    env["SP_PROFILE"] = profile
    env["SP_CONTAINER"] = container
    if "SP_SNAKEMAKE_ENV" not in env and Path(FIR_DEFAULTS["snakemake_env"]).exists():
        env["SP_SNAKEMAKE_ENV"] = FIR_DEFAULTS["snakemake_env"]

    # sp run snapshots code into $SP_STATE_DIR/code; targeted stages run from the live repo.
    pypath = str(state_dir / "code/src") if stage == "all" else str(repo / "src")
    bind_paths = sorted({
        str(repo.resolve()),
        str(inputs_dir.resolve()),
        str(masks_dir.resolve()),
        str(out_dir.resolve()),
        str(work_dir.resolve()),
    })
    binds = [
        "--cleanenv",
        "--env", "OMP_NUM_THREADS=1",
        "--env", "MALLOC_ARENA_MAX=2",
        "--env", f"PYTHONPATH={pypath}",
        "--env", f"SP_EXP_STORE={work_dir / 'exp-store'}",
        "--home", str(Path.home()),
    ]
    for bp in bind_paths:
        binds.extend(["--bind", bp])
    binds.extend(["--bind", f"{work_dir / 'tile-store'}:/local/scratch"])
    if env.get("SP_STAR_COLOURS"):
        sc_path = Path(env["SP_STAR_COLOURS"]).resolve()
        binds.extend(["--env", f"SP_STAR_COLOURS={sc_path}", "--bind", str(sc_path.parent)])
    apptainer_args = " ".join(binds)

    psf_threads = min(8, max(1, cores // max(1, min(3, len(tiles) * 3))))
    common = [
        "--executor", "local",
        "--cores", str(cores),
        "--apptainer-args", apptainer_args,
        "--set-threads", f"exp_split=4", f"exp_psf={psf_threads}", "tile_merge_headers=1",
        "--set-resources", "exp_psf:mem_mb=10000*attempt", "tile_ngmix:mem_mb=4000*attempt",
    ]

    log_path = out_dir / "logs/pipeline.log"
    t0 = time.time()
    with open(log_path, "w") as logf:
        def invoke(cmd: list[str], extra_env: dict | None = None) -> int:
            call_env = env if extra_env is None else {**env, **extra_env}
            logf.write(f"$ {' '.join(cmd)}\n")
            logf.flush()
            p = subprocess.run(cmd, env=call_env, stdout=logf, stderr=subprocess.STDOUT)
            return p.returncode

        if stage == "all":
            rc = invoke([str(sp_bin), "run", "-c", str(cfg_path), *common])
        else:
            rc = invoke(
                [str(sp_bin), "-c", str(cfg_path), "prepare_all_tiles", *common],
                extra_env={"SP_PHASE": "prepare"},
            )
            if rc == 0 and stage != "prepare":
                targets = stage_targets(stage, tiles, work_dir, out_dir / "products")
                rc = invoke(
                    [str(sp_bin), "-c", str(cfg_path), *targets, *common],
                    extra_env={"SP_PHASE": "compute"},
                )

    elapsed = time.time() - t0
    if rc != 0:
        print(f"[verify] FAILED (exit {rc}) after {elapsed:.1f}s — see {log_path}", file=sys.stderr)
        tail = log_path.read_text(errors="replace").splitlines()[-40:]
        print("\n".join(tail), file=sys.stderr)
        sys.exit(rc)

    return {
        "stage": stage,
        "tiles": tiles,
        "cores": cores,
        "seconds": round(elapsed, 1),
        "core_hours": round(elapsed * cores / 3600.0, 3),
        "log": str(log_path),
    }


def stage_targets(stage: str, tiles: list[str], work_dir: Path, prod_dir: Path) -> list[str]:
    if stage == "psf":
        return [str(prod_dir / "full_starcat_verify.hdf5")]
    if stage == "detect":
        return [str(work_dir / f"tiles/{t[:2]}/{t}/manifests/tile_detect.json") for t in tiles]
    if stage == "vignets":
        return [
            str(prod_dir / "full_starcat_verify.hdf5"),
            *[str(work_dir / f"tiles/{t[:2]}/{t}/manifests/tile_vignets.json") for t in tiles],
        ]
    raise ValueError(f"Unknown targeted stage: {stage}")


# --- loading ----------------------------------------------------------------

def find_starcat(run_dir: Path) -> Path | None:
    hits = sorted((run_dir / "products").glob("full_starcat_*.hdf5"))
    return hits[0] if hits else None


def find_final_cat(run_dir: Path, tile: str) -> Path | None:
    hits = sorted((run_dir / "products/tiles").glob(f"*/{tile}/final_cat-{tile}.hdf5"))
    return hits[0] if hits else None


def find_sexcat(run_dir: Path, tile: str, extra_work: Path | None = None) -> Path | None:
    dashed = tile.replace(".", "-")
    roots = [r for r in (extra_work, run_dir / "work", run_dir / "products") if r is not None]
    for root in roots:
        p = root / f"tiles/{tile[:2]}/{tile}/output/run_sp_tile_Sx/sextractor_runner/output/sexcat-{dashed}.fits"
        if p.exists():
            return p
    return None


def load_exposures(run_dir: Path, tiles: list[str]) -> dict[str, list[str]]:
    db = run_dir / "products/index/run_index.sqlite"
    if not db.exists():
        return {}
    con = sqlite3.connect(db)
    out = {}
    for t in tiles:
        rows = con.execute("SELECT exp_id FROM tile_exposures WHERE tile_id=? ORDER BY exp_id", (t,)).fetchall()
        out[t] = [r[0] for r in rows]
    con.close()
    return out


def load_stars(run_dir: Path, exposures: list[str] | None = None) -> dict[str, dict[str, np.ndarray]]:
    path = find_starcat(run_dir)
    if path is None:
        return {}
    out = {}
    with h5py.File(path, "r") as f:
        keys = exposures if exposures is not None else sorted(f["exposures"].keys())
        for e in keys:
            if str(e) not in f["exposures"]:
                continue
            d = f["exposures"][str(e)][()]
            out[str(e)] = {name: np.asarray(d[name]) for name in d.dtype.names}
    return out


def load_final_cat(run_dir: Path, tile: str) -> dict[str, np.ndarray]:
    path = find_final_cat(run_dir, tile)
    if path is None:
        return {}
    with h5py.File(path, "r") as f:
        return {k: f[k][()] for k in f.keys() if f[k].ndim == 1}


def load_sexcat(run_dir: Path, tile: str, extra_work: Path | None = None) -> tuple[dict[str, np.ndarray], dict[str, tuple]]:
    from astropy.io import fits

    path = find_sexcat(run_dir, tile, extra_work=extra_work)
    if path is None:
        return {}, {}
    scalars, arrays = {}, {}
    with fits.open(path) as hdul:
        data = hdul[2].data
        for name in data.columns.names:
            arr = np.asarray(data[name])
            if arr.ndim == 1:
                scalars[name] = arr
            else:
                arrays[name] = arr.shape
    return scalars, arrays


# --- metrics & comparisons --------------------------------------------------

def tile_of(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    y = np.rint((dec + 90) * 2).astype(int)
    x = np.rint(ra * np.cos(np.radians(y / 2 - 90)) * 2).astype(int) % 720
    return x * 1000 + y


def tile_key(tile_id: np.ndarray) -> np.ndarray:
    if tile_id.dtype.kind in "SUO":
        vals = np.array([float(v.decode() if isinstance(v, bytes) else v) for v in tile_id])
    else:
        vals = tile_id.astype(float)
    return np.rint(vals * 1000).astype(int)


def select_galaxies(cat: dict[str, np.ndarray]) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = cat["NGMIX_FLUX_NOSHEAR"] / cat["NGMIX_FLUX_ERR_NOSHEAR"]
        rel = cat["NGMIX_T_NOSHEAR"] / cat["NGMIX_T_PSF_RECONV_NOSHEAR"]
    own = tile_of(cat["XWIN_WORLD"], cat["YWIN_WORLD"]) == tile_key(cat["TILE_ID"])
    return (
        ((cat["FLAGS"].astype(int) & ~16) <= 2)
        & own
        & ~cat["MASK_4_Stars"].astype(bool)
        & ~cat["MASK_8_Manual"].astype(bool)
        & ~cat["MASK_1024_Maximask"].astype(bool)
        & ~cat["MASK_64_r"].astype(bool)
        & (cat["N_EPOCH"] >= 2)
        & (cat["NGMIX_MCAL_TYPES_FAIL"] == 0)
        & (cat["NGMIX_N_EPOCH"] >= 1)
        & (cat["NGMIX_G1_PSF_ORIG_NOSHEAR"] != -10)
        & (cat["NGMIX_G2_PSF_ORIG_NOSHEAR"] != -10)
        & (cat["NGMIX_FLAGS_NOSHEAR"] == 0)
        & (snr > 10) & (snr < 500)
        & (rel > 0.707) & (rel < 9)
    )


def compare_columns(a: dict[str, np.ndarray], b: dict[str, np.ndarray], key: str = "NUMBER") -> dict:
    """Compare two catalogues matched by `key` (or by row index if key absent)."""
    if key in a and key in b:
        ka, kb = a[key], b[key]
        common_keys, ia, ib = np.intersect1d(ka, kb, return_indices=True)
        n_matched = len(common_keys)
    else:
        n_matched = min(len(next(iter(a.values()))), len(next(iter(b.values()))))
        ia = ib = np.arange(n_matched)

    cols_a, cols_b = set(a.keys()), set(b.keys())
    shared = sorted(cols_a & cols_b - {key})
    changed, identical = [], 0
    for c in shared:
        xa, xb = a[c][ia], b[c][ib]
        if xa.dtype.kind not in "fiub" or xb.dtype.kind not in "fiub":
            continue
        fa, fb = xa.astype(float), xb.astype(float)
        same = (fa == fb) | (np.isnan(fa) & np.isnan(fb))
        if np.all(same):
            identical += 1
        else:
            diff = np.abs(fa - fb)[~same]
            changed.append({
                "column": c,
                "n_diff": int((~same).sum()),
                "frac_identical": round(float(same.mean()), 6),
                "max_abs_diff": float(np.nanmax(diff)) if len(diff) else 0.0,
                "median_abs_diff": float(np.nanmedian(diff)) if len(diff) else 0.0,
            })

    changed.sort(key=lambda r: (r["frac_identical"], -r["max_abs_diff"]))
    return {
        "n_run": int(len(next(iter(a.values())))) if a else 0,
        "n_base": int(len(next(iter(b.values())))) if b else 0,
        "n_matched": int(n_matched),
        "cols_only_in_run": sorted(cols_a - cols_b),
        "cols_only_in_base": sorted(cols_b - cols_a),
        "n_numeric_identical": identical,
        "n_numeric_changed": len(changed),
        "changed_columns": changed,
    }


def star_metrics(stars_by_exp: dict[str, dict[str, np.ndarray]]) -> dict:
    per_exp = {}
    for e, d in sorted(stars_by_exp.items()):
        ok = (d["HSM_FLAG_PSF"] == 0) & (d["HSM_FLAG_STAR"] == 0)
        n = int(ok.sum())
        if n == 0:
            per_exp[e] = {"n_good": 0, "n_total": int(len(ok))}
            continue
        t_s, t_p = d["HSM_T_STAR"][ok], d["HSM_T_PSF"][ok]
        dT = (t_s - t_p) / t_s
        de1 = d["HSM_G1_STAR"][ok] - d["HSM_G1_PSF"][ok]
        de2 = d["HSM_G2_STAR"][ok] - d["HSM_G2_PSF"][ok]
        fwhm = float(SIGMA_TO_FWHM * np.sqrt(np.median(t_p) / 2.0))
        per_exp[e] = {
            "n_total": int(len(ok)),
            "n_good": n,
            "fwhm_arcsec": round(fwhm, 3),
            "dT_mean": [float(dT.mean()), float(dT.std() / np.sqrt(n))],
            "de1_mean": [float(de1.mean()), float(de1.std() / np.sqrt(n))],
            "de2_mean": [float(de2.mean()), float(de2.std() / np.sqrt(n))],
        }
    return per_exp


def compare_stars(run_stars: dict, base_stars: dict) -> dict:
    shared_exps = sorted(set(run_stars) & set(base_stars))
    per_exp = {}
    any_diff = False
    for e in shared_exps:
        ra, rb = run_stars[e], base_stars[e]
        # Match by (CCD_NB, RA, DEC) ordering via row index when counts match, else report count diff
        n_a, n_b = len(ra["RA"]), len(rb["RA"])
        if n_a == n_b:
            col_diff = compare_columns(ra, rb, key="__index__")
        else:
            col_diff = {"n_run": n_a, "n_base": n_b, "n_matched": 0, "n_numeric_changed": -1, "changed_columns": []}
        if col_diff["n_numeric_changed"] != 0 or n_a != n_b:
            any_diff = True
        per_exp[e] = col_diff
    return {"shared_exposures": shared_exps, "any_diff": any_diff, "per_exposure": per_exp}


# --- figures ----------------------------------------------------------------

def plot_psf(run_stars: dict, base_stars: dict | None, out_path: Path) -> None:
    exps = sorted(run_stars.keys())
    if not exps:
        return
    fig, ax = plt.subplots(3, 2, figsize=(12, 9))
    metrics = [("dT", "(T★ − T_PSF) / T★"), ("de1", "g1★ − g1_PSF"), ("de2", "g2★ − g2_PSF")]
    mag_edges = np.arange(18.0, 22.01, 0.5)
    mag_mid = 0.5 * (mag_edges[:-1] + mag_edges[1:])

    def extract(stars_dict):
        mags, dT, de1, de2 = [], [], [], []
        per_e = {}
        for e in exps:
            if e not in stars_dict:
                continue
            d = stars_dict[e]
            ok = (d["HSM_FLAG_PSF"] == 0) & (d["HSM_FLAG_STAR"] == 0)
            v_dT = (d["HSM_T_STAR"][ok] - d["HSM_T_PSF"][ok]) / d["HSM_T_STAR"][ok]
            v_de1 = d["HSM_G1_STAR"][ok] - d["HSM_G1_PSF"][ok]
            v_de2 = d["HSM_G2_STAR"][ok] - d["HSM_G2_PSF"][ok]
            per_e[e] = {"dT": v_dT, "de1": v_de1, "de2": v_de2}
            mags.append(d["MAG"][ok])
            dT.append(v_dT); de1.append(v_de1); de2.append(v_de2)
        if not mags:
            return None, {}
        return {
            "MAG": np.concatenate(mags),
            "dT": np.concatenate(dT),
            "de1": np.concatenate(de1),
            "de2": np.concatenate(de2),
        }, per_e

    series = [("current", run_stars, -0.06, "C0")]
    if base_stars:
        series.insert(0, ("baseline", base_stars, 0.06, "C1"))

    x_idx = np.arange(len(exps))
    exp_labels = []
    sm = star_metrics(run_stars)
    for e in exps:
        fwhm = sm[e].get("fwhm_arcsec", float("nan"))
        exp_labels.append(f"{e}\n({fwhm:.2f}\")")

    for label, sdict, dx, color in series:
        cat, per_e = extract(sdict)
        if cat is None:
            continue
        for i, (key, ylabel) in enumerate(metrics):
            # Left: per exposure
            means, errs = [], []
            for e in exps:
                v = per_e.get(e, {}).get(key, np.array([]))
                means.append(v.mean() if len(v) else np.nan)
                errs.append(v.std() / np.sqrt(len(v)) if len(v) > 1 else np.nan)
            ax[i, 0].errorbar(x_idx + dx, means, errs, fmt="o", ms=5, color=color, label=label)
            ax[i, 0].axhline(0, color="gray", lw=0.6)
            ax[i, 0].set_ylabel(ylabel)
            ax[i, 0].set_xticks(x_idx)
            ax[i, 0].set_xticklabels(exp_labels)

            # Right: vs MAG
            b_m, b_e = [], []
            idx = np.digitize(cat["MAG"], mag_edges) - 1
            for k in range(len(mag_mid)):
                v = cat[key][idx == k]
                b_m.append(v.mean() if len(v) > 5 else np.nan)
                b_e.append(v.std() / np.sqrt(len(v)) if len(v) > 5 else np.nan)
            ax[i, 1].errorbar(mag_mid + dx * 0.5, b_m, b_e, fmt="o", ms=4, color=color, label=label)
            ax[i, 1].axhline(0, color="gray", lw=0.6)

    ax[2, 0].set_xlabel("Exposure (median FWHM)")
    ax[2, 1].set_xlabel("Star r mag")
    ax[0, 0].legend(frameon=False)
    fig.suptitle("ShapePipe verify — reserve-star PSF residuals")
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def plot_all(
    run_stars: dict,
    base_stars: dict | None,
    run_gal: dict[str, np.ndarray],
    base_gal: dict[str, np.ndarray] | None,
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))

    # Top row: reserve-star residuals per exposure
    exps = sorted(run_stars.keys())
    x_idx = np.arange(len(exps))
    sm = star_metrics(run_stars)
    exp_labels = [f"{e}\n({sm[e].get('fwhm_arcsec', float('nan')):.2f}\")" for e in exps]
    series = [("current", run_stars, -0.07, "C0")]
    if base_stars:
        series.insert(0, ("baseline", base_stars, 0.07, "C1"))
    for j, (key, title) in enumerate([("dT", "(T★ − T_PSF) / T★"), ("de1", "g1★ − g1_PSF"), ("de2", "g2★ − g2_PSF")]):
        for label, sdict, dx, color in series:
            m_dict = star_metrics(sdict)
            ys = [m_dict.get(e, {}).get(f"{key}_mean", [np.nan, np.nan])[0] for e in exps]
            es = [m_dict.get(e, {}).get(f"{key}_mean", [np.nan, np.nan])[1] for e in exps]
            ax[0, j].errorbar(x_idx + dx, ys, es, fmt="o", ms=5, color=color, label=label)
        ax[0, j].axhline(0, color="gray", lw=0.6)
        ax[0, j].set_title(title)
        ax[0, j].set_xticks(x_idx)
        ax[0, j].set_xticklabels(exp_labels)
    ax[0, 0].legend(frameon=False)

    # Bottom row: galaxy selection distributions + object-matched Δg
    g_series = [("current", run_gal, "C0")]
    if base_gal:
        g_series.insert(0, ("baseline", base_gal, "C1"))
    for label, cat, color in g_series:
        sel = select_galaxies(cat)
        snr = cat["NGMIX_FLUX_NOSHEAR"][sel] / cat["NGMIX_FLUX_ERR_NOSHEAR"][sel]
        rel = cat["NGMIX_T_NOSHEAR"][sel] / cat["NGMIX_T_PSF_RECONV_NOSHEAR"][sel]
        ax[1, 0].hist(np.log10(snr), bins=35, histtype="step", color=color, label=f"{label} ({int(sel.sum()):,})")
        ax[1, 1].hist(rel, bins=np.linspace(0.6, 5.0, 35), histtype="step", color=color, label=label)
    ax[1, 0].set_xlabel("log10 S/N (selected galaxies)"); ax[1, 0].legend(frameon=False)
    ax[1, 1].set_xlabel("T / T_PSF (selected galaxies)")

    if base_gal:
        sel_r = select_galaxies(run_gal)
        sel_b = select_galaxies(base_gal)
        _, ia, ib = np.intersect1d(run_gal["NUMBER"][sel_r], base_gal["NUMBER"][sel_b], return_indices=True)
        dg1 = run_gal["NGMIX_G1_NOSHEAR"][sel_r][ia] - base_gal["NGMIX_G1_NOSHEAR"][sel_b][ib]
        dg2 = run_gal["NGMIX_G2_NOSHEAR"][sel_r][ia] - base_gal["NGMIX_G2_NOSHEAR"][sel_b][ib]
        if np.any(dg1 != 0) or np.any(dg2 != 0):
            lim = max(0.02, float(np.percentile(np.abs(np.concatenate([dg1, dg2])), 99.5)))
            bins = np.linspace(-lim, lim, 61)
            ax[1, 2].hist(dg1, bins=bins, histtype="step", label=f"Δg1 (rms={np.sqrt(np.mean(dg1**2)):.2e})")
            ax[1, 2].hist(dg2, bins=bins, histtype="step", label=f"Δg2 (rms={np.sqrt(np.mean(dg2**2)):.2e})")
            ax[1, 2].set_yscale("log")
            ax[1, 2].legend(frameon=False)
        else:
            ax[1, 2].text(0.5, 0.5, f"Identical on all {len(ia):,}\nselected galaxies",
                          ha="center", va="center", transform=ax[1, 2].transAxes, fontsize=11)
        ax[1, 2].set_xlabel("Object-matched Δg (current − baseline)")
    else:
        sel_r = select_galaxies(run_gal)
        ax[1, 2].hist(run_gal["MAG_AUTO"][sel_r], bins=np.linspace(18, 26, 35), histtype="step", color="C0")
        ax[1, 2].set_xlabel("MAG_AUTO (selected galaxies)")

    fig.suptitle("ShapePipe verify — end-to-end star & galaxy summary")
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def plot_detect(run_sx: dict, base_sx: dict | None, out_path: Path) -> None:
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    series = [("current", run_sx, "C0")]
    if base_sx:
        series.insert(0, ("baseline", base_sx, "C1"))
    for label, cat, color in series:
        mag = cat["MAG_AUTO"]
        fr = cat["FLUX_RADIUS"]
        ax[0].hist(mag[(mag > 14) & (mag < 30)], bins=40, histtype="step", color=color,
                   label=f"{label} ({len(mag):,})")
        ax[1].hist(fr[(fr > 0) & (fr < 15)], bins=40, histtype="step", color=color, label=label)
        flags = cat["FLAGS"].astype(int)
        vals, cnts = np.unique(flags, return_counts=True)
        top = np.argsort(-cnts)[:8]
        ax[2].plot(vals[top], cnts[top], "o", color=color, label=label)
    ax[0].set_xlabel("MAG_AUTO"); ax[0].legend(frameon=False)
    ax[1].set_xlabel("FLUX_RADIUS [px]")
    ax[2].set_xlabel("FLAGS"); ax[2].set_yscale("log")
    fig.suptitle("ShapePipe verify — tile_detect catalogue")
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


# --- reporting --------------------------------------------------------------

def inspect_and_report(
    stage: str,
    tiles: list[str],
    run_dir: Path,
    base_dir: Path | None,
    out_dir: Path,
    run_meta: dict,
    git_meta: dict,
    work_dir: Path | None = None,
) -> tuple[dict, bool]:
    summary: dict = {"git": git_meta, "run": run_meta, "baseline": str(base_dir) if base_dir else None}
    any_diff = False

    # 1. Prepare / exposure index
    exps_by_tile = load_exposures(run_dir, tiles)
    all_exps = sorted({e for es in exps_by_tile.values() for e in es})
    summary["exposures_by_tile"] = exps_by_tile
    if base_dir:
        base_exps = load_exposures(base_dir, tiles)
        if base_exps and base_exps != exps_by_tile:
            any_diff = True
            summary["exposures_diff"] = {"run": exps_by_tile, "base": base_exps}

    # 2. Detect
    if stage in ("detect", "vignets"):
        tile = tiles[0]
        run_sx, run_arr = load_sexcat(run_dir, tile, extra_work=work_dir)
        base_sx, base_arr = load_sexcat(base_dir, tile) if base_dir else ({}, {})
        summary["detect"] = {
            "tile": tile,
            "n_objects": int(len(run_sx.get("NUMBER", []))),
            "n_flags_zero": int((run_sx["FLAGS"] == 0).sum()) if "FLAGS" in run_sx else 0,
            "array_columns": {k: list(v) for k, v in run_arr.items()},
        }
        if base_sx:
            dcmp = compare_columns(run_sx, base_sx, key="NUMBER")
            summary["detect"]["comparison"] = dcmp
            if dcmp["n_numeric_changed"] > 0 or dcmp["cols_only_in_run"] or dcmp["cols_only_in_base"] or run_arr != base_arr:
                any_diff = True
        if run_sx and stage == "detect":
            fig_path = out_dir / "verify_detect.png"
            plot_detect(run_sx, base_sx or None, fig_path)
            summary["figure"] = str(fig_path)

    # 3. PSF / stars
    run_stars = load_stars(run_dir, all_exps or None)
    base_stars = load_stars(base_dir, sorted(run_stars.keys())) if (base_dir and run_stars) else {}
    if run_stars:
        summary["stars"] = star_metrics(run_stars)
        if base_stars:
            summary["stars_baseline"] = star_metrics(base_stars)
            scmp = compare_stars(run_stars, base_stars)
            summary["stars_comparison"] = scmp
            if scmp["any_diff"]:
                any_diff = True
        if stage in ("psf", "vignets"):
            fig_path = out_dir / f"verify_{stage}.png"
            plot_psf(run_stars, base_stars or None, fig_path)
            summary["figure"] = str(fig_path)

    # 4. Final galaxy catalogue (stage == all)
    if stage == "all":
        tile = tiles[0]
        run_gal = load_final_cat(run_dir, tile)
        base_gal = load_final_cat(base_dir, tile) if base_dir else {}
        if run_gal:
            sel_r = select_galaxies(run_gal)
            r11 = float(((run_gal["NGMIX_G1_1P"][sel_r] - run_gal["NGMIX_G1_1M"][sel_r]) / 0.02).mean())
            r22 = float(((run_gal["NGMIX_G2_2P"][sel_r] - run_gal["NGMIX_G2_2M"][sel_r]) / 0.02).mean())
            snr_r = run_gal["NGMIX_FLUX_NOSHEAR"][sel_r] / run_gal["NGMIX_FLUX_ERR_NOSHEAR"][sel_r]
            summary["galaxies"] = {
                "tile": tile,
                "n_objects": int(len(run_gal["NUMBER"])),
                "n_columns": len(run_gal),
                "n_selected": int(sel_r.sum()),
                "R11": round(r11, 5),
                "R22": round(r22, 5),
                "snr_median": round(float(np.median(snr_r)), 2) if sel_r.any() else None,
            }
            if base_gal:
                gcmp = compare_columns(run_gal, base_gal, key="NUMBER")
                sel_b = select_galaxies(base_gal)
                _, ia, ib = np.intersect1d(run_gal["NUMBER"][sel_r], base_gal["NUMBER"][sel_b], return_indices=True)
                dg1 = run_gal["NGMIX_G1_NOSHEAR"][sel_r][ia] - base_gal["NGMIX_G1_NOSHEAR"][sel_b][ib]
                dg2 = run_gal["NGMIX_G2_NOSHEAR"][sel_r][ia] - base_gal["NGMIX_G2_NOSHEAR"][sel_b][ib]
                gcmp["selected_baseline"] = int(sel_b.sum())
                gcmp["selected_matched"] = int(len(ia))
                gcmp["mean_dg1"] = float(dg1.mean()) if len(dg1) else 0.0
                gcmp["mean_dg2"] = float(dg2.mean()) if len(dg2) else 0.0
                gcmp["rms_dg"] = float(np.sqrt((dg1**2 + dg2**2).mean() / 2.0)) if len(dg1) else 0.0
                summary["galaxies_comparison"] = gcmp
                if gcmp["n_numeric_changed"] > 0 or gcmp["cols_only_in_run"] or gcmp["cols_only_in_base"]:
                    any_diff = True
            fig_path = out_dir / "verify_all.png"
            plot_all(run_stars, base_stars or None, run_gal, base_gal or None, fig_path)
            summary["figure"] = str(fig_path)

    summary["verdict"] = (
        "NO_BASELINE" if base_dir is None
        else ("DIFFERENCES_DETECTED" if any_diff else "NO_DIFFERENCE")
    )
    (out_dir / "verify_summary.json").write_text(json.dumps(summary, indent=2))
    print_report(summary)
    return summary, any_diff


def print_report(s: dict) -> None:
    print("=" * 72)
    g = s["git"]
    r = s["run"]
    print(f"ShapePipe verify | stage={r.get('stage')} | tiles={','.join(r.get('tiles', []))} | "
          f"branch={g.get('branch')} ({g.get('commit', '')[:8]}{'*' if g.get('dirty') else ''})")
    if "seconds" in r:
        print(f"Runtime: {r['seconds']}s on {r['cores']} cores ({r['core_hours']} core-h)")
    print(f"Verdict: {s['verdict']}" + (f" (baseline: {s['baseline']})" if s["baseline"] else ""))
    print("-" * 72)

    if "detect" in s:
        d = s["detect"]
        print(f"tile_detect ({d['tile']}): {d['n_objects']:,} objects ({d['n_flags_zero']:,} FLAGS==0)")
        if "comparison" in d:
            c = d["comparison"]
            print(f"  vs baseline: matched={c['n_matched']:,}  identical_cols={c['n_numeric_identical']}  "
                  f"changed_cols={c['n_numeric_changed']}")
            for row in c["changed_columns"][:8]:
                print(f"    {row['column']:28s} frac_id={row['frac_identical']:.5f} max_diff={row['max_abs_diff']:.3e}")

    if "stars" in s:
        print("Reserve stars (per exposure):")
        sb = s.get("stars_baseline", {})
        for e, m in s["stars"].items():
            line = (f"  exp {e} (FWHM={m.get('fwhm_arcsec', 0):.2f}\", n={m['n_good']:,}): "
                    f"dT={m['dT_mean'][0]:+.5f}±{m['dT_mean'][1]:.5f}  "
                    f"de1={m['de1_mean'][0]:+.5f}±{m['de1_mean'][1]:.5f}  "
                    f"de2={m['de2_mean'][0]:+.5f}±{m['de2_mean'][1]:.5f}")
            if e in sb:
                bm = sb[e]
                line += (f"  [base: dT={bm['dT_mean'][0]:+.5f} "
                         f"de1={bm['de1_mean'][0]:+.5f} de2={bm['de2_mean'][0]:+.5f}]")
            print(line)
        if "stars_comparison" in s:
            for e, c in s["stars_comparison"]["per_exposure"].items():
                print(f"  exp {e} columns: identical={c.get('n_numeric_identical', 0)} changed={c['n_numeric_changed']}")
                for row in c["changed_columns"][:5]:
                    print(f"    {row['column']:28s} frac_id={row['frac_identical']:.5f} max_diff={row['max_abs_diff']:.3e}")

    if "galaxies" in s:
        gm = s["galaxies"]
        print(f"Final catalogue ({gm['tile']}): {gm['n_objects']:,} objects, {gm['n_columns']} cols, "
              f"{gm['n_selected']:,} selected (R11={gm['R11']:.4f}, R22={gm['R22']:.4f}, med S/N={gm['snr_median']})")
        if "galaxies_comparison" in s:
            c = s["galaxies_comparison"]
            print(f"  vs baseline: matched={c['n_matched']:,}  identical_cols={c['n_numeric_identical']}  "
                  f"changed_cols={c['n_numeric_changed']}  "
                  f"selected={gm['n_selected']:,} vs {c['selected_baseline']:,}")
            if c["cols_only_in_run"]:
                print(f"  cols only in run : {', '.join(c['cols_only_in_run'])}")
            if c["cols_only_in_base"]:
                print(f"  cols only in base: {', '.join(c['cols_only_in_base'])}")
            print(f"  selected paired: mean Δg=({c['mean_dg1']:+.3e}, {c['mean_dg2']:+.3e})  rms(Δg)={c['rms_dg']:.3e}")
            for row in c["changed_columns"][:10]:
                print(f"    {row['column']:28s} frac_id={row['frac_identical']:.5f} "
                      f"n_diff={row['n_diff']:5d} max_diff={row['max_abs_diff']:.3e}")

    if "figure" in s:
        print(f"Figure: {s['figure']}")
    print("=" * 72)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=STAGES, default="all", help="Stage to verify (default: all)")
    ap.add_argument("--tiles", nargs="+", default=DEFAULT_TILES, help="Fixture tile ID(s) (default: 592.232)")
    ap.add_argument("--inputs", default=default_path("SP_VERIFY_INPUTS", "inputs"), help="Staged inputs root")
    ap.add_argument("--masks", default=default_path("SP_VERIFY_MASKS", "masks"), help="External sky masks dir")
    ap.add_argument("--container", default=default_path("SP_CONTAINER", "container"), help="ShapePipe SIF image")
    ap.add_argument("--baseline", default=default_path("SP_VERIFY_BASELINE", "baseline"), help="Known-good run dir")
    ap.add_argument("--no-baseline", action="store_true", help="Skip baseline comparison (e.g. when creating a baseline)")
    ap.add_argument("--out", default=None, help="Output directory")
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]), help="ShapePipe checkout to run")
    ap.add_argument("--overrides", default=None, help="Optional YAML overrides for run config")
    ap.add_argument("--cores", type=int, default=None, help="Cores for local execution")
    ap.add_argument("--ngmix-chunks", type=int, default=None, help="Scatter chunks for tile_ngmix")
    ap.add_argument("--profile", default=os.environ.get("SP_PROFILE", "nibi"), help="Snakemake profile")
    ap.add_argument("--keep-work", action="store_true", help="Keep scratch work directory after run")
    ap.add_argument("--compare-only", default=None, help="Skip running; compare existing run dir against --baseline")
    ap.add_argument("--strict", action="store_true", help="Exit 1 if any difference vs baseline is detected")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    base_dir = None if args.no_baseline else (Path(args.baseline).resolve() if args.baseline else None)

    if args.compare_only:
        run_dir = Path(args.compare_only).resolve()
        out_dir = Path(args.out).resolve() if args.out else run_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        _, any_diff = inspect_and_report(
            args.stage, args.tiles, run_dir, base_dir, out_dir,
            {"stage": args.stage, "tiles": args.tiles, "compare_only": True},
            git_info(repo),
        )
        sys.exit(1 if (args.strict and any_diff) else 0)

    for name, val in (("inputs", args.inputs), ("masks", args.masks), ("container", args.container)):
        if not val or not Path(val).exists():
            sys.exit(f"verify: --{name} not found ({val!r}); pass --{name} or set the env var.")

    cores = args.cores or min(32, int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 8))
    ngmix_chunks = args.ngmix_chunks or min(16, cores)
    tmp_root = Path(os.environ.get("SLURM_TMPDIR") or "/tmp")
    out_dir = Path(args.out).resolve() if args.out else (tmp_root / f"sp-verify-{args.stage}-{int(time.time())}")
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = out_dir / "work" if args.keep_work else (tmp_root / f"sp-verify-work-{os.getpid()}")
    work_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / "tiles.txt").write_text("\n".join(args.tiles) + "\n")
    cfg_path = build_config(
        out_dir, work_dir, Path(args.inputs), Path(args.masks),
        args.container, ngmix_chunks,
        Path(args.overrides) if args.overrides else None, args.profile,
    )

    run_meta = run_launcher(
        args.stage, args.tiles, repo, out_dir, work_dir,
        Path(args.inputs), Path(args.masks), args.container,
        cfg_path, cores, args.profile,
    )
    _, any_diff = inspect_and_report(
        args.stage, args.tiles, out_dir, base_dir, out_dir, run_meta, git_info(repo), work_dir=work_dir,
    )

    if not args.keep_work and work_dir != out_dir / "work":
        shutil.rmtree(work_dir, ignore_errors=True)

    sys.exit(1 if (args.strict and any_diff) else 0)


if __name__ == "__main__":
    main()
