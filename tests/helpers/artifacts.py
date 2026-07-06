"""Artifact emission for guardrail tests — the GitHub Pages publish SEAM.

A guardrail test calls :func:`emit_star_response_artifacts` to drop three
files into ``tests/_artifacts/``:

* ``<name>.png``   — the R1/R2 distribution + jackknife-error figure.
* ``<name>.json``  — machine-readable status (metric values, errors, pass/fail).
* ``<name>.md``    — a short human-readable status summary.

This module *only emits*. A later CI step reads ``tests/_artifacts/`` and
publishes to Pages — that deploy is deliberately NOT built here. Keeping the
emit and the publish separate means the test can run (and the artifacts are
inspectable) with zero web infrastructure.

Plotting follows the repo's seaborn-first convention but degrades gracefully:
if seaborn/matplotlib are unavailable the figure is skipped and the JSON/MD
status still lands, so the seam never blocks a result from being recorded.
"""

import json
from datetime import datetime, timezone
from pathlib import Path


def _plot_star_response(png_path, R1, R2, summary):
    """R1/R2 histograms with jackknife mean +/- error bands. Best-effort."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        try:
            import seaborn as sns

            sns.set_theme(style="whitegrid")
            sns.set_palette("husl", 2)
        except ImportError:
            pass

        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for ax, R, label, mean, err in (
            (axes[0], R1, "R1", summary["R1_mean"], summary["R1_jk_err"]),
            (axes[1], R2, "R2", summary["R2_mean"], summary["R2_jk_err"]),
        ):
            R = np.asarray(R)
            # Clip the view to the informative core: per-object R has heavy
            # tails (fit failures reach |R| ~ 80) that would crush the bulk.
            lo, hi = np.percentile(R, [1, 99])
            ax.hist(R, bins=80, range=(lo, hi), histtype="stepfilled", alpha=0.5)
            ax.axvline(0.0, color="k", lw=1, ls="--", label="target 0")
            ax.axvline(mean, color="C3", lw=2, label=f"mean {label} = {mean:.4f}")
            ax.axvspan(mean - err, mean + err, color="C3", alpha=0.2,
                       label=f"jk +/- {err:.4f}")
            ax.set_xlim(lo, hi)
            ax.set_xlabel(label)
            ax.set_title(f"mean {label} = {mean:.4f} +/- {err:.4f}")
            ax.legend(fontsize=8)
        axes[0].set_ylabel("count")
        fig.suptitle(
            f"Star shear response  (N={summary['n_obj']}, "
            f"tol={summary['tolerance']})  -- {summary['status']}"
        )
        fig.tight_layout()
        fig.savefig(png_path, dpi=120)
        plt.close(fig)
        return str(png_path)
    except Exception as exc:  # noqa: BLE001 — emit must never crash the test
        return f"plot skipped: {exc!r}"


def emit_star_response_artifacts(artifacts_dir, summary, R1=None, R2=None,
                                 name="star_shear_response"):
    """Write JSON + markdown status and (when data given) the figure.

    Parameters
    ----------
    artifacts_dir : Path
        Target directory (the ``artifacts_dir`` fixture).
    summary : dict
        Must carry ``R1_mean``, ``R2_mean``, ``R1_jk_err``, ``R2_jk_err``,
        ``n_obj``, ``tolerance``, ``status`` (``"PASS"``/``"FAIL"``), and
        ``base_dir``.
    R1, R2 : array-like, optional
        Per-object responses for the histogram. Plot skipped if omitted.
    name : str
        Artifact basename.

    Returns
    -------
    dict
        Paths of the artifacts written (``json``, ``md``, ``png``).
    """
    artifacts_dir = Path(artifacts_dir)
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    summary = {**summary, "generated_utc": datetime.now(timezone.utc).isoformat()}

    # Status JSON is scalar-only: the per-object R1/R2 arrays are for the plot,
    # not the machine-readable summary a Pages step ingests.
    status_payload = {
        k: v for k, v in summary.items() if k not in ("R1", "R2")
    }
    json_path = artifacts_dir / f"{name}.json"
    json_path.write_text(json.dumps(status_payload, indent=2, default=str) + "\n")

    md_path = artifacts_dir / f"{name}.md"
    md_path.write_text(
        f"# Star shear-response guardrail — {summary['status']}\n\n"
        f"- generated: `{summary['generated_utc']}`\n"
        f"- outputs:   `{summary['base_dir']}`\n"
        f"- N objects: {summary['n_obj']}\n"
        f"- tolerance: |R| < {summary['tolerance']}\n\n"
        f"| metric | value | jackknife error |\n"
        f"|--------|-------|-----------------|\n"
        f"| `<R1>` | {summary['R1_mean']:.6f} | {summary['R1_jk_err']:.6f} |\n"
        f"| `<R2>` | {summary['R2_mean']:.6f} | {summary['R2_jk_err']:.6f} |\n\n"
        f"Target (Fabian, 2026-06-12): `R1 = R2 = 0 +/- 0.03`. "
        f"Deconvolution should leave no net star shear response.\n"
    )

    png_path = artifacts_dir / f"{name}.png"
    plot_result = (
        _plot_star_response(png_path, R1, R2, summary)
        if R1 is not None and R2 is not None
        else "plot skipped: no per-object data passed"
    )

    return {"json": str(json_path), "md": str(md_path), "png": plot_result}


def emit_mbias_artifacts(artifacts_dir, summary, name="mbias"):
    """Write the multiplicative-bias guardrail status (JSON + markdown).

    The sibling of :func:`emit_star_response_artifacts` for the Level-1a m-bias
    test. The metric is a single scalar per run (no per-object array), so this
    emits the scalar status only — the status-board figure for this rung is drawn
    separately from the recorded values. Same seam: files land in
    ``tests/_artifacts/`` on pass *and* fail, and a later Pages step publishes.

    Parameters
    ----------
    artifacts_dir : Path
        Target directory (typically ``tests/_artifacts/``).
    summary : dict
        Must carry ``injected_g1``, ``recovered_g1``, ``metacal_response_R11``,
        ``multiplicative_bias_m``, ``tolerance``, and ``status``
        (``"PASS"``/``"FAIL"``).
    name : str
        Artifact basename.

    Returns
    -------
    dict
        Paths of the artifacts written (``json``, ``md``).
    """
    artifacts_dir = Path(artifacts_dir)
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    summary = {**summary, "generated_utc": datetime.now(timezone.utc).isoformat()}

    json_path = artifacts_dir / f"{name}.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str) + "\n")

    md_path = artifacts_dir / f"{name}.md"
    md_path.write_text(
        f"# Multiplicative-bias guardrail — {summary['status']}\n\n"
        f"- generated: `{summary['generated_utc']}`\n"
        f"- injected g1: {summary['injected_g1']}\n"
        f"- metacal response R11: {summary['metacal_response_R11']:.6f}\n\n"
        f"| metric | value | tolerance |\n"
        f"|--------|-------|-----------|\n"
        f"| `m` | {summary['multiplicative_bias_m']:+.6f} | "
        f"|m| < {summary['tolerance']} |\n\n"
        f"Ideal-limit multiplicative bias `m = g_rec/g_inj − 1` after the metacal "
        f"response correction. Recovered g1 = {summary['recovered_g1']:.6f} vs "
        f"injected {summary['injected_g1']}.\n"
    )

    return {"json": str(json_path), "md": str(md_path)}
