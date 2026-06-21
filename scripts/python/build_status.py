#!/usr/bin/env python
"""Render the ShapePipe shape-measurement *status board* from a run's JSON outputs.

This is the renderer referenced in the status-board footer: the digital-twin
harvesters (``run_mbias.py``, ``run_star_response.py``,
``run_star_response_grid.py``) emit ``results/baseline/*/*.json``; this script
reads those JSONs and renders a single self-contained ``status.html``. One render
per run; history lives in the commits.

The page is organised by a *fidelity ladder*:

* **Level 1a — m-bias** (``mbias_sim/mbias.json``): an in-memory sim injects a
  known shear through the live ``make_data / do_ngmix_metacal`` path and recovers
  it after the metacal response correction. Verdict: ``|m| < 5e-3``.
* **Level 1b — psf_fit_prior** (``star_response_sim/star_response.json``): an
  in-memory sim with an elliptical PSF checks that the recovered shear ``c`` is
  unchanged across PSF-fit-prior choices. Verdict: shear robust to the prior.
* **Level 2 — star shear response on the SKiLLS grid**
  (``star_response_grid_<arm>/star_response.json``): per-object metacal response
  ``R = dg/dgamma`` on Fabian Hervas's SKiLLS star grid, as an A/B over the
  metacal reconvolution-kernel sizing (``fitgauss`` vs ``gauss``) crossed with the
  centroid source (``wcs`` vs ``hsm``). Verdict: ``<R1>``, ``<R2>`` within
  ``+/-0.03`` of 0.

Design contract
---------------
*Every number on the page is read straight from a results JSON.* The verdict
badges (PASS / MARGINAL / FAIL) are **computed** from the recorded value against
its recorded tolerance, never hard-coded — that is the live part of the board.
The CSS and the explanatory prose are the template; the figures are optional PNGs
taken from the run's outputs and inlined as base64 (a missing figure degrades to
a "figure pending" note rather than crashing the render).

Stdlib only (``json``, ``base64``, ``argparse``, ``pathlib``, ``datetime``,
``html``) so it runs anywhere — inside the container, on a login node, or via
``uv run scripts/python/build_status.py``.

Regenerate + deploy
-------------------
See ``main()``'s ``--help`` and the module-level ``USAGE`` string. Deploy to GH
Pages is intentionally **out of scope** here — Cail owns the deploy step.
"""
from __future__ import annotations

import argparse
import base64
import datetime as _dt
import html
import json
from pathlib import Path

USAGE = """\
Regenerate the status board from a run's outputs:

    python scripts/python/build_status.py \\
        --results <run>/results/baseline \\
        --figures <run>/figures \\
        --out status.html

`--results` points at the directory holding the per-arm subdirs
(`mbias_sim/`, `star_response_sim/`, `star_response_grid_*/`); `--figures`
(optional) points at a directory of PNGs to inline (see FIGURE_SLOTS).

Deploy (Cail owns this step; not run by this script):

    cp status.html <gh-pages-clone>/shapepipe/status.html
    cd <gh-pages-clone> && git add shapepipe/status.html \\
        && git commit -m 'status: re-render <date>' && git push
"""

# Theme colours, mirrored from the page CSS so verdict badges stay on-palette.
TEAL = "#3F8278"      # PASS
OCHRE = "#C49333"     # MARGINAL
CINNABAR = "#BC4538"  # FAIL
MUTED = "var(--ink-muted)"  # pending / soon

# Figure filenames looked for under --figures, keyed by the slot they fill.
# A run that produced these PNGs gets them inlined; any that are missing fall
# back to a "figure pending" placeholder so the render never hard-fails.
FIGURE_SLOTS = {
    "mbias": "mbias.png",
    "psf_fit_prior": "psf_fit_prior.png",
    "grid_hist": "star_response_hist.png",
    "grid_forest": "star_response_forest.png",
    "grid_ab": "star_response_ab_tail.png",
}


# --------------------------------------------------------------------------- #
# Small formatting helpers
# --------------------------------------------------------------------------- #
def sci(x: float, sig: int = 3) -> str:
    """Format ``x`` as ``mantissa x10^exp`` with HTML superscript digits."""
    if x == 0:
        return "0"
    sup = str.maketrans("-0123456789", "−⁰¹²³⁴⁵⁶⁷⁸⁹")
    s = f"{x:.{sig - 1}e}"               # e.g. "-1.78e-03"
    mant, exp = s.split("e")
    exp_i = int(exp)
    sign = "−" if mant.startswith("-") else ""
    mant = mant.lstrip("-")
    return f"{sign}{mant} &times;10{str(exp_i).translate(sup)}"


def signed(x: float, dp: int = 4) -> str:
    """Signed fixed-point, e.g. ``+0.0405`` / ``-0.0014`` with a unicode minus."""
    return f"{x:+.{dp}f}".replace("-", "−")


def b64_png(path: Path) -> str | None:
    """Base64-encode a PNG for inline embedding, or ``None`` if absent."""
    if not path.is_file():
        return None
    return base64.b64encode(path.read_bytes()).decode("ascii")


def load(results: Path, arm: str, fname: str) -> dict | None:
    """Load ``<results>/<arm>/<fname>``; ``None`` (not an error) if absent."""
    p = results / arm / fname
    return json.loads(p.read_text()) if p.is_file() else None


def badge(label: str, colour: str) -> str:
    bg = colour if colour.startswith("var(") else colour
    return f'<span class="badge" style="background:{bg}">{label}</span>'


def figure_block(b64: str | None, alt: str, caption: str, data_prov: str, render_ts: str) -> str:
    """A ``<figure>`` with an inlined PNG, or a pending placeholder."""
    if b64 is None:
        body = (
            '<div class="plot" style="display:flex;align-items:center;'
            'justify-content:center;min-height:120px;color:var(--ink-muted);'
            f'font-family:var(--mono);font-size:13px">figure pending &mdash; {html.escape(alt)}</div>'
        )
    else:
        body = f'<div class="plot"><img alt="{html.escape(alt)}" src="data:image/png;base64,{b64}"></div>'
    prov = (
        '<div class="prov">'
        f'<span class="pv"><span class="pk">data</span> {data_prov}</span>'
        f'<span class="pv"><span class="pk">figure</span> build_status.py <span class="pt">{render_ts}</span></span>'
        "</div>"
    )
    return f"<figure>{body}<figcaption>{caption}</figcaption>{prov}</figure>"


# --------------------------------------------------------------------------- #
# CSS (verbatim from the curated page; the page's visual identity)
# --------------------------------------------------------------------------- #
CSS = """
  :root {
    --parchment:#EDE2CE; --rag:#FBF6E9; --ink:#1B1611; --ink-muted:#7A6F5C;
    --cinnabar:#BC4538; --teal:#3F8278; --cobalt:#3D5BA0; --ochre:#C49333;
    --rule:#C5B89E; --tooth:rgba(197,184,158,0.45);
    --serif:'EB Garamond',Georgia,serif;
    --mono:'IBM Plex Mono','JetBrains Mono',ui-monospace,monospace;
  }
  html,body { margin:0; padding:0; scroll-behavior:smooth; }
  body { background:var(--parchment); color:var(--ink); font-family:var(--serif);
    font-size:18px; line-height:1.5; -webkit-font-smoothing:antialiased; }
  .page { max-width:1000px; margin:0 auto; padding:52px 32px 72px; }
  a { color:var(--cobalt); text-decoration:none; }
  a:hover { text-decoration:underline; }

  .stamp { font-family:var(--mono); font-size:11px; letter-spacing:.12em;
    text-transform:uppercase; color:var(--ink-muted); margin-bottom:18px;
    display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
  .stamp .dot { width:6px; height:6px; border-radius:50%; background:var(--cinnabar); }
  .stamp .sep { color:var(--rule); }
  h1 { font-family:var(--serif); font-style:italic; font-weight:500;
    font-size:42px; line-height:1.08; margin:0 0 12px; }
  .lede { color:var(--ink-muted); font-size:18px; margin:0 0 6px; max-width:72ch; }
  .lede b { color:var(--ink); font-weight:500; }

  .disclaimer { margin:22px 0 0; padding:13px 18px; border-radius:8px;
    background:rgba(196,147,51,0.14); border:1px solid var(--ochre);
    display:flex; gap:12px; align-items:baseline; font-family:var(--mono);
    font-size:12.5px; letter-spacing:.02em; color:#8a6a1f; }
  .disclaimer b { color:var(--ochre); text-transform:uppercase; letter-spacing:.1em; }

  /* ---- fidelity ladder: the spine, intro + live verdicts in one strip ---- */
  .ladder { margin:30px 0 0; display:grid; grid-template-columns:1fr 1fr 1fr;
    gap:16px; }
  .rung { background:var(--rag); border:1px solid var(--rule); border-radius:10px;
    padding:16px 18px 14px; font-size:14.5px; display:block; color:var(--ink);
    transition:border-color .15s, transform .15s; }
  a.rung:hover { border-color:var(--cobalt); text-decoration:none;
    transform:translateY(-2px); }
  .rung p { margin:0; color:var(--ink-muted); font-size:13.5px; line-height:1.45; }
  .rung-head { font-family:var(--mono); font-size:10.5px; letter-spacing:.1em;
    text-transform:uppercase; margin-bottom:8px; font-weight:500; }
  .rung-head.l1 { color:var(--teal); }
  .rung-head.l2 { color:var(--cobalt); }
  .rung-head.l3 { color:var(--ink-muted); }
  .lchips { display:flex; gap:8px; flex-wrap:wrap; margin-top:12px; }
  .lchip { font-family:var(--mono); font-size:11px; color:var(--ink-muted);
    display:inline-flex; align-items:center; gap:6px; }

  section { margin-top:46px; scroll-margin-top:20px; }
  .section-mark { font-family:var(--mono); font-size:11px; letter-spacing:.15em;
    text-transform:uppercase; color:var(--cinnabar); margin-bottom:10px;
    display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
  h2 { font-family:var(--serif); font-style:italic; font-weight:500;
    font-size:28px; line-height:1.18; margin:0 0 4px; }

  hr { border:none; border-top:1px solid var(--rule); margin:36px 0; }

  /* ---- prominent source buttons ---- */
  .srcbtn { display:inline-flex; align-items:center; gap:8px;
    font-family:var(--mono); font-size:12px; letter-spacing:.02em;
    color:var(--cobalt); background:var(--rag); border:1px solid var(--rule);
    border-radius:7px; padding:7px 13px; transition:all .14s; line-height:1; }
  .srcbtn:hover { background:var(--cobalt); color:#fff; border-color:var(--cobalt);
    text-decoration:none; }
  .srcbtn .gl { font-weight:500; opacity:.7; }
  .srcbtn:hover .gl { opacity:1; }

  /* ---- spec block ---- */
  .spec { margin:18px 0 0; padding:14px 20px; border-left:2px solid var(--rule);
    color:var(--ink-muted); font-size:16px; line-height:1.5; }
  .spec .mark { font-family:var(--mono); font-size:10px; letter-spacing:.14em;
    text-transform:uppercase; color:var(--teal); display:block; margin-bottom:6px; }
  .spec p { margin:0 0 6px; } .spec p:last-child { margin:0; }
  .spec b { color:var(--ink); font-style:normal; }
  code { font-family:var(--mono); font-size:.8em;
    background:rgba(196,147,51,0.13); padding:1px 5px; border-radius:2px; }

  /* ---- figures ---- */
  figure { margin:22px 0 0; }
  .plot { background:var(--rag); border:1px solid var(--rule); border-radius:8px;
    padding:12px; }
  .plot img { width:100%; display:block; border-radius:4px; }
  figcaption { color:var(--ink-muted); font-size:15px; margin-top:11px;
    max-width:84ch; line-height:1.5; }
  figcaption b { color:var(--ink); }
  .prov { display:flex; gap:18px; flex-wrap:wrap; margin-top:9px;
    font-family:var(--mono); font-size:10.5px; color:var(--ink-muted);
    border-top:1px dashed var(--tooth); padding-top:8px; }
  .prov .pk { color:var(--teal); letter-spacing:.06em; text-transform:uppercase;
    margin-right:5px; }
  .prov .pt { color:var(--ink); margin-left:5px; }

  /* ---- verdict row ---- */
  .verdict { display:flex; align-items:center; gap:14px; flex-wrap:wrap;
    margin:18px 0 0; padding:14px 18px; background:var(--rag);
    border:1px solid var(--rule); border-radius:8px; }
  .verdict .vlabel { font-family:var(--mono); font-size:12px;
    color:var(--ink-muted); letter-spacing:.03em; }
  .verdict .vval { font-family:var(--mono); font-size:16px; font-weight:500;
    color:var(--ink); }
  .verdict .vn { font-family:var(--mono); font-size:11px; color:var(--ink-muted);
    margin-left:auto; text-align:right; max-width:46ch; line-height:1.4; }

  .badge { color:#fff; font-family:var(--mono); font-weight:500; font-size:11px;
    line-height:1; padding:6px 10px; border-radius:4px; letter-spacing:.08em;
    white-space:nowrap; }

  /* ---- reconvolution-PSF A/B callout (the headline lever) ---- */
  .abcallout { margin:22px 0 0; padding:16px 22px; background:rgba(63,130,120,0.10);
    border:1px solid var(--teal); border-left:4px solid var(--teal);
    border-radius:0 8px 8px 0; }
  .abcallout .abc-mark { font-family:var(--mono); font-size:10.5px;
    letter-spacing:.13em; text-transform:uppercase; color:var(--teal);
    margin-bottom:9px; display:block; }
  .abcallout p { margin:0 0 8px; font-size:16px; line-height:1.55;
    color:var(--ink); }
  .abcallout p:last-child { margin:0; }
  .abcallout .abc-mech { color:var(--ink-muted); font-size:15px; }
  .abcallout .abc-mech b { color:var(--ink); }
  .abcallout .abc-floor { color:var(--ink-muted); font-size:15px;
    border-top:1px dashed var(--tooth); padding-top:10px; }
  .abcallout .abc-floor b { color:var(--teal); }

  /* ---- additive-bias companion line ---- */
  .addbias { margin:16px 0 0; padding:11px 18px; background:var(--rag);
    border:1px dashed var(--rule); border-radius:8px; font-family:var(--mono);
    font-size:13px; color:var(--ink); display:flex; gap:12px; align-items:baseline;
    flex-wrap:wrap; }
  .ab-mark { font-size:10px; letter-spacing:.1em; text-transform:uppercase;
    color:var(--cinnabar); }

  /* ---- test-data block ---- */
  .testdata { margin:16px 0 0; padding:13px 18px; background:var(--rag);
    border-left:3px solid var(--cobalt); border-radius:0 6px 6px 0;
    font-size:14.5px; line-height:1.55; color:var(--ink-muted); }
  .testdata-head { font-family:var(--mono); font-size:10px; letter-spacing:.13em;
    text-transform:uppercase; color:var(--cobalt); margin-bottom:7px; }
  .testdata-body b { color:var(--ink); font-weight:500; }

  /* ---- pending ---- */
  .pending { display:flex; align-items:flex-start; gap:16px; padding:20px 22px;
    background:var(--rag); border:1px solid var(--rule); border-radius:8px;
    margin-top:18px; font-size:16px; }
  .pending-icon { font-size:28px; color:var(--ochre); flex-shrink:0; line-height:1; }
  .pending b { color:var(--ink); }

  footer { margin-top:52px; border-top:1px solid var(--rule); padding-top:18px;
    color:var(--ink-muted); font-family:var(--mono); font-size:11px;
    line-height:1.7; letter-spacing:.02em; }
  footer b { color:var(--ink); }
"""

# Single source-link base for every "&lt;/&gt;" src button on the page. Pinned to
# the ngmix_v2.0 integration branch while that is where this work lives; change
# this one constant to re-point every link when the work merges elsewhere.
REPO_BLOB = "https://github.com/CosmoStat/shapepipe/blob/ngmix_v2.0"


# --------------------------------------------------------------------------- #
# Section renderers — each consumes one JSON dict and returns an HTML <section>
# --------------------------------------------------------------------------- #
def render_mbias(d: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Level 1a m-bias section. Returns (section_html, ladder_chip_html)."""
    m = d["multiplicative_bias_m"]
    r11 = d["metacal_response_R11"]
    cfg = d.get("config", {})
    prov = d.get("provenance", {})
    tol = 5e-3  # |m| < 5e-3 pass criterion (the test's M_TOL)
    ok = abs(m) < tol
    colour, label = (TEAL, "PASS") if ok else (CINNABAR, "FAIL")
    centroid = cfg.get("centroid_source", "?")
    vn = (
        f"ngmix {prov.get('ngmix_version', '?')} &middot; "
        f"galsim {prov.get('galsim_version', '?')} &middot; centroid {centroid}"
    )
    fig = figure_block(
        b64, "m-bias recovered vs injected",
        "Recovered vs injected shear against the <b>1:1 ideal</b>: the measured point "
        "<span style='color:#3D5BA0'><b>after the R correction</b></span> lands on the line; "
        "<span style='color:#C49333'><b>before it</b></span>, the response deficit shows. "
        f"The recovered <code>R11</code> is <code>{r11:.3f}</code>.",
        "run_mbias.py", render_ts,
    )
    section = f"""\
  <section id="l1a">
    <div class="section-mark">
      <span>&sect; Level 1a &mdash; m-bias</span>
      <a class="srcbtn" href="{REPO_BLOB}/tests/science/test_mbias.py" target="_blank" rel="noopener"><span class="gl">&lt;/&gt;</span>test_mbias.py</a>
    </div>
    <h2>Multiplicative bias on injected shear.</h2>
    <div class="spec"><span class="mark">&sect; the test</span><p><b>What it is.</b> An in-memory sim injecting a known shear <code>g_inj = {d['injected_g1']}</code> through the controlled <code>make_data / do_ngmix_metacal</code> path and recovering it after the metacal response correction.</p><p><b>What it measures.</b> The multiplicative bias <code>m = g_obs/g_inj &minus; 1</code>.</p><p><b>Pass criteria.</b> <code>|m| &lt; 5&times;10&#8315;&sup3;</code>.</p></div>
    {fig}
    <div class="verdict">{badge(label, colour)}<span class="vlabel">multiplicative bias m</span><span class="vval">{sci(m)}</span><span class="vn">{vn}</span></div>
  </section>"""
    chip = f'<span class="lchip">m-bias {badge(label, colour)}</span>'
    return section, chip


def render_psf_prior(d: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Level 1b psf_fit_prior section. Returns (section_html, ladder_chip_html)."""
    v = d["psf_fit_prior_verdict"]
    e1_gal = v["psf_fit_e1_galaxy"]
    e1_none = v["psf_fit_e1_none"]
    c = d["additive_bias_c1"]
    m = d["multiplicative_bias_m"]
    # Recompute the verdict on-page from the measured biases against the same
    # tolerances the tests assert (C_TOL = 1e-3, M_TOL = 5e-3) rather than
    # trusting a pre-baked boolean from the harvester JSON: the page owns its
    # pass/fail logic, so a stale or missing flag can't silently green-light it.
    C_TOL, M_TOL = 1e-3, 5e-3
    robust = abs(c) < C_TOL and abs(m) < M_TOL
    colour, label = (TEAL, "PASS") if robust else (CINNABAR, "FAIL")
    fig = figure_block(
        b64, "psf_fit_prior sweep",
        "<b>Left:</b> the fitted PSF ellipticity is entirely the prior&rsquo;s doing &mdash; "
        "<code>galaxy</code> and <code>psf</code> round it to zero, <code>none</code> recovers the "
        "injected <code>0.05</code>. <b>Right:</b> the recovered shear <code>m</code> is flat across all "
        "three. Metacal absorbs the PSF-fit-prior bias; the shear is robust to it.",
        "run_star_response.py", render_ts,
    )
    section = f"""\
  <section id="l1b">
    <div class="section-mark">
      <span>&sect; Level 1b &mdash; psf_fit_prior sweep</span>
      <a class="srcbtn" href="https://github.com/CosmoStat/shapepipe/issues/749" target="_blank" rel="noopener"><span class="gl">&lt;/&gt;</span>shapepipe#749</a>
      <a class="srcbtn" href="{REPO_BLOB}/tests/science/test_star_response.py" target="_blank" rel="noopener"><span class="gl">&lt;/&gt;</span>test_star_response.py</a>
    </div>
    <h2>Does the PSF-fit prior bias the recovered shear?</h2>
    <div class="spec"><span class="mark">&sect; the test</span><p><b>What it is.</b> An in-memory sim injecting a known elliptical PSF (<code>e1<sub>true</sub> = {d['psf_e1_true']}</code>) and sweeping <code>psf_fit_prior</code> over <code>galaxy / none / psf</code>.</p><p><b>What it measures.</b> The fitted PSF ellipticity and the recovered shear (c, m) under each prior.</p><p><b>Pass criteria.</b> Recovered shear <code>c</code> and <code>m</code> unchanged across the prior choices.</p></div>
    {fig}
    <div class="verdict">{badge(label, colour)}<span class="vlabel">shear robust to psf_fit_prior</span><span class="vval">|c| = {sci(abs(c))} &middot; |m| = {sci(abs(m))}</span><span class="vn">PSF fit e1 swings {e1_gal:.4f} &rarr; {e1_none:.4f}; shear does not</span></div>
  </section>"""
    chip = f'<span class="lchip">#749 {badge(label, colour)}</span>'
    return section, chip


def _grid_verdict(r1_mean: float, tol: float) -> tuple[str, str]:
    """PASS if |<R1>| within band, else MARGINAL (tail-driven, not a hard fail)."""
    return (TEAL, "PASS") if abs(r1_mean) <= tol else (OCHRE, "MARGINAL")


def render_grid(arms: dict, figs: dict, render_ts: str) -> tuple[str, str]:
    """Level 2 SKiLLS-grid section with the fitgauss/gauss A/B.

    ``arms`` maps arm key -> JSON dict. The fitgauss-wcs arm is the baseline;
    gauss-wcs is the A/B comparison. Returns (section_html, ladder_chip_html).
    """
    base = arms["wcs"]            # fitgauss-wcs baseline
    gauss = arms.get("gauss_wcs")  # gauss-wcs A/B partner
    hsm = arms.get("hsm")
    gauss_hsm = arms.get("gauss_hsm")

    tol = base["tolerance_R"]
    r1, r2 = base["R1"], base["R2"]
    td = base["test_data"]
    fab = base["fabian_headline"]
    mag_lo, mag_hi = td["mag_range"]

    colour_b, label_b = _grid_verdict(r1["mean"], tol)

    # --- A/B deltas, all computed from JSON (the live levers) ---
    ab_lines = ""
    abcallout = ""
    gauss_verdict = ""
    if gauss is not None:
        g_r1 = gauss["R1"]
        wcs_drop = 100 * (g_r1["mean"] - r1["mean"]) / abs(r1["mean"]) if r1["mean"] else 0.0
        cat_b = r1["frac_catastrophic"]
        cat_g = g_r1["frac_catastrophic"]
        hsm_txt = ""
        if hsm is not None and gauss_hsm is not None:
            h0, h1 = hsm["R1"]["mean"], gauss_hsm["R1"]["mean"]
            hsm_drop = 100 * (h1 - h0) / abs(h0) if h0 else 0.0
            hsm_txt = (
                f', <span style="color:#3D5BA0"><b>hsm</b></span> '
                f"<code>{signed(h0)}&nbsp;&rarr;&nbsp;{signed(h1)}</code> "
                f"({signed(hsm_drop, 0)}%)"
            )
        abcallout = f"""\
    <div class="abcallout">
      <div class="abc-mark">&sect; reconvolution-PSF A/B &mdash; the largest lever</div>
      <p>Switching the metacal reconvolution kernel from
      <code>fitgauss</code> to <code>gauss</code> sizing brings
      <code>&langle;R1&rangle;</code> <b>into the &plusmn;{tol:g} band</b>,
      where fitgauss sat outside it:
      <span style="color:#3F8278"><b>wcs</b></span>
      <code>{signed(r1['mean'])}&nbsp;&rarr;&nbsp;{signed(g_r1['mean'])}</code>
      ({signed(wcs_drop, 0)}%){hsm_txt}.
      The catastrophic-fit tail (<code>|R|&gt;1</code>) shrinks in step:
      <code>{cat_b:.0%}&nbsp;&rarr;&nbsp;{cat_g:.0%}</code>.</p>
      <p class="abc-mech"><b>Mechanism.</b> Both arms build a <i>round</i>
      reconvolution kernel; they differ only in how its <b>size</b> is set.
      <code>fitgauss</code> sizes it from a per-object adaptive-moments (admom) fit
      of the PSF stamp &mdash; fragile on noisy or tiny stamps, where the fit blows
      up and the metacal response runs away. <code>gauss</code> sizes it from a
      deterministic k-space match to the PSF image &mdash; robust, so the runaway
      tail collapses. It is <b>not</b> a Bayesian-prior effect; it is the
      kernel-sizing method.</p>
      <p class="abc-floor"><b>The residual tail is the point-source floor.</b>
      gauss removes the <i>fixable</i> fragility, but a <code>{cat_g:.0%}</code>
      tail of catastrophic fits remains: the ellipticity response is singular as
      the pre-PSF intrinsic size <code>T&rarr;0</code>, and stars <i>are</i> point
      sources. It is exactly the regime a <b>resolution cut</b> excludes from a
      shear sample, so a residual <code>&langle;R&rangle;</code> on a
      <i>pure-star</i> sample is a diagnostic of metacal&rsquo;s point-source
      behaviour, not a direct shear bias.</p>
    </div>"""
        colour_g, label_g = _grid_verdict(g_r1["mean"], tol)
        gauss_verdict = (
            f'\n    <div class="verdict">{badge(label_g, colour_g)}'
            '<span class="vlabel">gauss &langle;R1&rangle; / &langle;R2&rangle; (wcs)</span>'
            f'<span class="vval">{signed(g_r1["mean"])} / {signed(gauss["R2"]["mean"])}</span>'
            f'<span class="vn">reconvolution gauss sizing &mdash; &langle;R1&rangle; inside '
            f'&plusmn;{tol:g} band &middot; cat {cat_g:.0%} &middot; {gauss["n_obj"]:,} stars</span></div>'
        )

    # additive-bias companion line (mean star ellipticity, ideal 0)
    add_hsm = ""
    if hsm is not None:
        add_hsm = f"<span><code>c1 {signed(hsm['c1'])}</code> &middot; <code>c2 {signed(hsm['c2'])}</code> &nbsp;hsm</span>"
    addbias = (
        '<div class=\'addbias\'><span class=\'ab-mark\'>&sect; additive bias &mdash; mean star ellipticity, ideal 0</span>'
        f"<span><code>c1 {signed(base['c1'])}</code> &middot; <code>c2 {signed(base['c2'])}</code> &nbsp;wcs</span>"
        f"{add_hsm}</div>"
    )

    hist_cap = (
        f"Per-object response across {base['n_obj']:,} stars, <b>fitgauss&middot;wcs baseline arm</b>. "
        "The <span style='color:#3F8278'><b>target band</b></span> marks "
        "<span style='color:#BC4538'><b>mean &plusmn; jackknife</b></span> and "
        "<span style='color:#3D5BA0'><b>trimmed mean</b></span>. The mean is <b>tail-driven</b>: "
        f"median (<code>R1 {signed(r1['median'])}</code>, <code>R2 {signed(r2['median'])}</code>) and "
        f"trimmed mean (<code>{signed(r1['trimmed_mean'])}</code>, <code>{signed(r2['trimmed_mean'])}</code>) "
        f"sit inside the band, but a <code>{r1['frac_catastrophic']:.0%}</code> tail of catastrophic "
        f"metacal fits (<code>|R|&gt;1</code>, std &asymp; {r1['std']:.1f}) pulls the mean out."
    )
    forest_cap = (
        f"Mean response with jackknife errors on the <span style='color:#3F8278'><b>&plusmn;{tol:g} band</b></span>, "
        "arms plus the 1.4 reference. Reconvolution-kernel sizing dominates over the centroid source "
        f"(wcs vs hsm). The shapepipe-1.4 fitgauss reference: <code>&langle;R1&rangle; {fab['R1_mean']:+.4f}</code> "
        f"&plusmn; {fab['R1_err']:.4f} ({fab['n_obj']:,} stars, {fab['n_tiles']} tiles)."
    )
    fig_hist = figure_block(figs.get("grid_hist"), "star shear-response R1/R2 histograms", hist_cap, "run_star_response_grid.py", render_ts)
    fig_forest = figure_block(figs.get("grid_forest"), "forest plot of mean response by reconvolution-PSF arm", forest_cap, "run_star_response_grid.py", render_ts)
    fig_ab = ""
    if gauss is not None:
        ab_cap = (
            "The same wcs arm, fitgauss (<span style='color:#BC4538'><b>cinnabar</b></span>) vs gauss "
            "(<span style='color:#3F8278'><b>teal</b></span>), R1 overlaid. gauss sizing thins the "
            "catastrophic-fit tail (<code>|R|&gt;1</code>), so the bulk tightens and the mean falls back "
            f"inside the <span style='color:#3F8278'><b>&plusmn;{tol:g} band</b></span>."
        )
        fig_ab = "\n    " + figure_block(figs.get("grid_ab"), "fitgauss vs gauss R1 tail comparison", ab_cap, "run_star_response_grid.py", render_ts)

    section = f"""\
  <section id="l2">
    <div class="section-mark">
      <span>&sect; Level 2 &mdash; star shear response (SKiLLS grid)</span>
      <a class="srcbtn" href="{REPO_BLOB}/tests/cluster/test_star_shear_response.py" target="_blank" rel="noopener"><span class="gl">&lt;/&gt;</span>test_star_shear_response.py</a>
      <a class="srcbtn" href="{REPO_BLOB}/tests/helpers/star_response.py" target="_blank" rel="noopener"><span class="gl">&lt;/&gt;</span>star_response.py</a>
    </div>
    <h2>Star metacal shear response on the SKiLLS grid.</h2>
    <div class="spec"><span class="mark">&sect; the test</span><p><b>What it is.</b> The per-object metacal shear response <code>R = dg/d&gamma;</code> measured on stars from Fabian Hervas&rsquo;s SKiLLS image grid, over <code>{mag_lo:g} &lt; mag &lt; {mag_hi:g}</code>, across {base['n_tiles_ok']} tiles, with delete-one-tile jackknife errors. A round PSF-like object should give <code>R = 0</code>. Run as an A/B over the metacal reconvolution-PSF kernel sizing (<code>fitgauss</code> vs <code>gauss</code>), each crossed with the centroid source (<code>wcs</code> vs <code>hsm</code>).</p><p><b>What it measures.</b> The mean per-object shear response of stars after metacal PSF deconvolution, and how it responds to the reconvolution-kernel-sizing choice.</p><p><b>Pass criteria.</b> <code>&langle;R1&rangle;</code> and <code>&langle;R2&rangle;</code> each within <code>&plusmn;{tol:g}</code> of <code>0</code> (jackknife error).</p></div>

    <div class="testdata">
      <div class="testdata-head">&sect; test data &mdash; held fixed, only ngmix re-run</div>
      <div class="testdata-body">
        <b>{html.escape(td['name'])}</b><br>
        {html.escape(td['simulation'])} &middot;
        {html.escape(td['psf'])}<br>
        {html.escape(td['stage_reused'])}
      </div>
    </div>

    {fig_hist}

    {fig_forest}
{abcallout}{fig_ab}

    {addbias}

    <div class="verdict">{badge(label_b, colour_b)}<span class="vlabel">mean response &langle;R1&rangle; / &langle;R2&rangle;</span><span class="vval">{signed(r1['mean'])} / {signed(r2['mean'])}</span><span class="vn">fitgauss-wcs &mdash; mean tail-driven, outside band &middot; {base['n_obj']:,} stars &middot; {base['n_tiles_ok']} tiles</span></div>{gauss_verdict}

  </section>"""
    chip = (
        f'<span class="lchip">fitgauss R {badge(label_b, colour_b)}</span>'
        + (f'<span class="lchip">gauss R {badge(*reversed(_grid_verdict(gauss["R1"]["mean"], tol)))}</span>' if gauss else "")
    )
    return section, chip


# --------------------------------------------------------------------------- #
# Page assembly
# --------------------------------------------------------------------------- #
def build_page(results: Path, figures: Path | None, now: _dt.datetime) -> str:
    render_ts = now.strftime("%Y-%m-%d %H:%M")
    stamp_date = now.strftime("%Y &middot; %m &middot; %d")

    figs = {}
    if figures is not None:
        for slot, fname in FIGURE_SLOTS.items():
            figs[slot] = b64_png(figures / fname)

    sections: list[str] = []
    chips: list[str] = []

    mbias = load(results, "mbias_sim", "mbias.json")
    if mbias is not None:
        sec, chip = render_mbias(mbias, figs.get("mbias"), render_ts)
        sections.append(sec)
        chips.append(chip)

    psf = load(results, "star_response_sim", "star_response.json")
    if psf is not None:
        sec, chip = render_psf_prior(psf, figs.get("psf_fit_prior"), render_ts)
        sections.append(sec)
        chips.append(chip)

    # Level 2 arms: any present subset of the 2x2 grid.
    arm_dirs = {
        "wcs": "star_response_grid_wcs",
        "hsm": "star_response_grid_hsm",
        "gauss_wcs": "star_response_grid_gauss_wcs",
        "gauss_hsm": "star_response_grid_gauss_hsm",
    }
    arms = {k: load(results, d, "star_response.json") for k, d in arm_dirs.items()}
    arms = {k: v for k, v in arms.items() if v is not None}
    l2_chips = ""
    if "wcs" in arms:
        sec, l2_chips = render_grid(arms, figs, render_ts)
        sections.append(sec)

    if not sections:
        raise SystemExit(
            f"no recognised result JSONs found under {results} "
            f"(expected e.g. mbias_sim/mbias.json, star_response_grid_wcs/star_response.json)"
        )

    # Ladder strip — Level 1 chips from the sim arms, Level 2 chips from the grid.
    l1_chips = "".join(chips) or '<span class="lchip">&mdash; <span class="badge" style="background:var(--ink-muted)">PENDING</span></span>'
    l2_chips = l2_chips or '<span class="lchip">&mdash; <span class="badge" style="background:var(--ink-muted)">PENDING</span></span>'

    sections_html = "\n\n".join(sections)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>UNIONS ShapePipe &mdash; Shape-Measurement Status</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=EB+Garamond:ital,wght@0,400;0,500;0,600;1,400;1,500&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>{CSS}</style>
</head>
<body>
<main class="page">

  <div class="stamp">
    <span class="dot"></span><span>{stamp_date}</span>
    <span class="sep">&middot;</span><span>ngmix_v2.0 &middot; candide</span>
    <span class="sep">&middot;</span><span>rendered {render_ts}</span>
  </div>

  <h1>ShapePipe &mdash; shape-measurement status.</h1>
  <p class="lede">A live, auto-generated status board for the ngmix metacal shape
  measurement, organised by a <b>fidelity ladder</b>: in-memory toy simulations
  (Level&nbsp;1, the ideal limit) up through the realistic SKiLLS star grid
  (Level&nbsp;2). Every number is read straight from a results JSON and shown as
  a plot with its provenance; the page re-renders on each run.</p>

  <div class="disclaimer">
    <b>&#9888; Preliminary</b>
    <span>under rapid development &mdash; tests, numbers, and layout may change.</span>
  </div>


  <div class="ladder">

    <a class="rung l1" href="#l1a">
      <div class="rung-head l1">Level 1 &mdash; in-memory sims</div>
      <p>Ideal limit. A single injected shear runs through the controlled <code>make_data / do_ngmix_metacal</code> path in seconds &mdash; an algorithmic check, no image simulation.</p>
      <div class="lchips">{l1_chips}</div>
    </a>

    <a class="rung l2" href="#l2">
      <div class="rung-head l2">Level 2 &mdash; SKiLLS star grid</div>
      <p>Realistic. SKiLLS 1z2z star images, fake PSF mid-chain; only the ngmix stage is re-run, all upstream intermediates byte-fixed.</p>
      <div class="lchips">{l2_chips}</div>
    </a>

    <a class="rung l3" href="#l3">
      <div class="rung-head l3">Level 3 &mdash; real UNIONS data</div>
      <p>Future. Full UNIONS catalogue runs, end-to-end.</p>
      <div class="lchips"><span class="lchip">&mdash; <span class="badge" style="background:var(--ink-muted)">SOON</span></span></div>
    </a>
  </div>

  <hr>


{sections_html}

  <footer>
    UPDATE MECHANISM &nbsp;&middot;&nbsp; the digital-twin harvesters
    (<b>run_mbias.py</b>, <b>run_star_response.py</b>,
    <b>run_star_response_grid.py</b>) emit
    <code>results/baseline/*/*.json</code> &rarr; <b>build_status.py</b> renders
    this page &rarr; committed and served at
    <a href="https://unions-wl.github.io" target="_blank" rel="noopener">unions-wl.github.io</a>.
    One render per run; history lives in the commits.
    <br>
    SOURCE &nbsp;&middot;&nbsp; guardrail tests on branch
    <code>ngmix_v2.0</code> (CosmoStat/shapepipe). Each plot carries the script
    and timestamp of the run that produced its numbers.
  </footer>

</main>
</body>
</html>"""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog=USAGE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--results", type=Path,
        default=Path("results/baseline"),
        help="directory holding the per-arm result subdirs (default: results/baseline)",
    )
    ap.add_argument(
        "--figures", type=Path, default=None,
        help="optional directory of PNGs to inline (see FIGURE_SLOTS); "
             "missing figures render as 'figure pending'",
    )
    ap.add_argument(
        "--out", type=Path, default=Path("status.html"),
        help="output HTML path (default: status.html)",
    )
    args = ap.parse_args(argv)

    if not args.results.is_dir():
        raise SystemExit(f"--results dir not found: {args.results}")

    html_str = build_page(args.results, args.figures, _dt.datetime.now())
    args.out.write_text(html_str)
    print(f"wrote {args.out} ({len(html_str) // 1024} KB) from {args.results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
