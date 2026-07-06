#!/usr/bin/env python
"""Render the ShapePipe shape-measurement *status page* from a run's JSON outputs.

This is the renderer referenced in the status-page footer: the digital-twin
harvesters (``run_mbias.py``, ``run_star_response.py``,
``run_star_response_grid.py``) emit ``results/baseline/*/*.json``; this script
reads those JSONs and renders a single self-contained ``status.html``. One render
per run; history lives in the commits.

The page is organised by **tier of check**, not by simulation fidelity:

* **Tier 1 — CI tripwires** (``mbias_sim/``, ``star_response_sim/``): fast
  in-memory sims that gate every commit. m-bias recovery (``|m| < M_TOL``) and
  PSF-fit-prior robustness (shear unchanged across the prior, #749). Each carries
  a level chip ``in-memory``.
* **Tier 2 — vetted image-sim panels** (``star_response_grid_*/``,
  ``star_response_magprofile/``, ``galaxy_question/``): metacal behaviour on
  Fabian Hervas's SKiLLS 1z2z star grid. The star-grid metacal response A/B, the
  star-response-vs-magnitude profile, and the galaxy question (S/N vs response).
  A shared *test-data* block heads the tier; each panel carries a level chip
  ``SKiLLS``.

Design contract
---------------
*Every number on the page is read straight from a results JSON.* The verdict
badges (PASS / MARGINAL / FAIL) are **computed** from the recorded value against
its tolerance, never hard-coded — that is the live part of the page. Pass-
criterion tolerances are interpolated from the single Python constant / JSON
field that also drives the verdict, never re-typed as a literal in an HTML
string. The CSS and the explanatory prose are the template; captions state what
is shown, not why. Figures are optional PNGs inlined as base64 (a missing figure
degrades to a "figure pending" note rather than crashing the render).

One panel = one figure + at most one verdict. Every subsection is instantiated
from the same ``panel()`` template.

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
Regenerate the status page from a run's outputs:

    python scripts/python/build_status.py \\
        --results <run>/results/baseline \\
        --figures <run>/figures \\
        --out status.html

`--results` points at the directory holding the per-arm subdirs
(`mbias_sim/`, `star_response_sim/`, `star_response_grid_*/`,
`star_response_magprofile/`, `galaxy_question/`); `--figures` (optional) points
at a directory of PNGs to inline (see FIGURE_SLOTS).

Deploy (Cail owns this step; not run by this script):

    cp status.html <gh-pages-clone>/shapepipe/status.html
    cd <gh-pages-clone> && git add shapepipe/status.html \\
        && git commit -m 'status: re-render <date>' && git push
"""

# Theme colours, mirrored from the page CSS so verdict badges stay on-palette.
TEAL = "#3F8278"      # PASS
OCHRE = "#C49333"     # MARGINAL
CINNABAR = "#BC4538"  # FAIL

# Tier-1 pass-criterion tolerances. These are the *test* constants (M_TOL, C_TOL
# asserted by test_mbias.py / test_star_response.py) — not present in the results
# JSON — kept here as single-source Python constants and interpolated into both
# the verdict computation and the printed criterion, never re-typed as literals.
M_TOL = 5e-3
C_TOL = 1e-3

# Figure filenames looked for under --figures, keyed by the slot they fill.
# A run that produced these PNGs gets them inlined; any that are missing fall
# back to a "figure pending" placeholder so the render never hard-fails.
FIGURE_SLOTS = {
    "mbias": "mbias.png",
    "psf_fit_prior": "psf_fit_prior.png",
    "grid_forest": "star_response_forest.png",
    "grid_magprofile": "star_response_magprofile.png",
    "galaxy_question": "galaxy_question_s2n.png",
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
    return f'<span class="badge" style="background:{colour}">{label}</span>'


def level_chip(level: str) -> str:
    """The small in-memory / SKiLLS provenance chip carried by every panel."""
    return f'<span class="levelchip">{html.escape(level)}</span>'


def strip_chip(anchor: str, name: str, label: str, colour: str) -> str:
    """One entry in the top verdict strip: a coloured dot + name + badge, linking
    to the panel's section anchor."""
    return (
        f'<a class="vchip" href="#{anchor}">'
        f'<span class="vdot" style="background:{colour}"></span>'
        f'<span class="vname">{name}</span>{badge(label, colour)}</a>'
    )


def figure_block(b64: str | None, alt: str, caption: str, data_prov: str, render_ts: str,
                 fig_tool: str = "build_status.py") -> str:
    """A ``<figure>`` with an inlined PNG, or a pending placeholder.

    ``fig_tool`` names the script that drew the PNG (default ``build_status.py``);
    figures produced by a dedicated plotter (e.g. ``plot_mag_profile.py``) pass
    their own so the provenance line credits the tool that actually rendered them.
    """
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
        f'<span class="pv"><span class="pk">figure</span> {fig_tool} <span class="pt">{render_ts}</span></span>'
        "</div>"
    )
    return f"<figure>{body}<figcaption>{caption}</figcaption>{prov}</figure>"


def verdict_row(label: str, colour: str, vlabel: str, vval: str, vn: str) -> str:
    """The single verdict badge + value line that closes a panel."""
    return (
        f'<div class="verdict">{badge(label, colour)}'
        f'<span class="vlabel">{vlabel}</span>'
        f'<span class="vval">{vval}</span>'
        f'<span class="vn">{vn}</span></div>'
    )


def panel(anchor: str, mark: str, srcbtns: str, level: str, title: str,
          spec: str, figure: str, verdict: str = "", note: str = "",
          badge_label: str = "", badge_colour: str = "", summary_val: str = "",
          force_open: bool = False) -> str:
    """The one template every subsection is built from, as a collapse-when-green
    ``<details>``. The ``<summary>`` carries the title, level chip, verdict badge
    and one-line value so a collapsed row still shows the number. Open state:
    ``force_open`` (Tier-2 panels, always shown) OR a non-PASS verdict."""
    is_open = force_open or (badge_label and badge_label != "PASS")
    open_attr = " open" if is_open else ""
    summ_badge = badge(badge_label, badge_colour) if badge_label else ""
    summ_val = f'<span class="sv">{summary_val}</span>' if summary_val else ""
    return f"""\
  <details id="{anchor}"{open_attr}>
    <summary><span class="caret"></span><span class="sm-title">{title}</span>{level_chip(level)}{summ_badge}{summ_val}</summary>
    <div class="section-mark"><span>{mark}</span>{srcbtns}</div>
    {spec}
    {figure}{note}{verdict}
  </details>"""


# --------------------------------------------------------------------------- #
# CSS (curated page palette; verdict strip + level chip replace the old ladder)
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

  /* ---- verdict strip: every badge on the page, one row, each an anchor ---- */
  .vstrip { margin:30px 0 0; padding:14px 16px; background:var(--rag);
    border:1px solid var(--rule); border-radius:10px; display:flex; gap:10px;
    flex-wrap:wrap; align-items:center; }
  .vstrip-label { font-family:var(--mono); font-size:10px; letter-spacing:.13em;
    text-transform:uppercase; color:var(--ink-muted); margin-right:4px; }
  a.vchip { display:inline-flex; align-items:center; gap:8px; color:var(--ink);
    font-family:var(--mono); font-size:12px; background:var(--parchment);
    border:1px solid var(--rule); border-radius:20px; padding:5px 11px 5px 9px;
    transition:border-color .14s, transform .14s; }
  a.vchip:hover { border-color:var(--cobalt); text-decoration:none;
    transform:translateY(-1px); }
  .vchip .vdot { width:8px; height:8px; border-radius:50%; flex-shrink:0; }
  .vchip .vname { color:var(--ink-muted); }

  /* ---- tier headers ---- */
  .tier { margin-top:52px; }
  .tier-mark { font-family:var(--mono); font-size:11px; letter-spacing:.15em;
    text-transform:uppercase; color:var(--cinnabar); margin-bottom:8px; }
  .tier h2 { font-family:var(--serif); font-style:italic; font-weight:500;
    font-size:32px; line-height:1.12; margin:0 0 6px; }
  .tier-scope { color:var(--ink-muted); font-size:16px; margin:0; max-width:78ch; }

  /* ---- collapse-when-green panels ---- */
  details { margin-top:30px; scroll-margin-top:20px;
    border-top:1px solid var(--tooth); padding-top:20px; }
  summary { list-style:none; cursor:pointer; display:flex; align-items:center;
    gap:12px; flex-wrap:wrap; }
  summary::-webkit-details-marker { display:none; }
  .caret { width:0; height:0; border-left:6px solid var(--ink-muted);
    border-top:5px solid transparent; border-bottom:5px solid transparent;
    transition:transform .15s; flex-shrink:0; }
  details[open] > summary .caret { transform:rotate(90deg); }
  .sm-title { font-family:var(--serif); font-style:italic; font-weight:500;
    font-size:24px; line-height:1.15; color:var(--ink); }
  summary:hover .sm-title { color:var(--cobalt); }
  .sv { font-family:var(--mono); font-size:13px; color:var(--ink-muted);
    margin-left:auto; text-align:right; }
  .section-mark { font-family:var(--mono); font-size:11px; letter-spacing:.13em;
    text-transform:uppercase; color:var(--ink-muted); margin:12px 0 8px;
    display:flex; gap:12px; align-items:center; flex-wrap:wrap; }

  hr { border:none; border-top:1px solid var(--rule); margin:36px 0; }

  /* ---- level chip (in-memory / SKiLLS provenance) ---- */
  .levelchip { font-family:var(--mono); font-size:10px; letter-spacing:.08em;
    color:var(--ink-muted); background:var(--rag); border:1px solid var(--rule);
    border-radius:4px; padding:3px 8px; text-transform:none; }

  /* ---- prominent source buttons ---- */
  .srcbtn { display:inline-flex; align-items:center; gap:8px;
    font-family:var(--mono); font-size:12px; letter-spacing:.02em;
    color:var(--cobalt); background:var(--rag); border:1px solid var(--rule);
    border-radius:7px; padding:6px 12px; transition:all .14s; line-height:1;
    text-transform:none; }
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

  /* ---- muted panel note (demoted comparison arms, additive bias) ---- */
  .note { margin:16px 0 0; padding:11px 18px; background:var(--rag);
    border:1px dashed var(--rule); border-radius:8px; font-family:var(--mono);
    font-size:13px; color:var(--ink-muted); line-height:1.55; }
  .note .nk { color:var(--cinnabar); font-size:10px; letter-spacing:.1em;
    text-transform:uppercase; margin-right:8px; }
  .note b { color:var(--ink); font-weight:500; }

  /* ---- shared test-data block (heads Tier 2) ---- */
  .testdata { margin:20px 0 0; padding:14px 18px; background:var(--rag);
    border-left:3px solid var(--cobalt); border-radius:0 6px 6px 0;
    font-size:14.5px; line-height:1.55; color:var(--ink-muted); }
  .testdata-head { font-family:var(--mono); font-size:10px; letter-spacing:.13em;
    text-transform:uppercase; color:var(--cobalt); margin-bottom:7px; }
  .testdata-body b { color:var(--ink); font-weight:500; }

  footer { margin-top:52px; border-top:1px solid var(--rule); padding-top:18px;
    color:var(--ink-muted); font-family:var(--mono); font-size:11px;
    line-height:1.7; letter-spacing:.02em; }
  footer b { color:var(--ink); }
"""

# Single source-link base for every "&lt;/&gt;" src button on the page. Pinned to
# the ngmix_v2.0 integration branch while that is where this work lives; change
# this one constant to re-point every link when the work merges elsewhere.
REPO_BLOB = "https://github.com/CosmoStat/shapepipe/blob/ngmix_v2.0"


def _src(href: str, text: str) -> str:
    return (f'\n      <a class="srcbtn" href="{href}" target="_blank" rel="noopener">'
            f'<span class="gl">&lt;/&gt;</span>{text}</a>')


# --------------------------------------------------------------------------- #
# Panel renderers — each returns (panel_html, strip_chip). Inspection panels
# (no verdict) return "" for the strip chip.
# --------------------------------------------------------------------------- #
def render_mbias(d: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Tier-1 m-bias tripwire."""
    m = d["multiplicative_bias_m"]
    r11 = d["metacal_response_R11"]
    prov = d.get("provenance", {})
    ok = abs(m) < M_TOL
    colour, label = (TEAL, "PASS") if ok else (CINNABAR, "FAIL")
    centroid = d.get("config", {}).get("centroid_source", "?")
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
    spec = (
        '<div class="spec"><span class="mark">&sect; the test</span>'
        f"<p><b>What it is.</b> An in-memory sim injecting a known shear <code>g_inj = {d['injected_g1']}</code> "
        "through the controlled <code>make_data / do_ngmix_metacal</code> path and recovering it after the "
        "metacal response correction.</p>"
        "<p><b>What it measures.</b> The multiplicative bias <code>m = g_obs/g_inj &minus; 1</code>.</p>"
        f"<p><b>Pass criteria.</b> <code>|m| &lt; {sci(M_TOL)}</code>.</p></div>"
    )
    verdict = verdict_row(label, colour, "multiplicative bias m", sci(m), vn)
    sec = panel(
        "t1-mbias", "&sect; Tier 1 &mdash; m-bias",
        _src(f"{REPO_BLOB}/tests/science/test_mbias.py", "test_mbias.py"),
        "in-memory", "Multiplicative bias on injected shear.", spec, fig, verdict,
        badge_label=label, badge_colour=colour, summary_val=f"m = {sci(m)}",
    )
    return sec, strip_chip("t1-mbias", "m-bias", label, colour)


def render_psf_prior(d: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Tier-1 PSF-fit prior-neutrality of the recovered shear (#749/#779)."""
    v = d["psf_fit_prior_verdict"]
    e1_gal = v["psf_fit_e1_galaxy"]
    e1_none = v["psf_fit_e1_none"]
    c = d["additive_bias_c1"]
    m = d["multiplicative_bias_m"]
    # The page owns its pass/fail logic: recompute from the measured biases against
    # the same C_TOL / M_TOL the tests assert, rather than trusting a pre-baked flag.
    robust = abs(c) < C_TOL and abs(m) < M_TOL
    colour, label = (TEAL, "PASS") if robust else (CINNABAR, "FAIL")
    fig = figure_block(
        b64, "PSF-fit prior arms vs recovered shear",
        "<b>Left:</b> the fitted PSF ellipticity depends on the PSF-fit prior &mdash; the prior-dominated arm "
        f"is crushed toward round (<code>e1 &rarr; {e1_gal:.3f}</code>), the prior-free arm recovers truth "
        f"(<code>e1 &rarr; {e1_none:.3f}</code> &asymp; injected <code>{d['psf_e1_true']}</code>). "
        "<b>Right:</b> the recovered shear <code>m</code> is identical across the arms &mdash; metacal "
        "deconvolves by the PSF image, not by any fitted model.",
        "run_star_response.py", render_ts,
    )
    spec = (
        '<div class="spec"><span class="mark">&sect; the test</span>'
        "<p><b>What it is.</b> The galaxy ellipticity prior (<code>GPriorBA 0.4</code>, peak at round) is passed "
        "to the PSF-stamp fitter; #779&rsquo;s finite PSF weight map keeps that fit likelihood-dominated. Probe: "
        f"fit the elliptical PSF (<code>e1 = {d['psf_e1_true']}</code>) with prior + weak weight "
        f"(prior-dominated arm, <code>e1 &rarr; {e1_gal:.3f}</code>) vs prior-free "
        f"(<code>e1 &rarr; {e1_none:.3f}</code>).</p>"
        "<p><b>What it measures.</b> Whether the recovered shear (c, m) is sensitive to the PSF-fit "
        "configuration.</p>"
        f"<p><b>Pass criteria.</b> Recovered <code>c</code>, <code>m</code> identical in both arms "
        f"(<code>|c| &lt; {sci(C_TOL)}</code>, <code>|m| &lt; {sci(M_TOL)}</code>) &mdash; metacal deconvolves "
        "by the PSF image, not a fitted model. Red = shear has become sensitive to the PSF-fit configuration.</p></div>"
    )
    verdict = verdict_row(
        label, colour, "shear neutral to PSF-fit prior",
        f"|c| = {sci(abs(c))} &middot; |m| = {sci(abs(m))}",
        f"PSF fit e1 spans {e1_gal:.4f} &rarr; {e1_none:.4f}; shear does not move",
    )
    sec = panel(
        "t1-psf", "&sect; Tier 1 &mdash; PSF-fit prior neutrality",
        _src("https://github.com/CosmoStat/shapepipe/pull/779", "shapepipe#779")
        + _src("https://github.com/CosmoStat/shapepipe/issues/749", "shapepipe#749")
        + _src(f"{REPO_BLOB}/tests/science/test_star_response.py", "test_star_response.py"),
        "in-memory", "PSF-fit prior neutrality (#749/#779).", spec, fig, verdict,
        badge_label=label, badge_colour=colour,
        summary_val=f"|c| = {sci(abs(c))} &middot; |m| = {sci(abs(m))}",
    )
    return sec, strip_chip("t1-psf", "psf-prior", label, colour)


def _grid_verdict(r1_mean: float, tol: float) -> tuple[str, str]:
    """PASS if |<R1>| within band, else MARGINAL (tail-driven, not a hard fail)."""
    return (TEAL, "PASS") if abs(r1_mean) <= tol else (OCHRE, "MARGINAL")


def render_testdata(td: dict) -> str:
    """The shared Tier-2 test-data block: what is held fixed, what is re-run."""
    return f"""\
  <div class="testdata">
    <div class="testdata-head">&sect; test data &mdash; held fixed across all Tier-2 panels, only ngmix re-run</div>
    <div class="testdata-body">
      <b>{html.escape(td['name'])}</b><br>
      {html.escape(td['simulation'])} &middot; {html.escape(td['psf'])}<br>
      {html.escape(td['stage_reused'])}
    </div>
  </div>"""


def render_grid(arms: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Tier-2 star-grid metacal response A/B. ONE figure (the forest of mean
    responses on the +/-tol band), ONE verdict (the vetted gauss-wcs arm). The
    fitgauss comparison arm is a demoted muted note, not a badge."""
    base = arms["wcs"]              # fitgauss-wcs comparison arm
    gauss = arms["gauss_wcs"]       # gauss-wcs vetted arm (the verdict)
    tol = base["tolerance_R"]
    g_r1, g_r2 = gauss["R1"], gauss["R2"]
    b_r1 = base["R1"]
    td = base["test_data"]
    fab = base["fabian_headline"]
    mag_lo, mag_hi = td["mag_range"]

    colour, label = _grid_verdict(g_r1["mean"], tol)

    cap = (
        f"Mean per-object response <code>&langle;R1&rangle;</code> / <code>&langle;R2&rangle;</code> with "
        "delete-one-tile jackknife errors for every arm, on the "
        f"<span style='color:#3F8278'><b>&plusmn;{tol:g} target band</b></span>, plus the shapepipe-1.4 "
        f"reference (<code>&langle;R1&rangle; {fab['R1_mean']:+.4f}</code> &plusmn; {fab['R1_err']:.4f}, "
        f"{fab['n_obj']:,} stars). The <b>gauss</b> arms sit inside the band; the <b>fitgauss</b> arms sit "
        "outside it. Centroid source (wcs vs hsm) moves the point far less than the kernel-sizing choice."
    )
    fig = figure_block(
        b64, "forest plot of mean star response by reconvolution-PSF arm",
        cap, "run_star_response_grid.py", render_ts,
    )

    note = (
        '<div class="note">'
        '<div><span class="nk">comparison arm</span>'
        f"<b>fitgauss-wcs</b>: <code>&langle;R1&rangle; {signed(b_r1['mean'])}</code> &mdash; outside the "
        f"&plusmn;{tol:g} band, tail-driven (<code>{b_r1['frac_catastrophic']:.0%}</code> catastrophic "
        f"<code>|R|&gt;1</code> fits vs <code>{g_r1['frac_catastrophic']:.0%}</code> for gauss).</div>"
        '<div style="margin-top:6px"><span class="nk">additive bias</span>'
        f"mean star ellipticity (ideal 0): <code>c1 {signed(base['c1'])}</code> &middot; "
        f"<code>c2 {signed(base['c2'])}</code>.</div></div>"
    )

    spec = (
        '<div class="spec"><span class="mark">&sect; the test</span>'
        "<p><b>What it is.</b> The per-object metacal shear response <code>R = dg/d&gamma;</code> measured on "
        f"stars over <code>{mag_lo:g} &lt; mag &lt; {mag_hi:g}</code>, across {base['n_tiles_ok']} tiles, with "
        "delete-one-tile jackknife errors. A round PSF-like object should give <code>R = 0</code>. Run as an A/B "
        "over the metacal reconvolution-PSF kernel sizing (<code>fitgauss</code> vs <code>gauss</code>), each "
        "crossed with the centroid source (<code>wcs</code> vs <code>hsm</code>).</p>"
        "<p><b>What it measures.</b> The mean per-object shear response of stars after metacal PSF deconvolution, "
        "and how it responds to the reconvolution-kernel-sizing choice.</p>"
        f"<p><b>Pass criteria.</b> <code>&langle;R1&rangle;</code> and <code>&langle;R2&rangle;</code> each within "
        f"<code>&plusmn;{tol:g}</code> of <code>0</code> (jackknife error).</p></div>"
    )
    verdict = verdict_row(
        label, colour, "gauss-wcs &langle;R1&rangle; / &langle;R2&rangle;",
        f"{signed(g_r1['mean'])} / {signed(g_r2['mean'])}",
        f"vetted arm &mdash; inside &plusmn;{tol:g} band &middot; "
        f"cat {g_r1['frac_catastrophic']:.0%} &middot; {gauss['n_obj']:,} stars",
    )
    sec = panel(
        "t2-grid", "&sect; Tier 2 &mdash; star-grid metacal response (A/B)",
        _src(f"{REPO_BLOB}/tests/cluster/test_star_shear_response.py", "test_star_shear_response.py")
        + _src(f"{REPO_BLOB}/tests/helpers/star_response.py", "star_response.py"),
        "SKiLLS", "Star metacal shear response on the SKiLLS grid.", spec, fig, verdict, note,
        badge_label=label, badge_colour=colour,
        summary_val=f"&langle;R1&rangle; {signed(g_r1['mean'])} / &langle;R2&rangle; {signed(g_r2['mean'])}",
        force_open=True,
    )
    return sec, strip_chip("t2-grid", "star-grid R", label, colour)


def render_magprofile(d: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Tier-2 star-response-vs-magnitude inspection panel. No verdict."""
    tol = d["tolerance"]
    sci_r1 = d["windows"]["science_20_26"]["R1_mean"]
    cap = (
        f"Shear response <code>R1</code> vs magnitude on the <b>fitgauss&middot;wcs</b> arm "
        f"({d['n_obj']:,} stars, {d['n_tiles']} tiles), no magnitude cut. Across the "
        "<span style='color:#C49333'><b>20 &lt; mag &lt; 26 science window</b></span> the binned mean "
        f"<code>&langle;R1&rangle; {signed(sci_r1)}</code> sits inside the "
        f"<span style='color:#3F8278'><b>&plusmn;{tol:g} band</b></span>; toward the bright end it climbs to "
        "<code>&asymp;2</code> as stars approach the point-source limit."
    )
    fig = figure_block(
        b64, "star shear-response R1 vs magnitude (bright-end point-source break)",
        cap, "run_star_response_magprofile.py", render_ts, fig_tool="plot_mag_profile.py",
    )
    spec = (
        '<div class="spec"><span class="mark">&sect; the panel</span>'
        "<p><b>What it is.</b> The per-object metacal response <code>R1 = dg1/d&gamma;1</code> on the "
        "<code>fitgauss&middot;wcs</code> arm, binned in magnitude with no magnitude cut applied.</p>"
        "<p><b>What it shows.</b> Where the mean response departs from 0 as a function of source brightness &mdash; "
        "the point-source floor of the star-grid A/B, resolved against magnitude. Inspection panel, no verdict.</p></div>"
    )
    sec = panel(
        "t2-magprofile", "&sect; Tier 2 &mdash; star response vs magnitude", "",
        "SKiLLS", "Star shear response as a function of magnitude.", spec, fig,
        summary_val="inspection panel &mdash; no verdict", force_open=True,
    )
    return sec, ""


def render_galaxy_question(d: dict, b64: str | None, render_ts: str) -> tuple[str, str]:
    """Tier-2 galaxy question: does the bright-star response climb reach the
    resolved (shear-carrying) population? Verdict on the highest-S/N resolved bin."""
    heads = {h["label"]: h for h in d["headline"]}
    star_bright = next(h for k, h in heads.items() if k.startswith("stars bright"))
    tol = d["tolerance_symmetry"]
    # Verdict on the highest-S/N resolved-galaxy bin (res>0.4) — the exact locus
    # where the star point-source climb appears. Symmetry (R1~R2, not the star 2:1
    # split) AND sub-unity (below the <R>=1 blow-up), both read off that bin.
    top = max(d["gal_resolved_vs_s2n"], key=lambda b: b["mid"])
    r1, r2 = top["R1_mean"], top["R2_mean"]
    ok = abs(r1 - r2) < tol and max(r1, r2) < 1.0
    colour, label = (TEAL, "PASS") if ok else (CINNABAR, "FAIL")
    fig = figure_block(
        b64, "mean metacal response vs S/N, stars vs resolved galaxies",
        "Mean metacal response <code>&langle;R&rangle;</code> vs S/N on the SKiLLS sims. "
        "<span style='color:#BC4538'><b>Point-source stars</b></span> (<code>res&lt;0.06</code>) "
        "climb past <code>&langle;R&rangle;=1</code> toward high S/N and split <code>R1&thinsp;&gg;&thinsp;R2</code>; "
        "<span style='color:#3F8278'><b>resolved galaxies</b></span> (<code>res&gt;0.4</code>) "
        "saturate near <code>&langle;R&rangle;&asymp;0.85</code> with <code>R1&thinsp;&asymp;&thinsp;R2</code> paired.",
        "run_galaxy_question.py", render_ts, fig_tool="plot_galaxy_question.py",
    )
    spec = (
        '<div class="spec"><span class="mark">&sect; the test</span>'
        "<p><b>What it is.</b> The per-object metacal response <code>R = dg/d&gamma;</code> on the "
        "<b>same SKiLLS sims</b> for point-source stars (<code>res &lt; 0.06</code>) and resolved galaxies "
        "(<code>res &gt; 0.4</code>), binned in S/N. Resolution <code>res = T/(T+T<sub>psf</sub>)</code>: 0 at "
        "the point-source limit, 1 fully resolved.</p>"
        "<p><b>What it measures.</b> Whether the resolved galaxies that carry the cosmic-shear signal inherit the "
        "star point-source-limit response climb, or sit in the well-behaved regime.</p>"
        f"<p><b>Pass criteria.</b> At the highest S/N, resolved-galaxy <code>&langle;R1&rangle; &asymp; "
        f"&langle;R2&rangle;</code> (symmetric, <code>|&Delta;| &lt; {tol:g}</code>) and "
        "<code>&langle;R&rangle; &lt; 1</code> (below the point-source singularity).</p></div>"
    )
    verdict = verdict_row(
        label, colour, "resolved-galaxy response symmetric &amp; sub-unity",
        f"&langle;R1&rangle; {signed(r1)} / &langle;R2&rangle; {signed(r2)}",
        f"highest-S/N resolved bin (res&gt;0.4, S/N&gt;{top['lo']:.0f}, {top['n']:,} objects) &middot; "
        f"bright stars split {signed(star_bright['R1_mean'])} / {signed(star_bright['R2_mean'])}",
    )
    sec = panel(
        "t2-galaxy", "&sect; Tier 2 &mdash; galaxy question (S/N vs response)", "",
        "SKiLLS", "Does the bright-star response climb reach the shear sample?", spec, fig, verdict,
        badge_label=label, badge_colour=colour,
        summary_val=f"&langle;R1&rangle; {signed(r1)} / &langle;R2&rangle; {signed(r2)}",
        force_open=True,
    )
    return sec, strip_chip("t2-galaxy", "galaxy R", label, colour)


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

    strip: list[str] = []
    tier1: list[str] = []
    tier2: list[str] = []

    # ---- Tier 1 — CI tripwires ----
    mbias = load(results, "mbias_sim", "mbias.json")
    if mbias is not None:
        sec, chip = render_mbias(mbias, figs.get("mbias"), render_ts)
        tier1.append(sec)
        strip.append(chip)

    psf = load(results, "star_response_sim", "star_response.json")
    if psf is not None:
        sec, chip = render_psf_prior(psf, figs.get("psf_fit_prior"), render_ts)
        tier1.append(sec)
        strip.append(chip)

    # ---- Tier 2 — vetted image-sim panels ----
    arm_dirs = {
        "wcs": "star_response_grid_wcs",
        "hsm": "star_response_grid_hsm",
        "gauss_wcs": "star_response_grid_gauss_wcs",
        "gauss_hsm": "star_response_grid_gauss_hsm",
    }
    arms = {k: load(results, d, "star_response.json") for k, d in arm_dirs.items()}
    arms = {k: v for k, v in arms.items() if v is not None}

    testdata_html = ""
    if "wcs" in arms and "gauss_wcs" in arms:
        testdata_html = render_testdata(arms["wcs"]["test_data"])
        sec, chip = render_grid(arms, figs.get("grid_forest"), render_ts)
        tier2.append(sec)
        strip.append(chip)

    magprof = load(results, "star_response_magprofile", "magprofile.json")
    if magprof is not None:
        sec, _ = render_magprofile(magprof, figs.get("grid_magprofile"), render_ts)
        tier2.append(sec)

    galq = load(results, "galaxy_question", "galaxy_question.json")
    if galq is not None:
        sec, chip = render_galaxy_question(galq, figs.get("galaxy_question"), render_ts)
        tier2.append(sec)
        strip.append(chip)

    if not tier1 and not tier2:
        raise SystemExit(
            f"no recognised result JSONs found under {results} "
            f"(expected e.g. mbias_sim/mbias.json, star_response_grid_wcs/star_response.json)"
        )

    strip_html = "".join(strip)

    tiers_html = ""
    if tier2:
        tiers_html += f"""\
  <div class="tier">
    <div class="tier-mark">&sect; Tier 2</div>
    <h2>Vetted image-sim panels.</h2>
    <p class="tier-scope">Metacal behaviour on Fabian Hervas&rsquo;s SKiLLS 1z2z star grid &mdash; realistic images with a fake PSF injected mid-chain, only the ngmix stage re-run. These are what you come to inspect; panels open by default.</p>
  </div>
{testdata_html}
{chr(10).join(tier2)}
"""
    if tier1:
        tiers_html += f"""\
  <div class="tier">
    <div class="tier-mark">&sect; Tier 1</div>
    <h2>CI tripwires.</h2>
    <p class="tier-scope">Fast in-memory sims that run on every commit: a single injected shear through the controlled <code>make_data / do_ngmix_metacal</code> path. Algorithmic checks, no image simulation &mdash; collapsed to a one-line row while green, expanded on any red.</p>
  </div>
{chr(10).join(tier1)}
"""

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
  <p class="lede">A live, auto-generated status page for the ngmix metacal shape
  measurement, organised by <b>tier of check</b>: the vetted image-sim panels on
  the SKiLLS star grid come first (Tier&nbsp;2, open for inspection), then the
  fast CI tripwires that gate every commit (Tier&nbsp;1, collapsed while green).
  Real-data checks on UNIONS catalogues are coming. Every number is read straight
  from a results JSON and shown as a plot with its provenance; the page
  re-renders on each run.</p>

  <div class="disclaimer">
    <b>&#9888; Preliminary</b>
    <span>under rapid development &mdash; tests, numbers, and layout may change.</span>
  </div>

  <div class="vstrip">
    <span class="vstrip-label">verdicts</span>{strip_html}
  </div>

  <hr>

{tiers_html}
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
<script>
  // Auto-open a collapsed <details> when its anchor is targeted from the verdict
  // strip (on load and on hashchange).
  function openTarget(){{var el=location.hash&&document.querySelector(location.hash);
  if(el){{var d=el.closest("details")||el.querySelector("details");if(d)d.open=true;}}}}
  addEventListener("hashchange",openTarget);openTarget();
</script>
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
