"""Reclamation must not change what ``sp report`` says about a unit.

``workflow/scripts/run_report.py`` reads a unit's verdict from two places with
the same shape and different lifetimes: the records on disk
(``manifests/*.json`` + ``logs/*.json``), and — once ``clean_tile`` or
``clean_exposure`` has deleted those — the copies absorbed into the unit's
``cleaned.json``. Both readers face the same collision, because several files
map to one ``(unit, stage)``: a stage's manifest and its byte-identical log, and
the eight ngmix chunks, which all carry ``stage: "tile_ngmix"``.

If the two readers resolve that collision differently, RECLAMATION SILENTLY
EDITS THE REPORT. It did: absorption used first-key-wins, so only
``tile_ngmix_1`` survived and a ``warn`` on any other chunk became ``complete``
the moment the tile was cleaned — and the tile dropped out of the "tiles not
complete" table. This module pins the invariant that catches that whole class:
**the tombstone path and the on-disk path agree, stage for stage.**

Deliberately container-free — run_report.py is stdlib-only, so nothing here
imports astropy, numpy, or shapepipe.
"""

import importlib.util
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "run_report.py"


def _load():
    """Import the script by path — ``workflow/scripts`` is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; bin/sp calls it by path"
    spec = importlib.util.spec_from_file_location("_run_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


run_report = _load()

TILE = "186.307"
NGMIX_CHUNKS = 8


def _manifest(stage, status="complete", found=1, expect=1):
    """One completeness verdict, in the shape completeness.py writes."""
    return {
        "stage": stage, "level": "tile", "unit": TILE, "status": status,
        "runners": {f"{stage}_runner": {
            "found": found, "expect": expect,
            "status": status, "warn": status == "warn"}},
        "failures": [],
    }


def _records(ngmix_status=None):
    """A finished tile's manifests, keyed by file stem as the store keys them.

    ``ngmix_status`` is ``{chunk_number: status}``; unlisted chunks are complete.
    """
    out = {s: _manifest(s) for s in run_report.TILE_STAGES if s != "tile_ngmix"}
    for k in range(1, NGMIX_CHUNKS + 1):
        out[f"tile_ngmix_{k}"] = _manifest(
            "tile_ngmix", (ngmix_status or {}).get(k, "complete"))
    return out


def _cleaned_store(root, records, survivors=("tile_vignets", "tile_find_exposures")):
    """A reclaimed tile store on disk: the survivors, and the tombstone.

    Exactly what ``clean_tile`` leaves behind — ``logs/`` gone, ``manifests/``
    holding only the manifests other mechanisms own, every record absorbed into
    ``cleaned.json``.
    """
    tile_dir = root / "tiles" / TILE[:2] / TILE
    (tile_dir / "manifests").mkdir(parents=True)
    for stem in survivors:
        (tile_dir / "manifests" / f"{stem}.json").write_text(
            json.dumps(records[stem]))
    (tile_dir / "cleaned.json").write_text(json.dumps({
        "tile": TILE, "manifests": records, "benchmarks": {}}))
    return tile_dir


def _live_store(root, records):
    """The same tile BEFORE reclamation: every manifest, and its identical log."""
    tile_dir = root / "tiles" / TILE[:2] / TILE
    for sub in ("manifests", "logs"):
        (tile_dir / sub).mkdir(parents=True)
        for stem, m in records.items():
            (tile_dir / sub / f"{stem}.json").write_text(json.dumps(m))
    return tile_dir


def _stage_status(run_dir, cleaned_survivors=None):
    """``{stage: status}`` as the report would tally it, from a store on disk."""
    manifests = run_report.load_manifests(run_dir, "tiles")
    if cleaned_survivors is not None:
        run_report.absorb_tombstones(run_dir, "tiles", manifests,
                                     cleaned_survivors)
    return {s: m.get("status") for s, m in manifests[TILE].items()}


@pytest.mark.parametrize("chunk", range(1, NGMIX_CHUNKS + 1))
def test_a_warning_chunk_survives_reclamation(tmp_path, chunk):
    """A warn on ANY chunk reads the same before and after the tile is cleaned.

    Parametrised over all eight because the bug was invisible on chunk 1: that
    is the stem first-key-wins happened to keep, so the one chunk anybody would
    test by hand was the one chunk that worked.
    """
    records = _records({chunk: "warn"})

    live = tmp_path / "live"
    _live_store(live, records)
    before = _stage_status(live)

    cleaned = tmp_path / "cleaned"
    _cleaned_store(cleaned, records)
    after = _stage_status(cleaned, run_report.SURVIVING_TILE_STAGES)

    assert before["tile_ngmix"] == "warn"
    assert after == before


def test_a_failed_chunk_survives_reclamation(tmp_path):
    """Worst-status-wins, not merely warn: failed beats warn beats complete."""
    records = _records({3: "warn", 6: "failed"})
    cleaned = tmp_path / "cleaned"
    _cleaned_store(cleaned, records)
    after = _stage_status(cleaned, run_report.SURVIVING_TILE_STAGES)
    assert after["tile_ngmix"] == "failed"


def test_reclaimed_tile_reports_every_stage(tmp_path):
    """The survivors must not read as a rebuilt chain.

    ``clean_tile`` cannot empty ``manifests/`` — two manifests are currency
    other mechanisms own — so without the survivor set the "records on disk mean
    a rebuilt chain" guard fires and a reclaimed tile reports as one that ran
    two stages and stopped.
    """
    cleaned = tmp_path / "cleaned"
    _cleaned_store(cleaned, _records())
    manifests = run_report.load_manifests(cleaned, "tiles")
    got = run_report.absorb_tombstones(cleaned, "tiles", manifests,
                                       run_report.SURVIVING_TILE_STAGES)
    assert got == {TILE}
    assert sorted(manifests[TILE]) == sorted(run_report.TILE_STAGES)


def test_rebuilt_chain_keeps_the_guard(tmp_path):
    """A unit with a NON-survivor record on disk is a rebuilt chain: skip it.

    This is the default (empty survivor set) that the exposure side uses, and it
    must keep behaving exactly as it did before tiles needed an argument here.
    """
    cleaned = tmp_path / "cleaned"
    tile_dir = _cleaned_store(cleaned, _records())
    (tile_dir / "manifests" / "tile_detect.json").write_text(
        json.dumps(_manifest("tile_detect", "failed")))

    manifests = run_report.load_manifests(cleaned, "tiles")
    got = run_report.absorb_tombstones(cleaned, "tiles", manifests,
                                       run_report.SURVIVING_TILE_STAGES)
    assert got == set()
    assert manifests[TILE]["tile_detect"]["status"] == "failed"
    assert "tile_ngmix" not in manifests[TILE]
