"""Join diagnostics reach the manifest published beside each tile catalogue."""

import importlib
import json
from pathlib import Path
import subprocess

import pytest

from shapepipe.modules.sextractor_package import match_catalogue as mc
from tests.module.test_match_catalogue import _external, _grid, _sexcat


@pytest.mark.parametrize("n_ext,n_pair", [(845, 774), (2400, 2296), (100, 100), (30, 0)])
def test_match_record_reaches_published_manifest(tmp_path, monkeypatch, n_ext, n_pair):
    scripts = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
    monkeypatch.syspath_prepend(str(scripts))
    completeness = importlib.import_module("completeness")
    unit = "058.220"
    run = tmp_path / unit
    detect_out = run / "output/run_sp_tile_Sx/sextractor_runner/output"
    detect_out.mkdir(parents=True)
    x, y = _grid(n_ext)
    cat = _sexcat(detect_out / "sexcat-058-220.fits", x[:n_pair], y[:n_pair])
    ext = _external(tmp_path / "external.cat", range(1, n_ext + 1), x, y)
    expected = mc.match_catalogue(cat, ext)
    (detect_out / "background.fits").touch()
    detected, ok = completeness.build_manifest("tile_detect", run, unit)
    assert ok
    assert detected["catalogue_match"] == expected
    mdir = run / "manifests"
    mdir.mkdir()
    (mdir / "tile_detect.json").write_text(json.dumps(detected))
    # Chunk aggregation is orthogonal; exercise the real make-cat manifest
    # builder with its already-tested epoch-cuts reader replaced here.
    monkeypatch.setattr(completeness, "gather_manifests", lambda *args: {})
    monkeypatch.setenv("NGMIX_N_CHUNKS", "1")
    output = run / "output/run_sp_tile_Mc/make_cat_runner/output"
    output.mkdir(parents=True)
    (output / "final_cat-058-220.hdf5").touch()
    made, ok = completeness.build_manifest("tile_make_cat", run, unit)
    assert ok
    manifest = mdir / "tile_make_cat.json"
    manifest.write_text(json.dumps(made))
    published = tmp_path / "products/tiles/05/058.220/tile_make_cat.json"
    published.parent.mkdir(parents=True)
    # Execute the publication commands from the actual rule, not a test copy.
    rules = (scripts.parent / "rules/tile.smk").read_text()
    start = rules.index("'  cp -f {output.manifest}")
    stop = rules.index("'", rules.index("'  mv -f", start) + 1)
    import ast
    commands = "".join(ast.literal_eval(line.strip())
                       for line in rules[start:stop + 1].splitlines())
    commands = commands.replace("{output.manifest}", str(manifest))
    commands = commands.replace("{output.stats}", str(published))
    subprocess.run(["bash", "-ec", commands], check=True)
    record = json.loads(published.read_text())
    assert record["unit"] == unit
    assert record["catalogue_match"] == expected
    assert float(record["unit"]) == 58.220  # final catalogue TILE_ID
