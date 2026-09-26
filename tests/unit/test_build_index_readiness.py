"""``build_index``: a tile is ready only while its exposure list exists.

The index keeps a tile's edges after its exposure list disappears, because
clean_exposure's consumer sets must still see the tile that read an exposure.
Readiness must not ride on those edges: a tile named in ``missing.json`` is not
ready, for the merges (``campaign_tiles``) or for the Snakefile's
``TILES_READY``, and comes back when its list does.
"""

import importlib.util
import json
import re
import sqlite3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
TILE, OTHER = "210.282", "211.282"


def _load():
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(
            "_build_index", SCRIPTS / "build_index.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


bi = _load()


def _list(run, tile, names):
    path = bi.exp_list_path(run, tile)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{n}\n" for n in names))
    return path


def test_a_missing_tile_is_not_ready_but_keeps_its_edges(tmp_path):
    run, db = tmp_path / "run", tmp_path / "index.sqlite"
    listed = _list(run, TILE, ["2605805p"])
    _list(run, OTHER, ["2605806p"])
    tiles = tmp_path / "tiles.txt"
    tiles.write_text(f"{TILE}\n{OTHER}\n")
    bi.build([TILE, OTHER], run, db)
    assert bi.campaign_tiles(tiles, db) == [TILE, OTHER]

    listed.unlink()
    bi.build([TILE, OTHER], run, db, missing_threshold=1.0)
    assert json.loads((db.parent / "missing.json").read_text()) == [TILE]
    assert bi.campaign_tiles(tiles, db) == [OTHER]
    assert bi.campaign_exposures(tiles, db) == ["2605806"]
    with sqlite3.connect(db) as con:
        edges = con.execute("SELECT exp_id FROM tile_exposures "
                            "WHERE tile_id = ?", (TILE,)).fetchall()
    assert edges == [("2605805",)], "the cleanup consumer edge was dropped"

    _list(run, TILE, ["2605805p"])
    bi.build([TILE, OTHER], run, db)
    assert bi.campaign_tiles(tiles, db) == [TILE, OTHER]


def test_the_snakefile_reads_readiness_from_build_index():
    snakefile = (REPO_ROOT / "workflow" / "Snakefile").read_text()
    line = re.search(r"^TILES_READY = .*$", snakefile, re.M).group(0)
    assert "build_index.ready_tiles(" in snakefile and "_READY_INDEXED" in line
