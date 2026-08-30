"""``clean_tile`` deletes a tile's store without ever following a symlink out of it.

``workflow/scripts/clean_tile.py`` reclaims a finished tile by deleting its whole
scratch directory bar four survivors. A finished tile holds NINE symlinks in two
classes, and the second is why this module exists:

  * ``exp_forest/<shard>/<exp>/output`` — 7 links into the exposure stores, each
    shared with 7-10 other tiles;
  * ``output/run_sp_tile_Git/get_images_runner/output/CFIS_{image,weight}-*`` —
    2 links into ``/project/def-mjhudson/unions-wl/tiles``, 621 GB of staged
    survey imaging across 2,536 files on the BACKED-UP, GROUP-SHARED filesystem.

Both are handed wholesale to ``shutil.rmtree`` rather than caught by ``prune``'s
own ``is_symlink()`` test, so the safety argument rests on rmtree's semantics at
a depth ``prune`` never inspects. That is a fine thing to rely on and a bad thing
to leave unpinned: the failure mode is silent, immediate and unrecoverable, and
it would be introduced by an innocent-looking edit to ``prune``.

Deliberately container-free — clean_tile.py is stdlib-only.
"""

import importlib.util
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "clean_tile.py"

TILE = "186.307"
SURVIVORS = (
    "cleaned.json",
    "manifests/tile_vignets.json",
    "manifests/tile_find_exposures.json",
    "output/run_sp_tile_Fe/find_exposures_runner/output/exp_numbers-186-307.txt",
)


def _load():
    """Import the script by path — ``workflow/scripts`` is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    spec = importlib.util.spec_from_file_location("_clean_tile", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


clean_tile = _load()


@pytest.fixture
def store(tmp_path):
    """A finished tile's store, plus the two off-store trees it links into.

    Mirrors the real layout of 186.307 (smk-g4) at the depths that matter: the
    forest links sit three levels down, the Git links four, and ``prune``
    descends into ``output/`` only because the Fe survivor lives there.
    """
    precious = {}
    for name, rel in (("exposures", "exp/21/2114045/output"),
                      ("project_tiles", "unions-wl/tiles")):
        d = tmp_path / name / rel
        d.mkdir(parents=True)
        (d / "DO_NOT_DELETE").write_text("survey data\n")
        precious[name] = d

    tile_dir = tmp_path / "run" / "tiles" / TILE[:2] / TILE
    for rel in SURVIVORS[1:]:                       # the tombstone is written by the job
        p = tile_dir / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("{}" if rel.endswith(".json") else "2114045p\n")

    # bulk that must go
    for rel in ("output/run_sp_tile_Sx/sextractor_runner/output/sexcat.fits",
                "output/run_sp_tile_Uz/uncompress_fits_runner/output/image.fits",
                "logs/tile_detect.json", "manifests/tile_detect.json",
                "tile_numbers.txt"):
        p = tile_dir / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("bulk")

    # class 1: the exposure forest (rmtree'd as a top-level entry)
    forest = tile_dir / "exp_forest" / "21" / "2114045"
    forest.mkdir(parents=True)
    (forest / "output").symlink_to(precious["exposures"])

    # class 2: the get_images links into /project (rmtree'd inside output/)
    git = tile_dir / "output" / "run_sp_tile_Git" / "get_images_runner" / "output"
    git.mkdir(parents=True)
    (git / "CFIS_image-186-307.fits").symlink_to(
        precious["project_tiles"] / "CFIS.186.307.r.fits")
    (git / "CFIS_weight-186-307.fitsfz").symlink_to(
        precious["project_tiles"] / "CFIS.186.307.r.weight.fits.fz")

    return tile_dir, precious


def test_prune_never_follows_a_symlink_out_of_the_store(store):
    """The whole point: reclaiming one tile touches nothing outside that tile."""
    tile_dir, precious = store
    before = {k: sorted(p.rglob("*")) for k, p in precious.items()}

    clean_tile.prune(tile_dir, {tile_dir / r for r in SURVIVORS}, [])

    for name, d in precious.items():
        assert d.is_dir(), f"{name} tree was removed"
        assert sorted(d.rglob("*")) == before[name], f"{name} tree was modified"
        assert (d / "DO_NOT_DELETE").read_text() == "survey data\n"


def test_prune_keeps_exactly_the_survivors(store):
    """Ten inodes: the four survivors and the directories that carry them."""
    tile_dir, _ = store
    clean_tile.prune(tile_dir, {tile_dir / r for r in SURVIVORS}, [])
    left = sorted(str(p.relative_to(tile_dir)) for p in tile_dir.rglob("*"))
    assert left == sorted([
        "manifests", "manifests/tile_vignets.json",
        "manifests/tile_find_exposures.json",
        "output", "output/run_sp_tile_Fe",
        "output/run_sp_tile_Fe/find_exposures_runner",
        "output/run_sp_tile_Fe/find_exposures_runner/output",
        "output/run_sp_tile_Fe/find_exposures_runner/output/exp_numbers-186-307.txt",
    ])


def test_prune_unlinks_a_dangling_link(store):
    """A forest link whose exposure clean_exposure already reclaimed.

    ``exists()`` follows the link and is False for a dangling one, so an
    exists()-first test would skip it and leave the link behind.
    """
    tile_dir, precious = store
    import shutil
    shutil.rmtree(precious["exposures"])
    clean_tile.prune(tile_dir, {tile_dir / r for r in SURVIVORS}, [])
    assert not (tile_dir / "exp_forest").exists()


def test_survivor_contract_refuses_before_deleting_anything(store):
    """A missing survivor aborts with the store intact — it is a contract."""
    tile_dir, _ = store
    (tile_dir / "manifests" / "tile_vignets.json").unlink()
    survivors = clean_tile.survivor_paths(tile_dir, TILE)
    with pytest.raises(SystemExit) as exc:
        clean_tile.require_survivors(survivors, TILE)
    assert "tile_vignets.json" in str(exc.value)
    assert (tile_dir / "output" / "run_sp_tile_Sx").is_dir()


def test_absorption_is_additive_over_an_existing_tombstone(tmp_path):
    """A second clean must not blank the record the first one saved.

    ``script_hash`` reruns this job on any edit to clean_tile.py, so a re-clean
    over an already-pruned store is routine — and it finds only the surviving
    manifests on disk.
    """
    tomb = tmp_path / "cleaned.json"
    tomb.write_text(json.dumps({
        "tile": TILE,
        "manifests": {f"tile_ngmix_{k}": {"stage": "tile_ngmix"} for k in range(1, 9)},
        "benchmarks": {"tile_ngmix_1.benchmark.tsv": {"s": "6852.09"}},
    }))
    manifests, benchmarks = clean_tile.previous_record(tomb)
    assert len(manifests) == 8
    assert benchmarks["tile_ngmix_1.benchmark.tsv"]["s"] == "6852.09"
