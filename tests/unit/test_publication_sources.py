"""Publication waits for durable success or explicit abandonment."""

import importlib.util
from pathlib import Path

import pytest


SPEC = importlib.util.spec_from_file_location(
    "publication_sources", Path(__file__).resolve().parents[2]
    / "workflow/scripts/publication_sources.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


@pytest.mark.parametrize("kind,column,manifest", [
    ("tiles", "tile_id", "tiles/12/123.456/tile_make_cat.json"),
    ("stars", "exp_id", "exp/12/123.456/manifests/exp_persist.json"),
    ("maps", "exp_id", "exp/12/123.456/manifests/exp_maps.json"),
])
def test_missing_manifest_blocks_until_abandoned(tmp_path, kind, column, manifest):
    listing = tmp_path / "list.csv"
    listing.write_text(f"pixel,{column}\n7,123.456\n7,222.333\n8,999.999\n")
    success = tmp_path / "products" / manifest
    success.parent.mkdir(parents=True)
    success.write_text("{}")
    with pytest.raises(RuntimeError, match="222.333"):
        module.finished_sources(7, kind, listing, tmp_path / "products", tmp_path)
    (tmp_path / "abandoned.txt").write_text("222.333 # failed permanently\n")
    assert module.finished_sources(7, kind, listing, tmp_path / "products", tmp_path) == ["123.456"]
    (tmp_path / "abandoned.txt").write_text("123.456\n222.333\n")
    assert module.finished_sources(7, kind, listing, tmp_path / "products", tmp_path) == []


def test_duplicate_list_entry_is_rejected(tmp_path):
    listing = tmp_path / "list.csv"
    listing.write_text("pixel,tile_id\n7,123.456\n7,123.456\n")
    with pytest.raises(ValueError, match="Duplicate"):
        module.finished_sources(7, "tiles", listing, tmp_path, tmp_path)
