"""Epoch-cut diagnostics survive cleaning without retaining module logs."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unions
SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
TILE = "259.283"
COUNTS = dict(
    considered=100,
    masked_fraction=7,
    central_veto=3,
    failed=2,
    objects_emptied=1,
)


def load(name):
    """Import a workflow script with its sibling imports available."""
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(
            name, SCRIPTS / f"{name}.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.remove(str(SCRIPTS))


epoch = load("epoch_cuts")
completeness = load("completeness")
clean = load("clean_tile")
backfill = load("backfill_epoch_cuts")


def summary(counts=COUNTS):
    """Render the production end-of-loop summary."""
    return (
        "10/10/2026 07:47:51 epoch cuts: "
        + " ".join(f"{k}={v}" for k, v in counts.items())
        + "\n"
    )


def write(path, text):
    """Create one fixture file and its parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def verdict(stage, **extra):
    """Build a successful completeness verdict."""
    return dict(
        stage=stage,
        unit=TILE,
        status="complete",
        failures=[],
        runners={},
        **extra,
    )


def test_parse_timestamp_and_unrelated_lines():
    """Only the summary supplies counts."""
    assert (
        epoch.parse_epoch_cuts(
            "ngmix loop finished\n" + summary() + "all done\n"
        )
        == COUNTS
    )


def test_parse_pr922_exact_line():
    """Accept the exact field names and spacing emitted by PR #922."""
    assert epoch.parse_epoch_cuts(
        "epoch cuts: considered=8181 masked_fraction=235 "
        "central_veto=271 failed=0 objects_emptied=16"
    ) == dict(
        considered=8181,
        masked_fraction=235,
        central_veto=271,
        failed=0,
        objects_emptied=16,
    )


def test_zero_is_data_not_missing():
    """An explicit zero summary is valid; absent evidence is not."""
    zeros = dict.fromkeys(epoch.FIELDS, 0)
    assert epoch.parse_epoch_cuts(summary(zeros)) == zeros
    with pytest.raises(ValueError, match="found 0"):
        epoch.parse_epoch_cuts("no summary\n")


@pytest.mark.parametrize(
    "text",
    [
        summary() * 2,
        "epoch cuts: considered=1\n",
        summary().replace("failed=2", "failed=-2"),
        summary().replace("failed=2", "failed=2.5"),
        summary().replace("failed=2", "failed=2 failed=2"),
        summary().replace("failed=2", "failed=2 extra=1"),
        summary().replace("considered=100", "considered=1"),
        summary().replace("objects_emptied=1", "objects_emptied=101"),
    ],
)
def test_malformed_and_ambiguous_counts_refused(text):
    """Never salvage partial records or combine retry summaries."""
    with pytest.raises(ValueError):
        epoch.parse_epoch_cuts(text)


@pytest.mark.parametrize("bad", [True, -1, 1.5, "1"])
def test_counts_require_actual_nonnegative_integers(bad):
    """JSON booleans and strings are not counts."""
    with pytest.raises(ValueError):
        epoch.validate_counts(dict(COUNTS, failed=bad))


def test_aggregate_all_fields_once_in_numeric_chunk_order():
    """All five counters sum over exactly the configured chunks."""
    chunks = {str(k): dict(COUNTS) for k in reversed(range(1, 17))}
    record = epoch.aggregate_chunks(chunks, 16)
    assert record["schema_version"] == 1
    assert list(record["chunks"]) == [str(k) for k in range(1, 17)]
    assert record["totals"] == {k: 16 * v for k, v in COUNTS.items()}
    assert record["n_chunks"] == 16


@pytest.mark.parametrize(
    "chunks,n",
    [
        ({"1": COUNTS}, 2),
        ({"1": COUNTS, "3": COUNTS}, 2),
        ({"1": COUNTS, "2": COUNTS}, 1),
        ({}, 0),
        ({}, None),
        ({"1": dict(COUNTS, unknown=0)}, 1),
    ],
)
def test_partial_extra_or_invalid_chunks_never_become_zeros(chunks, n):
    """The chunk set and each count record must be complete."""
    with pytest.raises(ValueError):
        epoch.aggregate_chunks(chunks, n)


def test_chunk_check_reads_worker_not_duplicate_shape_pipe_log(tmp_path):
    """Read one worker, not its stage's duplicate aggregate log."""
    root = tmp_path / TILE
    stage = root / "output" / "run_sp_tile_ngmix_Ng1u"
    write(stage / "ngmix_runner/output/ngmix.fits", "fits fixture")
    write(stage / "ngmix_runner/logs/process-259-283.log", summary())
    write(stage / "logs/log_sp.log", summary())
    record, ok = completeness.build_manifest(
        "tile_ngmix", root, TILE, stage.name
    )
    assert ok and record["epoch_cuts"] == COUNTS
    write(stage / "ngmix_runner/logs/process-other.log", summary())
    record, ok = completeness.build_manifest(
        "tile_ngmix", root, TILE, stage.name
    )
    assert not ok and record["status"] == "failed"
    assert "epoch_cuts" not in record


def test_make_cat_gathers_structured_records_and_preserves_mtime(
    tmp_path, monkeypatch
):
    """Make-cat needs no logs, and identical manifests keep their mtime."""
    root = tmp_path / TILE
    stage = root / "output" / "run_sp_tile_Mc"
    write(stage / "make_cat_runner/output/final_cat.fits", "fits fixture")
    for k in (1, 2):
        write(
            root / f"manifests/tile_ngmix_{k}.json",
            json.dumps(verdict("tile_ngmix", epoch_cuts=COUNTS)),
        )
    monkeypatch.setenv("NGMIX_N_CHUNKS", "2")
    record, ok = completeness.build_manifest("tile_make_cat", root, TILE)
    assert ok
    assert record["epoch_cuts"]["totals"] == {
        k: 2 * v for k, v in COUNTS.items()
    }
    path = tmp_path / "persistent.json"
    text = json.dumps(record, sort_keys=True)
    completeness.write_if_changed(path, text)
    before = path.stat().st_mtime_ns
    completeness.write_if_changed(path, text)
    assert path.stat().st_mtime_ns == before
    write(
        root / "manifests/tile_ngmix_2.json",
        json.dumps(
            dict(verdict("tile_ngmix", epoch_cuts=COUNTS), status="failed")
        ),
    )
    record, ok = completeness.build_manifest("tile_make_cat", root, TILE)
    assert not ok and "epoch_cuts" not in record


@pytest.fixture
def campaign(tmp_path):
    """Create an idle two-chunk campaign with retained logs and symlinks."""
    run = tmp_path / "run"
    products = tmp_path / "products"
    tile = run / "tiles" / TILE[:2] / TILE
    final_cat = write(
        products / "tiles" / TILE[:2] / TILE / f"final_cat-{TILE}.hdf5",
        "persistent catalogue",
    )
    for path in clean.survivor_paths(tile, TILE).values():
        write(path, "{}" if path.suffix == ".json" else "123456p\n")
    write(
        tile / "manifests/tile_make_cat.json",
        json.dumps(verdict("tile_make_cat")),
    )
    for k in (1, 2):
        write(
            tile / f"manifests/tile_ngmix_{k}.json",
            json.dumps(verdict("tile_ngmix")),
        )
        name = f"run_sp_tile_ngmix_Ng{k}u"
        write(
            tile
            / f"logs/modules/{name}/ngmix_runner/logs/process-259-283.log",
            summary(),
        )
        # Original and ShapePipe log copies must not count twice.
        write(
            tile / f"output/{name}/ngmix_runner/logs/process-259-283.log",
            summary(),
        )
        write(tile / f"logs/modules/{name}/logs/log_sp.log", summary())
    write(tile / "output/bulk/data", "x" * 1024)
    precious = write(tmp_path / "outside/precious", "x" * 5000)
    (tile / "output/bulk/link").symlink_to(
        precious.parent, target_is_directory=True
    )
    (tile / "output/bulk/dangling").symlink_to(tmp_path / "absent")
    return run, products, tile, final_cat, precious


def tree_bytes(root):
    """Snapshot fixture bytes and mtimes without reading symlink targets."""
    return {
        str(p): (p.read_bytes(), p.stat().st_mtime_ns)
        for p in root.rglob("*")
        if p.is_file() and not p.is_symlink()
    }


def test_dry_run_is_read_only_and_uses_cleaner_keep_rules(
    campaign, monkeypatch, capsys
):
    """Preview exact whitelist bytes without modifying either filesystem."""
    run, products, tile, cat, precious = campaign
    before = tree_bytes(run.parent)
    monkeypatch.setattr(
        backfill.clean_tile,
        "reclaim",
        lambda *a: pytest.fail("dry-run must never reclaim"),
    )
    assert (
        backfill.main(
            [
                "--run-dir",
                str(run),
                "--products-dir",
                str(products),
                "--n-chunks",
                "2",
                "--dry-run",
            ]
        )
        == 0
    )
    assert tree_bytes(run.parent) == before
    assert not (cat.parent / "tile_make_cat.json").exists()
    # Same traversal as actual cleaning; no symlink target or keep-file bytes.
    keep = set(clean.survivor_paths(tile, TILE).values()) | {
        tile / "cleaned.json"
    }
    expected = sum(
        len(data)
        for path, (data, _) in before.items()
        if Path(path).is_relative_to(tile) and Path(path) not in keep
    )
    output = capsys.readouterr().out
    assert f"WOULD FREE {expected} bytes" in output
    assert "1 ready, 0 skipped, 0 refused" in output
    assert precious.stat().st_size == 5000


def test_backfill_publishes_before_calling_shared_cleaner(
    campaign, monkeypatch
):
    """The real backfill calls only the shared cleaner, after publication."""
    run, products, tile, cat, _ = campaign
    called = []

    def reclaim(tile_dir, tile_id, tomb, final):
        record = epoch.require_durable_manifest(final, tile_id, 2)
        called.append((tile_dir, tomb, record))

    monkeypatch.setattr(backfill.clean_tile, "reclaim", reclaim)
    assert (
        backfill.main(
            [
                "--run-dir",
                str(run),
                "--products-dir",
                str(products),
                "--n-chunks",
                "2",
            ]
        )
        == 0
    )
    assert len(called) == 1
    assert called[0][2]["epoch_cuts"]["totals"] == {
        k: 2 * v for k, v in COUNTS.items()
    }
    path = cat.parent / "tile_make_cat.json"
    before = path.stat().st_mtime_ns
    # Retry reads the durable structured record, no longer requiring logs.
    monkeypatch.setattr(
        backfill, "read_chunk_logs", lambda *a: pytest.fail("already durable")
    )
    backfill.process_tile(tile, products, 2, dry_run=False)
    assert path.stat().st_mtime_ns == before


def test_missing_catalogue_is_skipped_without_reading_logs(
    campaign, monkeypatch, capsys
):
    """Never touch unfinished tiles."""
    run, products, tile, cat, _ = campaign
    # Point to an empty products root; leave the fixture catalogue intact.
    monkeypatch.setattr(
        backfill,
        "read_chunk_logs",
        lambda *a: pytest.fail("unfinished tile must not read logs"),
    )
    status, size, record = backfill.process_tile(
        tile, products / "absent", 2, True
    )
    assert (status, size, record) == ("skipped", 0, None)
    assert "SKIP no persisted final_cat" in capsys.readouterr().out


def test_missing_summary_refuses_without_publishing_or_cleaning(
    campaign, monkeypatch
):
    """A partial calibration record leaves all scratch data intact."""
    run, products, tile, cat, _ = campaign
    write(
        tile
        / "logs/modules/run_sp_tile_ngmix_Ng2u/ngmix_runner/logs"
        / "process-259-283.log",
        "no summary",
    )
    monkeypatch.setattr(
        backfill.clean_tile,
        "reclaim",
        lambda *a: pytest.fail("missing evidence must not clean"),
    )
    assert (
        backfill.main(
            [
                "--run-dir",
                str(run),
                "--products-dir",
                str(products),
                "--n-chunks",
                "2",
                "--dry-run",
            ]
        )
        == 1
    )
    assert not (cat.parent / "tile_make_cat.json").exists()
    assert (tile / "output/bulk/data").exists()


def test_cleaner_refuses_without_durable_counts(campaign):
    """A persisted catalogue alone no longer authorizes cleaning."""
    _, _, tile, cat, _ = campaign
    with pytest.raises(FileNotFoundError):
        clean.reclaim(tile, TILE, tile / "cleaned.json", cat)
    write(
        cat.parent / "tile_make_cat.json", json.dumps(verdict("tile_make_cat"))
    )
    with pytest.raises(ValueError):
        clean.reclaim(tile, TILE, tile / "cleaned.json", cat)
    assert not (tile / "cleaned.json").exists()
    assert (tile / "output/bulk/data").exists()


def test_durable_record_detects_tampered_total(campaign):
    """Validate stored totals against the complete per-chunk evidence."""
    _, _, tile, cat, _ = campaign
    record = backfill.backfill_manifest(tile, cat, 2)
    record["epoch_cuts"]["totals"]["failed"] += 1
    write(cat.parent / "tile_make_cat.json", json.dumps(record))
    with pytest.raises(ValueError, match="invalid epoch cuts aggregate"):
        epoch.require_durable_manifest(cat, TILE)


def test_one_root_refused_before_any_writes(campaign):
    """Products within scratch cannot be durable survivors."""
    run, _, _, _, _ = campaign
    with pytest.raises(SystemExit):
        backfill.main(
            [
                "--run-dir",
                str(run),
                "--products-dir",
                str(run),
                "--n-chunks",
                "2",
                "--dry-run",
            ]
        )
