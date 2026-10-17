"""Resume a tile after its node-local vignette manifest has been swept."""

import os
import shutil
from collections import Counter
from pathlib import Path

import pytest
from snakemake.io import is_flagged


def test_resume_edge_cut_requires_input_trigger_off(profile):
    """Changing a satisfied chunk's input set must not invalidate it."""
    assert "input" not in profile["rerun-triggers"]


def seed_outputs(campaign, resolve_dag):
    """Materialise satisfied outputs without executing the science pipeline."""
    with resolve_dag(campaign) as dag:
        for job in dag.jobs:
            for output in job.output:
                path = Path(output)
                path.parent.mkdir(parents=True, exist_ok=True)
                if is_flagged(output, "directory"):
                    path.mkdir(exist_ok=True)
                    stamp = path / ".snakemake_timestamp"
                    stamp.touch()
                    os.utime(stamp, (1, 1))
                else:
                    path.touch()
                os.utime(path, (1, 1))


@pytest.mark.parametrize("state", [
    "merged", "chunks_done", "partial", "missing_chunkdir",
])
def test_resume_after_vignets_swept(campaign, resolve_dag, state):
    """Only unfinished chunks read the rebuilt store; make_cat always does."""
    seed_outputs(campaign, resolve_dag)
    tile = next(iter(campaign.ready))
    root = campaign.tile_manifest(tile, "tile_vignets").parents[1]
    campaign.tile_manifest(tile, "tile_vignets").unlink()
    campaign.tile_manifest(tile, "tile_make_cat").unlink()
    campaign.final_cat(tile).unlink()
    if state == "merged":
        # merge_cats reclaims its temp inputs in production.
        for path in (root / "output").glob("run_sp_tile_ngmix_Ng*u"):
            shutil.rmtree(path)
    else:
        campaign.tile_manifest(tile, "tile_merge_cats").unlink()
    if state == "partial":
        campaign.tile_manifest(tile, "tile_ngmix_8").unlink()
    if state == "missing_chunkdir":
        shutil.rmtree(root / "output/run_sp_tile_ngmix_Ng8u")
    with resolve_dag(campaign) as dag:
        pending = list(dag.workflow.dag.needrun_jobs())
        counts = Counter(job.rule.name for job in pending)
        print(f"Resume dry-run jobs ({state}):", dict(sorted(counts.items())))
        assert counts["tile_vignets"] == 1
        assert counts["tile_make_cat"] == 1
        assert counts["tile_ngmix"] == (state in {"partial", "missing_chunkdir"})
        assert counts["tile_merge_cats"] == (state != "merged")
        vignets = next(j for j in pending if j.rule.name == "tile_vignets")
        for job in pending:
            if job.rule.name in {"tile_make_cat", "tile_ngmix"}:
                assert vignets in dag.workflow.dag.dependencies[job]


def test_finished_tile_after_store_only_rebuild(campaign, resolve_dag):
    """A resumed tile must stay finished on the next invocation."""
    seed_outputs(campaign, resolve_dag)
    tile = next(iter(campaign.ready))
    # The store-only rebuild is newer than the preserved chunk products;
    # make_cat then succeeds. alloc_resume skips tiles with a final catalogue.
    campaign.tile_manifest(tile, "tile_vignets").touch()
    with resolve_dag(campaign) as dag:
        outputs = [
            Path(output)
            for job in dag.jobs_for("tile_make_cat")
            if job.wildcards.tile == tile
            for output in job.output
        ]
    # Include any additional durable outputs, such as published tile stats.
    for path in outputs:
        path.touch()
    with resolve_dag(campaign) as dag:
        counts = Counter(j.rule.name for j in dag.workflow.dag.needrun_jobs())
        assert counts["tile_ngmix"] == 0
        assert counts["tile_merge_cats"] == 0
        assert counts["tile_vignets"] == 0
        assert counts["tile_make_cat"] == 0


def test_real_vignets_rebuild_still_invalidates_chunks(campaign, resolve_dag):
    """The storage-only exception must not hide upstream science changes."""
    seed_outputs(campaign, resolve_dag)
    tile = next(iter(campaign.ready))
    campaign.tile_manifest(tile, "tile_detect").touch()
    campaign.tile_manifest(tile, "tile_make_cat").unlink()
    campaign.final_cat(tile).unlink()
    with resolve_dag(campaign) as dag:
        pending = list(dag.workflow.dag.needrun_jobs())
        counts = Counter(j.rule.name for j in pending)
        assert counts["tile_vignets"] == 1
        assert counts["tile_ngmix"] == dag.namespace["NGMIX_CHUNKS"]
