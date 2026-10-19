"""Every resumed reader must wait for the rebuilt node-local vignette store."""

import os
import shutil
from collections import Counter
from pathlib import Path

import pytest
from snakemake.io import is_flagged


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
    "merged", "merged_reclaimed", "chunks_done", "partial", "missing_chunkdir",
])
@pytest.mark.parametrize("invalidate_detect", [False, True], ids=[
    "swept_only", "swept_and_upstream_invalidated",
])
def test_resume_after_vignets_swept(
    campaign, resolve_dag, state, invalidate_detect,
):
    """Sweep alone and sweep plus upstream invalidation both preserve ordering.

    The newer tile_detect input reproduces the review's counterexample:
    completed chunks can still be scheduled and must retain their writer edge.
    """
    seed_outputs(campaign, resolve_dag)
    tile = next(iter(campaign.ready))
    root = campaign.tile_manifest(tile, "tile_vignets").parents[1]
    campaign.tile_manifest(tile, "tile_vignets").unlink()
    campaign.tile_manifest(tile, "tile_make_cat").unlink()
    campaign.final_cat(tile).unlink()
    if state == "merged_reclaimed":
        # merge_cats reclaims its temp inputs in production.
        for path in (root / "output").glob("run_sp_tile_ngmix_Ng*u"):
            shutil.rmtree(path)
    elif state != "merged":
        campaign.tile_manifest(tile, "tile_merge_cats").unlink()
    if state == "partial":
        campaign.tile_manifest(tile, "tile_ngmix_8").unlink()
    if state == "missing_chunkdir":
        shutil.rmtree(root / "output/run_sp_tile_ngmix_Ng8u")
    if invalidate_detect:
        campaign.tile_manifest(tile, "tile_detect").touch()
    with resolve_dag(campaign) as dag:
        pending = list(dag.workflow.dag.needrun_jobs())
        counts = Counter(job.rule.name for job in pending)
        print(
            f"Resume dry-run jobs ({state}, invalidate_detect={invalidate_detect}):",
            dict(sorted(counts.items())),
        )
        assert counts["tile_vignets"] == 1
        assert counts["tile_ngmix"] == dag.namespace["NGMIX_CHUNKS"]
        assert counts["tile_merge_cats"] == 1
        assert counts["tile_make_cat"] == 1
        writer = next(j for j in pending if j.rule.name == "tile_vignets")
        readers = [
            job for job in pending
            if job.rule.name in {"tile_make_cat", "tile_ngmix"}
        ]
        assert len(readers) == dag.namespace["NGMIX_CHUNKS"] + 1
        for job in readers:
            assert job.wildcards.tile == tile
            assert writer in dag.workflow.dag.dependencies[job], (
                f"pending reader {job.rule.name} {job.wildcards} has no writer "
                f"dependency; inputs={list(job.input)}"
            )
