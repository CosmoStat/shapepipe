"""Resolved-job checks for campaign scope, product paths, and PSF custody."""

import re
from collections import Counter
from pathlib import Path

import pytest
from snakemake.exceptions import WorkflowError

from tests.workflow.harness import Campaign

BASE_RULES = {
    "all", "tile_get_images", "tile_uncompress", "tile_find_exposures",
    "exp_get_images", "exp_split", "exp_psf", "clean_exposure",
    "tile_exp_forest", "tile_merge_headers", "tile_detect", "tile_vignets",
    "tile_ngmix", "tile_merge_cats", "tile_make_cat", "clean_tile",
    "final_cat_merge",
}
TILE_SHAPE_STORE_RULES = ("tile_vignets", "tile_ngmix", "tile_make_cat")
PSF_RULES = {"exp_persist", "star_cat_merge", "exp_maps", "exposure_maps"}
CATALOGUE_RULES = {"tile_get_catalogue"}


def test_rule_set_matches_input_mode(campaign, dag):
    """A PSF gate cannot remove real-PSF products or add them to fake PSFs;
    the UNIONS-catalogue detection adds its fetch rule and nothing else."""
    expected = BASE_RULES.copy()
    if campaign.psf_model != "fake":
        expected |= PSF_RULES
    if campaign.tile_detection == "unions_catalogue":
        expected |= CATALOGUE_RULES
    assert dag.rule_names == expected
    assert "merge_final_cats" not in dag.declared_rule_names


def test_tile_detect_joins_the_catalogue_iff_unions(campaign, dag):
    """Data's tile_detect waits on the fetch and exports its catalogue as
    SP_MATCH_CATALOGUE; the image-simulation prologue exports it empty."""
    for job in dag.jobs_for("tile_detect"):
        tile = job.wildcards.tile
        inputs = {str(f) for f in job.input}
        fetch = str(campaign.tile_manifest(tile, "tile_get_catalogue"))
        if campaign.tile_detection == "unions_catalogue":
            gic = (campaign.run_dir / "tiles" / tile[:2] / tile / "output"
                   / "run_sp_tile_Gic" / "get_images_runner" / "output")
            cat = gic / f"CFIS_cat-{tile.replace('.', '-')}.cat"
            assert fetch in inputs
            assert f"export SP_MATCH_CATALOGUE='{cat}'" in job.params.pre
        else:
            assert fetch not in inputs
            assert "export SP_MATCH_CATALOGUE=''" in job.params.pre.split("\n")


def test_clean_exposure_waits_on_persist_iff_psf(campaign, dag):
    """Reclamation waits for persistence and exactly its in-scope readers."""
    jobs = dag.jobs_for("clean_exposure")
    assert {job.wildcards.exp for job in jobs} == set(campaign.exposures)
    for job in jobs:
        exp = job.wildcards.exp
        expected = [
            campaign.tile_manifest(tile, "tile_vignets")
            for tile, exposures in campaign.ready.items() if exp in exposures
        ]
        if campaign.psf_model != "fake":
            expected.append(campaign.persist_manifest(exp))
            expected.append(campaign.persist_manifest(exp, "exp_maps"))
        assert Counter(map(str, job.input)) == Counter(map(str, expected)), (
            "clean-exposure-waits-on-persist-iff-psf", exp, list(job.input)
        )


def test_final_cat_merge_reads_every_ready_tile(campaign, dag):
    """A merge cannot drop a ready tile or pull one from outside this batch."""
    jobs = dag.jobs_for("final_cat_merge")
    assert len(jobs) == 1
    assert Counter(map(str, jobs[0].input)) == Counter(
        str(campaign.final_cat(tile)) for tile in campaign.ready
    )


def test_products_use_products_dir_and_run_name(campaign, dag):
    """Neither scratch nor a directory basename can name durable products."""
    expected = {
        "tile_make_cat": {
            campaign.final_cat(tile) for tile in campaign.ready
        },
        "final_cat_merge": {
            campaign.products_dir / f"final_cat_{campaign.name}.hdf5"
        },
        "star_cat_merge": set(),
        "exp_persist": set(),
    }
    if campaign.psf_model != "fake":
        expected["star_cat_merge"] = {
            campaign.products_dir / f"full_starcat_{campaign.name}.hdf5"
        }
        expected["exp_persist"] = {
            campaign.persist_manifest(exp) for exp in campaign.exposures
        }
    output_names = {
        "tile_make_cat": "final_cat", "final_cat_merge": "merged",
        "star_cat_merge": "star_cat", "exp_persist": "manifest",
    }
    for rule, paths in expected.items():
        actual = {
            Path(getattr(job.output, output_names[rule]))
            for job in dag.jobs_for(rule)
        }
        assert actual == paths, (rule, actual, paths)
        assert all(
            path.is_relative_to(campaign.products_dir) for path in actual
        )
    for job in dag.jobs_for("exp_persist"):
        assert Path(job.params.dest) == (
            campaign.products_dir / "exp" / job.wildcards.exp[:2]
            / job.wildcards.exp / "psf"
        )
    for rule in ("final_cat_merge", "star_cat_merge"):
        for job in dag.jobs_for(rule):
            assert job.params.campaign == campaign.name
            assert f"--campaign '{campaign.name}'" in job.shellcmd
    assert dag.namespace["CAMPAIGN"] == campaign.name
    assert Path(dag.namespace["INDEX_DB"]) == campaign.index_db


def test_tile_store_is_unique_per_campaign(campaign, tmp_path, resolve_dag):
    """Concurrent campaigns over the same tiles get distinct node-local stores;
    within a campaign, every tile_shape member of a tile names the same one."""
    other = Campaign(tmp_path / "another-campaign", campaign.input_type,
                     campaign.psf_model)
    stores = []
    for each in (campaign, other):
        by_tile = {}
        with resolve_dag(each) as dag:
            for rule in TILE_SHAPE_STORE_RULES:
                for job in dag.jobs_for(rule):
                    match = re.search(r'^export SP_LOCAL="([^"]+)"$',
                                      job.params.pre, re.MULTILINE)
                    assert match, (rule, job.wildcards_dict)
                    by_tile.setdefault(job.wildcards.tile, set()).add(
                        match.group(1))
        assert by_tile.keys() == set(each.ready)
        assert all(len(paths) == 1 for paths in by_tile.values()), by_tile
        stores.append({tile: paths.pop() for tile, paths in by_tile.items()})
    first, second = stores
    assert all(first[tile] != second[tile] for tile in first), (first, second)


def test_module_logs_outlive_their_run_dirs(campaign, dag):
    """tile_vignets (node-local) and every ngmix chunk (temp()) copy their
    module logs into the tile's logs/modules/ on the shared root, after
    shapepipe_run and before the rule exits, whatever its rc."""
    sources = {"tile_vignets": lambda job: '"$NGMIX_VIGNET_DIR"',
               "tile_ngmix": lambda job: f'"{job.output.chunkdir}"'}
    for rule, source in sources.items():
        for job in dag.jobs_for(rule):
            name = ("run_sp_tile_PiViVi" if rule == "tile_vignets" else
                    f"run_sp_tile_ngmix_Ng{job.wildcards.chunk}u")
            lines = job.shellcmd.splitlines()
            keep = [i for i, line in enumerate(lines)
                    if f'"$SP_RUN/logs/modules/{name}"' in line]
            assert len(keep) == 1, (rule, job.wildcards_dict)
            line = lines[keep[0]]
            assert f"cd {source(job)}" in line
            assert '-path "*/logs/*"' in line and line.rstrip().endswith(">&2")
            run = next(i for i, l in enumerate(lines)
                       if l.startswith("shapepipe_run "))
            assert run < keep[0] < len(lines) - 1
            assert lines[-1] == "exit $rc"


def test_missing_run_fails_during_parse(campaign, resolve_dag):
    """Explicit paths cannot bypass the required campaign name diagnostic."""
    campaign.omit_run()
    with pytest.raises(WorkflowError, match=(
        r"for machine='candide', input_type="
        r".*: run\. Set them in your run config \(SP_RUN_CONFIG\)\."
    )):
        with resolve_dag(campaign):
            pytest.fail("a campaign without run: must fail at parse time")


def test_mccd_is_refused_during_parse(tmp_path, resolve_dag):
    """MCCD products are unreadable to persistence and the star merge."""
    campaign = Campaign(tmp_path / "campaign", "data", "mccd")
    with pytest.raises(WorkflowError, match=r"psf_model=mccd: PSF persistence"):
        with resolve_dag(campaign):
            pytest.fail("psf_model=mccd must be refused at parse time")
