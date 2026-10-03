"""Resolved-job checks for campaign scope, product paths, and PSF custody."""

import os
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
PSF_RULES = {"exp_persist", "star_cat_merge", "exp_footprint", "nexp_map"}
DEFECT_RULES = {"exp_defect_map", "defect_map_merge"}   # MAPS_DEFECTS: data only
CATALOGUE_RULES = {"tile_get_catalogue"}


def test_rule_set_matches_input_mode(campaign, dag):
    """A PSF gate cannot remove real-PSF products or add them to fake PSFs;
    the UNIONS-catalogue detection adds its fetch rule and nothing else."""
    expected = BASE_RULES.copy()
    if campaign.psf_model != "fake":
        expected |= PSF_RULES
    if campaign.input_type == "data":
        expected |= DEFECT_RULES
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
        if campaign.input_type == "data":
            expected.append(campaign.exp_manifest(exp, "exp_defect_map"))
        if campaign.psf_model != "fake":
            expected.append(campaign.exp_manifest(exp, "exp_footprint"))
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


def test_campaign_merges_are_local_and_containerized(campaign, dag):
    """Terminal campaign merges run locally with the workflow's image."""
    expected = {"final_cat_merge"}
    if campaign.psf_model != "fake":
        expected.add("star_cat_merge")
    assert {
        rule for rule in expected if dag.jobs_for(rule)
    } == expected
    for rule in expected:
        jobs = dag.jobs_for(rule)
        assert len(jobs) == 1
        job = jobs[0]
        assert job.is_local
        assert job.needs_singularity
        assert job.container_img_url == str(campaign.image)
        assert job.threads == 1
        assert job.rule.group is None


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


def test_missing_run_fails_during_parse(campaign, resolve_dag):
    """Explicit paths cannot bypass the required campaign name diagnostic."""
    campaign.omit_run()
    with pytest.raises(WorkflowError, match=(
        r"for machine='candide', input_type="
        r".*: run\. Set them in your run config \(SP_RUN_CONFIG\)\."
    )):
        with resolve_dag(campaign):
            pytest.fail("a campaign without run: must fail at parse time")


@pytest.mark.parametrize("key", ["coverage", "defect_map"])
def test_retired_map_keys_fail_during_parse(campaign, resolve_dag, key):
    """A run config still carrying a retired block is refused, naming its
    replacement, instead of parsing with the block silently ignored."""
    campaign.config[key] = {"enabled": False}
    campaign.write_config()
    with pytest.raises(WorkflowError, match=rf"{key} -> exposure_maps\."):
        with resolve_dag(campaign):
            pytest.fail(f"a run config with {key}: must fail at parse time")


def test_mccd_is_refused_during_parse(tmp_path, resolve_dag):
    """MCCD products are unreadable to persistence and the star merge."""
    campaign = Campaign(tmp_path / "campaign", "data", "mccd")
    with pytest.raises(WorkflowError, match=r"psf_model=mccd: PSF persistence"):
        with resolve_dag(campaign):
            pytest.fail("psf_model=mccd must be refused at parse time")


# --- blend_handling ---------------------------------------------------------

# The memory each stage adds under uberseg for SEG_VIGNET.
SEG_MEM_MB = {"tile_detect": 3000, "tile_ngmix": 500}
BLEND_EXPORTS = {"tile_detect": {"SP_SEG_VIGNET": "True"},
                 "tile_ngmix": {"SP_BLEND_HANDLING": "uberseg"}}
# The option each stage's committed ini reads the export through, and the
# value it must resolve to with and without it.
BLEND_OPTIONS = {"tile_detect": ("SEG_VIGNET", "False", "True"),
                 "tile_ngmix": ("BLEND_HANDLING", "noisefill", "uberseg")}


def _surface(campaign, resolve_dag):
    """Every job's prologue, shell and memory, keyed by rule and wildcards."""
    with resolve_dag(campaign) as dag:
        return dag.rule_names, {
            (job.rule.name, tuple(sorted(job.wildcards_dict.items()))): (
                getattr(job.params, "pre", None), job.shellcmd,
                job.resources.get("mem_mb"))
            for job in dag.jobs
        }


def _committed_value(shell, config_dir, option, env, monkeypatch):
    """``option`` of the module section of the ini ``shell`` runs, expanded
    under ``env`` as ShapePipe expands it."""
    from shapepipe.pipeline.config import CustomParser

    name = shell.split('shapepipe_run -c "$SP_CONFIG/')[1].split('"')[0]
    parser = CustomParser()
    parser.optionxform = str
    assert parser.read(config_dir / name)
    section = next(s for s in parser.sections() if s.endswith("_RUNNER"))
    for key in ("SP_SEG_VIGNET", "SP_BLEND_HANDLING"):
        monkeypatch.delenv(key, raising=False)
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    return parser.getexpanded(section, option)


@pytest.mark.parametrize("detection", ["sextractor", "unions_catalogue"])
def test_blend_handling_reaches_detection_and_ngmix_only_under_uberseg(
        tmp_path, resolve_dag, monkeypatch, detection):
    """A campaign without the knob plans exactly the uberseg campaign. Against
    explicit noisefill, uberseg only adds its two exports to tile_detect's and
    tile_ngmix's prologues and the seg stamps' memory to both, and each
    export turns its committed ini's option from the noise-fill default to
    the uberseg value."""
    surfaces = {}
    for blend in (None, "noisefill", "uberseg"):
        campaign = Campaign(tmp_path / str(blend), "data", "psfex")
        campaign.config["tile_detection"] = detection
        if blend is not None:
            campaign.config["blend_handling"] = blend
        campaign.write_config()
        surfaces[blend] = _surface(campaign, resolve_dag)

    def normalized(blend):
        rules, jobs = surfaces[blend]
        root = str(tmp_path / str(blend))
        return rules, {key: tuple(v.replace(root, "<ROOT>")
                                  if isinstance(v, str) else v for v in value)
                       for key, value in jobs.items()}

    assert normalized("uberseg") == normalized(None)
    rules, noisefill = normalized("noisefill")
    uberseg_rules, uberseg = normalized("uberseg")
    assert uberseg_rules == rules
    assert uberseg.keys() == noisefill.keys()
    config_dir = Path(__file__).parents[2] / "workflow" / "config" / "cfis"
    for key, (pre, shell, mem) in uberseg.items():
        rule = key[0]
        nf_pre, nf_shell, nf_mem = noisefill[key]
        exports = BLEND_EXPORTS.get(rule, {})
        lines = {f"export {name}='{value}'" for name, value in exports.items()}

        def without_exports(text):
            if text is None:
                return None
            return "\n".join(line for line in text.split("\n")
                             if line not in lines)

        # The rendered shell carries the prologue, so it moves with it and
        # nowhere else.
        assert without_exports(shell) == nf_shell, key
        if pre is None:
            assert nf_pre is None and not exports, key
            continue
        assert lines <= set(pre.split("\n")), key
        assert without_exports(pre) == nf_pre, key
        if exports:
            option, default, value = BLEND_OPTIONS[rule]
            assert _committed_value(shell, config_dir, option, {},
                                    monkeypatch) == default
            assert _committed_value(shell, config_dir, option, exports,
                                    monkeypatch) == value
        assert mem == nf_mem + SEG_MEM_MB.get(rule, 0), key


def test_unknown_blend_handling_fails_during_parse(tmp_path, resolve_dag):
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    campaign.config["blend_handling"] = "mof"
    campaign.write_config()
    with pytest.raises(WorkflowError, match=r"Invalid blend_handling='mof'"):
        with resolve_dag(campaign):
            pytest.fail("an unknown blend_handling must fail at parse time")


INHERITED_BLEND_ENV = {
    f"{prefix}{name}": value
    for prefix in ("", "APPTAINERENV_", "SINGULARITYENV_")
    for name, value in (("SP_SEG_VIGNET", "True"),
                        ("SP_BLEND_HANDLING", "uberseg"))
}


def test_noisefill_ignores_blend_variables_in_the_launch_shell(
        tmp_path, resolve_dag):
    """A noisefill campaign launched from a shell that still exports the
    uberseg variables plans exactly the campaign launched from a clean one,
    and the parse leaves none of them in the environment jobs inherit (the
    slurm executor submits with --export=ALL from this process)."""
    campaigns = {}
    for name in ("clean", "dirty"):
        campaigns[name] = Campaign(tmp_path / name, "data", "psfex")
        campaigns[name].config["blend_handling"] = "noisefill"
        campaigns[name].write_config()
    clean, dirty = campaigns["clean"], campaigns["dirty"]
    _, clean_jobs = _surface(clean, resolve_dag)

    with resolve_dag(dirty, launch_env=INHERITED_BLEND_ENV) as dag:
        leaked = sorted(set(INHERITED_BLEND_ENV) & set(os.environ))
        dirty_jobs = {
            (job.rule.name, tuple(sorted(job.wildcards_dict.items()))): (
                getattr(job.params, "pre", None), job.shellcmd,
                job.resources.get("mem_mb"))
            for job in dag.jobs
        }
    assert leaked == []

    def normalized(jobs, root):
        return {key: tuple(v.replace(str(root), "<ROOT>")
                           if isinstance(v, str) else v for v in value)
                for key, value in jobs.items()}

    assert (normalized(dirty_jobs, tmp_path / "dirty")
            == normalized(clean_jobs, tmp_path / "clean"))
    for _, shell, _ in dirty_jobs.values():
        assert "SP_BLEND_HANDLING" not in (shell or "")
        assert "SP_SEG_VIGNET" not in (shell or "")
