"""Reclaimed fragments are products, not rasterization targets."""

import re
import shutil

import pytest

from tests.workflow import harness
from tests.workflow.harness import Campaign
from workflow.scripts import build_index


def test_clean_live_historical_exposure(tmp_path, resolve_dag, monkeypatch):
    """Cleanup must build a missing fragment outside the ready tiles' scope."""
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    exp = "2243883"
    assert exp not in campaign.exposures
    build_index.exp_list_path(campaign.run_dir, campaign.outside).write_text(
        exp + "p\n")
    build_index.build([campaign.outside], campaign.run_dir, campaign.index_db)
    persist = campaign.persist_manifest(exp)
    persist.parent.mkdir(parents=True)
    persist.write_text("{}")
    store = campaign.run_dir / "exp" / exp[:2] / exp
    (store / "output/run_sp_exp_Sp/split_exp_runner/output").mkdir(parents=True)
    (store / "manifests").mkdir()
    (store / "manifests/exp_psf.json").write_text("{}")
    fragment = persist.parent.parent / "maps" / f"maps-{exp}.hsp"
    assert not fragment.exists()
    with resolve_dag(campaign) as dag:
        clean, = [job for job in dag.jobs_for("clean_exposure")
                  if job.wildcards.exp == exp]
        maps, = [job for job in dag.jobs_for("exp_maps")
                 if job.wildcards.exp == exp]
        assert dag.workflow.dag.needrun(maps)
        assert str(fragment) in maps.output
        assert maps in dag.workflow.dag.dependencies[clean]
        merge, = dag.jobs_for("exposure_maps")
        assert str(fragment) not in merge.input
        # Workers must retain the historical producer even though only the
        # head schedules rule all's cleanup targets.
        monkeypatch.setattr(type(dag.workflow), "is_main_process",
                            property(lambda self: False))
        dag.namespace["clean_exposures"].cache_clear()
        assert dag.namespace["clean_targets"]() == []
        assert re.fullmatch(dag.namespace["exp_maps_live_pattern"](), exp)


def test_live_missing_fragment_is_built_before_merge_and_reclamation(
        tmp_path, resolve_dag):
    """A surviving map manifest cannot hide a missing live fragment."""
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    for exp in campaign.exposures:
        persist = campaign.persist_manifest(exp)
        persist.parent.mkdir(parents=True)
        persist.write_text('{}')
        campaign.persist_manifest(exp, "exp_maps").write_text('{}')
        store = campaign.run_dir / "exp" / exp[:2] / exp
        split = store / "output/run_sp_exp_Sp/split_exp_runner/output"
        split.mkdir(parents=True)
        (store / "manifests").mkdir()
        (store / "manifests/exp_psf.json").write_text('{}')
    with resolve_dag(campaign) as dag:
        merge, = dag.jobs_for("exposure_maps")
        jobs = dag.jobs_for("exp_maps")
        assert len(jobs) == len(campaign.exposures)
        for job in jobs:
            assert dag.workflow.dag.needrun(job)
            assert job in dag.workflow.dag.dependencies[merge]
            clean, = [j for j in dag.jobs_for("clean_exposure")
                      if j.wildcards.exp == job.wildcards.exp]
            assert job in dag.workflow.dag.dependencies[clean]
            assert "--image-mef" not in job.shellcmd


@pytest.mark.parametrize("change", ["script", "shell", "params"])
@pytest.mark.parametrize("map_manifest", [False, True])
def test_reclaimed_fragments_have_no_producer(
        tmp_path, monkeypatch, resolve_dag, change, map_manifest):
    """Code/params edits cannot schedule or delete irreplaceable fragments."""
    # Mutate an isolated workflow copy, never the checkout used by other tests.
    repo = tmp_path / "code"
    for directory in ("workflow", "profiles", "scripts", "src"):
        shutil.copytree(harness.REPO / directory, repo / directory,
                        ignore=shutil.ignore_patterns("__pycache__"))
    monkeypatch.setattr(harness, "REPO", repo)
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    fragments = []
    for exp in campaign.exposures:
        persist = campaign.persist_manifest(exp)
        persist.parent.mkdir(parents=True)
        persist.write_text('{}')
        if map_manifest:
            campaign.persist_manifest(exp, "exp_maps").write_text('{}')
        store = campaign.run_dir / "exp" / exp[:2] / exp
        store.mkdir(parents=True)
        (store / "cleaned.json").write_text('{}')
        fragment = persist.parent.parent / "maps" / f"maps-{exp}.hsp"
        fragment.parent.mkdir()
        fragment.write_bytes(f"irreplaceable fragment {exp}".encode())
        fragments.append(fragment)
    for tile in campaign.ready:
        final = campaign.final_cat(tile)
        final.parent.mkdir(parents=True)
        final.touch()
        campaign._manifest(tile, "tile_vignets")
    raw = campaign.root / "inputs" / campaign.name / "exposures"
    assert not raw.exists()
    before = [(path.read_bytes(), path.stat().st_mtime_ns)
              for path in fragments]

    def check():
        with resolve_dag(campaign) as dag:
            assert not dag.jobs_for("exp_maps")
            producer = dag.workflow.get_rule("exp_maps")
            merge, = dag.jobs_for("exposure_maps")
            for path in fragments:
                assert str(path) in set(map(str, merge.input))
                assert not producer.is_producer(str(path))
            assert dag.namespace["exposure_maps_exposures"]() == list(
                campaign.exposures)
        assert [(path.read_bytes(), path.stat().st_mtime_ns)
                for path in fragments] == before

    check()
    if change == "script":
        path = repo / "workflow/scripts/exp_maps.py"
        path.write_text(path.read_text() + "\n# Simulated rasterizer change\n")
    else:
        path = repo / "workflow/rules/exposure.smk"
        text = path.read_text()
        if change == "shell":
            text = text.replace(
                'f"python {SCRIPTS}/exp_maps.py"',
                '"echo changed rasterizer shell >&2\\n"\n'
                '        f"python {SCRIPTS}/exp_maps.py"')
        else:
            text = text.replace("script_hash = EXP_MAPS_HASH",
                                "script_hash = EXP_MAPS_HASH + '-changed'")
        assert text != path.read_text()
        path.write_text(text)
    check()
