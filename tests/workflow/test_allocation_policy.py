"""Refuse allocation reclamation before generation publication."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from snakemake.exceptions import WorkflowError

from tests.workflow.harness import Campaign, resolve


@pytest.mark.parametrize("selector", ["clean", "clean_tiles"])
def test_allocation_refuses_cleaner_override(tmp_path, monkeypatch, selector):
    """Refuse cleaner selectors that override allocation safety."""
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    campaign.config.update(clean=False, clean_tiles=False)
    campaign.config[selector] = "true"
    campaign.write_config()
    with pytest.raises(
        WorkflowError, match="Allocation mode requires clean=false"
    ):
        with resolve(
            campaign, monkeypatch, environment={"SP_IN_ALLOCATION": "1"}
        ):
            pytest.fail("a cleaner override must fail allocation planning")


def test_exposure_prologue_links_a_fresh_store(tmp_path, resolve_dag):
    """Create the exposure directory before linking local outputs."""
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    exp = sorted(campaign.exposures)[0]
    with resolve_dag(campaign) as dag:
        store = Path(dag.namespace["exp_dir"](exp))
        if store.exists():
            shutil.rmtree(store)
        prologue = dag.namespace["unit_pre"]("exp_get_images", exp)
        local = tmp_path / "node/exp-store"
        result = subprocess.run(
            ["bash", "-c", prologue],
            capture_output=True,
            text=True,
            env=dict(os.environ, SP_EXP_STORE=str(local)),
        )
        assert result.returncode == 0, result.stderr
        assert (store / "output").is_symlink()
        assert (store / "output").resolve() == local / exp
        assert (store / "manifests").is_dir()
