"""Campaign-boundary checks on the shell and params.pre rerun surface."""

import json
from pathlib import Path

from tests.workflow.harness import Campaign
from tests.workflow.params import params_pin

PIN = Path(__file__).with_name("params_pin.json")
BOUNDARY_MESSAGE = (
    "params.pre is a rerun trigger under both profiles; this change reruns "
    "every finished unit of a resumed campaign — land it at a campaign "
    "boundary and update the pin"
)


def test_unit_pre_changes_at_campaign_boundary(psfex_dag, pytestconfig):
    """A shared prologue or per-rule shell edit must move the reviewed pin."""
    actual = params_pin(psfex_dag)
    if pytestconfig.getoption("--update-params-pin"):
        PIN.write_text(json.dumps(actual, indent=2, sort_keys=True) + "\n")
    assert PIN.is_file(), BOUNDARY_MESSAGE
    assert actual == json.loads(PIN.read_text()), BOUNDARY_MESSAGE


def test_params_pin_ignores_fixture_root(tmp_path, resolve_dag):
    """A different temporary campaign directory cannot require a new pin."""
    pins = []
    for directory in ("first-root", "another-root"):
        campaign = Campaign(tmp_path / directory, "data", "psfex")
        with resolve_dag(campaign) as dag:
            pins.append(params_pin(dag))
    assert pins[0] == pins[1], "normalise fixture paths before hashing"


def test_params_is_a_rerun_trigger_under_both_profiles(profile):
    """Neither cluster profile can silently disable the protected trigger."""
    assert "params" in profile["rerun-triggers"], BOUNDARY_MESSAGE
