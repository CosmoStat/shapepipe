"""Isolated campaign fixtures for planning-only Snakemake tests."""

import pytest

from tests.workflow.harness import MODES, Campaign, resolve


@pytest.fixture(params=MODES, ids=[f"{mode}+{psf}" for mode, psf in MODES])
def campaign(request, tmp_path):
    """Create a disposable campaign for each supported input/PSF pair."""
    return Campaign(tmp_path / "campaign", *request.param)


@pytest.fixture
def resolve_dag(monkeypatch):
    """Expose the resolver so tests can also assert parse-time failures."""
    return lambda campaign: resolve(campaign, monkeypatch)


@pytest.fixture
def dag(campaign, resolve_dag):
    """Keep the API and its jobs alive for the duration of one check."""
    with resolve_dag(campaign) as resolved:
        yield resolved
