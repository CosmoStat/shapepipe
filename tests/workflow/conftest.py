"""Isolated campaign fixtures for planning-only Snakemake tests."""

import pytest

from tests.workflow.harness import MODES, Campaign, load_profile, resolve


def pytest_addoption(parser):
    """Expose the explicit campaign-boundary pin update switch."""
    parser.addoption(
        "--update-params-pin", action="store_true", default=False,
        help="Update the reviewed params.pre/shell pin at a campaign boundary",
    )


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


@pytest.fixture
def psfex_dag(tmp_path, resolve_dag):
    """Resolve the canonical data+psfex campaign for the prologue pin."""
    with resolve_dag(Campaign(tmp_path / "campaign", "data", "psfex")) as dag:
        yield dag


@pytest.fixture(params=["candide", "nibi"])
def profile(request):
    """Read each profile's actual rerun-trigger policy."""
    return load_profile(request.param)
