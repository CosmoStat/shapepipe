"""Shared pytest configuration for the ShapePipe test suite.

Markers, environment detection, and the candide skip policy live here so
every test module — wherever it sits in the tree — sees the same rules.

Markers (also declared in ``pyproject.toml`` so ``--strict-markers`` is on):

* ``slow``    — heavy compute (minutes), not part of the fast inner loop.
* ``candide`` — needs the candide cluster and/or its real on-disk data;
  meaningless (and auto-skipped) anywhere else.

A ``candide``-marked test only runs on a candide node. Everywhere else it
is skipped with a clear reason, so the same suite is green on a laptop, in
CI, and on the cluster — the cluster-only tests simply do not fire off it.
"""

import os
import re
import socket

import pytest

from hypothesis import settings


settings.register_profile("ci", derandomize=True, max_examples=50)
settings.register_profile("dev", max_examples=200)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "ci"))


# --------------------------------------------------------------------------- #
# Candide detection
# --------------------------------------------------------------------------- #

# Candide compute / login nodes are named c01, c03, n22..n36, etc. The login
# host this suite is most often driven from is ``c03``. We match the candide
# node-name families rather than a fixed list so new nodes are covered, and
# allow an explicit override for CI or odd hostnames.
_CANDIDE_HOST_RE = re.compile(r"^(c\d|n\d{2})", re.IGNORECASE)


def on_candide():
    """Return True when running on a candide node.

    The check is, in order: an explicit ``SHAPEPIPE_ON_CANDIDE`` override
    (``1``/``0``), then the hostname against the candide node-name families
    (``c0x`` login, ``nXX`` compute). Cheap, import-safe, no cluster calls.
    """
    override = os.environ.get("SHAPEPIPE_ON_CANDIDE")
    if override is not None:
        return override == "1"
    return bool(_CANDIDE_HOST_RE.match(socket.gethostname()))


# --------------------------------------------------------------------------- #
# Marker registration + skip policy
# --------------------------------------------------------------------------- #


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: heavy compute (minutes); excluded from the fast inner loop.",
    )
    config.addinivalue_line(
        "markers",
        "candide: needs the candide cluster and/or its real data; "
        "auto-skipped elsewhere.",
    )


def pytest_collection_modifyitems(config, items):
    """Skip ``candide`` tests off-cluster.

    Collection still happens everywhere — so ``pytest --collect-only`` shows
    the cluster tests exist — they are just marked skipped at run time when
    not on candide.
    """
    if on_candide():
        return
    skip_candide = pytest.mark.skip(
        reason="needs candide (set SHAPEPIPE_ON_CANDIDE=1 to force)"
    )
    for item in items:
        if "candide" in item.keywords:
            item.add_marker(skip_candide)


# --------------------------------------------------------------------------- #
# Shared fixtures
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="session")
def artifacts_dir():
    """Directory where guardrail tests drop plots + status summaries.

    The artifacts SEAM: a later GitHub Pages / status step publishes from
    here. Created on demand so a clean checkout has nothing to commit until
    a test actually emits.
    """
    from pathlib import Path

    path = Path(__file__).parent / "tests" / "_artifacts"
    path.mkdir(parents=True, exist_ok=True)
    return path
