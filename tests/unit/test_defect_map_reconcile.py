"""``merge_defect_map`` agrees the campaign's defect map, or rebuilds it.

The union of the per-exposure defect fragments (CosmoStat/shapepipe#878) is the
one campaign product whose reconciliation is ASYMMETRIC, and that asymmetry is
the reason this file exists. A union cannot be un-OR-ed: two exposures both set
a pixel and nothing in the map records which. So

  * a NEW fragment is OR-ed in on the spot — the cheap, common path;
  * a fragment that LEFT the campaign, or one that CHANGED on disk, forces a
    REBUILD from every fragment;
  * neither leaves the map UNTOUCHED, mtime included, because mtime is a
    rerun trigger and an unconditional rewrite makes every invocation look
    like a change.

Those are four branches of ``reconcile_plan`` whose failure mode is silent: an
edit that stopped treating a removal as a rebuild leaves the map carrying bits
from exposures the campaign no longer has, and nothing downstream — nothing in
this workflow reads the map at all — would ever notice. Hence the pins here.

The SIDECAR is pinned alongside, for the fifth case the plan cannot see: the
record carries two fields about the CAMPAIGN (how many exposures it has, which
of them have no fragment) that can move while the map cannot. The docstring
promises a short map says so on disk; that only holds if a no-op still refreshes
the record.

Fragments are made with ``HealSparseMap.make_empty`` and a handful of pixel ids
— the merge half never opens a FITS image or a WCS, so nothing here needs one.
Needs healsparse, so it runs inside the container and skips outside.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest


pytestmark = pytest.mark.unions

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
SCRIPT = SCRIPTS / "merge_defect_map.py"

NSIDE = 4096
NSIDE_COV = 32


def _load():
    """Import the rule's script by path — scripts/ is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location("_merge_defect_map",
                                                      SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


@pytest.fixture(scope="module")
def merge():
    pytest.importorskip("healsparse")
    pytest.importorskip("h5py")      # hdf5_reconcile, where the stamp lives
    return _load()


def _fragment(merge, root: Path, exp: str, pixels) -> Path:
    """Write one exposure's fragment where ``fragment_path`` expects it."""
    import numpy as np
    import healsparse as hsp

    path = merge.fragment_path(root, exp)
    path.parent.mkdir(parents=True, exist_ok=True)
    frag = hsp.HealSparseMap.make_empty(NSIDE_COV, NSIDE, np.bool_,
                                        bit_packed=True)
    frag[np.asarray(pixels, dtype=np.int64)] = True
    frag.write(str(path), clobber=True)
    return path


def _valid(path: Path):
    import healsparse as hsp
    return set(int(p) for p in hsp.HealSparseMap.read(str(path)).valid_pixels)


def _run(merge, root: Path, have: dict, missing=()):
    """plan + apply, the way ``main`` does, returning (plan, record)."""
    output = root / "defect_map_test.hsp"
    sidecar = root / "defect_map_test.json"
    plan = merge.reconcile_plan(output, sidecar, have)
    if plan.empty():
        return plan, json.loads(sidecar.read_text())
    record = merge.apply_plan(output, sidecar, plan, have, list(missing),
                              NSIDE_COV, NSIDE)
    return plan, record


@pytest.fixture
def campaign(merge, tmp_path):
    """Two exposures, disjoint pixels, merged once. The starting state."""
    have = {
        "2079612p": _fragment(merge, tmp_path, "2079612p", [10, 11, 12]),
        "2079613p": _fragment(merge, tmp_path, "2079613p", [20, 21]),
    }
    plan, _ = _run(merge, tmp_path, have)
    assert plan.rebuild, "the first merge has no map to append to"
    return tmp_path, have


def test_first_merge_is_the_union(merge, campaign):
    root, _ = campaign
    assert _valid(root / "defect_map_test.hsp") == {10, 11, 12, 20, 21}


def test_append_reads_only_the_new_fragment(merge, campaign):
    """A grown campaign is an APPEND, not a rebuild — that is the cheap path."""
    root, have = campaign
    have = dict(have)
    have["2079614p"] = _fragment(merge, root, "2079614p", [30])
    plan, record = _run(merge, root, have)
    assert plan.rebuild == [], plan.describe()
    assert plan.append == ["2079614p"]
    assert _valid(root / "defect_map_test.hsp") == {10, 11, 12, 20, 21, 30}
    assert set(record["exposures"]) == set(have)


def test_removal_forces_a_rebuild_and_drops_the_pixels(merge, campaign):
    """The case a union cannot do incrementally, and the reason for rebuild."""
    root, have = campaign
    have = {k: v for k, v in have.items() if k != "2079613p"}
    plan, _ = _run(merge, root, have)
    assert plan.append == [] and plan.rebuild == ["2079612p"], plan.describe()
    assert "left the campaign" in plan.reason
    assert _valid(root / "defect_map_test.hsp") == {10, 11, 12}


def test_changed_fragment_forces_a_rebuild(merge, campaign):
    """A restamped fragment is not trusted to be a superset of what went in."""
    root, have = campaign
    _fragment(merge, root, "2079613p", [20, 21, 22])
    import os
    os.utime(have["2079613p"], (0, 0))
    plan, _ = _run(merge, root, have)
    assert plan.rebuild == ["2079612p", "2079613p"], plan.describe()
    assert "changed on disk" in plan.reason
    assert _valid(root / "defect_map_test.hsp") == {10, 11, 12, 20, 21, 22}


def test_no_op_leaves_the_map_untouched(merge, campaign):
    """UNTOUCHED, not rewritten identically: mtime is a rerun trigger."""
    root, have = campaign
    output = root / "defect_map_test.hsp"
    before = output.stat().st_mtime_ns
    plan, _ = _run(merge, root, have)
    assert plan.empty(), plan.describe()
    assert output.stat().st_mtime_ns == before


def test_no_op_still_refreshes_a_stale_sidecar(merge, campaign):
    """The campaign moved, the map could not: the RECORD must still say so.

    Tiles whose exposures were all reclaimed by a workflow predating this rule
    add nothing to merge and nothing to remove — an empty plan — but they change
    what the campaign asked for. A sidecar that kept reporting the old counts
    would make a short map look complete on disk.
    """
    root, have = campaign
    sidecar = root / "defect_map_test.json"
    output = root / "defect_map_test.hsp"
    before_map = output.stat().st_mtime_ns
    plan, record = _run(merge, root, have)
    assert plan.empty()

    missing = ["2079999p"]
    plan = merge.reconcile_plan(output, sidecar, have)
    assert plan.empty(), "an exposure with no fragment is not in the plan"
    fresh = merge.build_record(
        have, missing, NSIDE_COV, NSIDE,
        record["n_pixels"], record["n_coverage_pixels"])
    merge.write_sidecar(sidecar, fresh)

    after = json.loads(sidecar.read_text())
    assert after["exposures_without_fragment"] == missing
    assert after["campaign_exposures"] == len(have) + 1
    assert output.stat().st_mtime_ns == before_map, "map must not move"


def test_nside_mismatch_is_an_error_not_an_upgrade(merge, campaign):
    """The ladder's resolution is a campaign decision, not a per-fragment one."""
    import numpy as np
    import healsparse as hsp

    root, have = campaign
    odd = merge.fragment_path(root, "2079615p")
    odd.parent.mkdir(parents=True, exist_ok=True)
    frag = hsp.HealSparseMap.make_empty(NSIDE_COV, NSIDE // 2, np.bool_,
                                        bit_packed=True)
    frag[np.asarray([5], dtype=np.int64)] = True
    frag.write(str(odd), clobber=True)

    have = dict(have)
    have["2079615p"] = odd
    with pytest.raises(SystemExit) as exc:
        _run(merge, root, have)
    assert "re-rasterize 2079615p" in str(exc.value)
