"""The workflow's coverage map: exposure footprint records -> a HealSparse nexp map.

``tests/module/test_coverage.py`` pins the nexp contract through the
``exp_ra_dec.txt`` CLI path. This module pins the SAME contract through the path
the Snakemake workflow actually takes — one ``exp_footprint.json`` per exposure
on the products root, globbed and fed to ``build_map`` as arrays — because that
path has two properties the text one does not and neither is visible in the map:

  * it is CAMPAIGN-CUMULATIVE by construction. The script globs every record
    under ``<products_dir>/exp/*/*/manifests/``, sharded, including exposures
    whose scratch stores were reclaimed. A regression that fed it only the
    declared rule inputs would build a map of one batch and look fine.
  * the records are RAW SKY. Unwrapping across RA=0 happens once, inside
    ``build_map``; a record written pre-unwrapped, or unwrapped twice, silently
    moves a footprint by 360 degrees.

Two exposures, one CCD each, offset so they overlap — the same fixture geometry
as the CLI test, so the two routes are comparable by eye.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "coverage_map.py"


def _load():
    """Import the script by path — ``workflow/scripts`` is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    spec = importlib.util.spec_from_file_location("_coverage_map", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


coverage_map = _load()


def write_footprint(products, exp, ccds):
    """One exposure's record, at the sharded path exp_footprint writes to."""
    path = (products / "exp" / exp[:2] / exp / "manifests"
            / "exp_footprint.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "stage": "exp_footprint", "level": "exp", "unit": exp,
        "status": "complete",
        "n_ccd_headers": len(ccds), "n_valid_psf": len(ccds),
        "ccds": [{"id": f"{exp}-{i}", "ra": ra, "dec": dec}
                 for i, (ra, dec) in enumerate(ccds)],
        "ccds_no_psf": [],
    }, indent=2, sort_keys=True))
    return path


def box(ra_min, ra_max, dec_min, dec_max):
    """A rectangular CCD footprint, corners counterclockwise from bottom-left."""
    return ([ra_min, ra_max, ra_max, ra_min],
            [dec_min, dec_min, dec_max, dec_max])


def run(products, tmp_path, monkeypatch, nside="1024"):
    out = tmp_path / "coverage" / "coverage.hsp"
    manifest = tmp_path / "coverage" / "manifests" / "coverage_map.json"
    monkeypatch.setattr(sys, "argv", [
        "coverage_map.py",
        "--products-dir", str(products),
        "--out", str(out), "--manifest", str(manifest),
        "--nside-coverage", "32", "--nside", nside])
    coverage_map.main()
    import healsparse as hsp
    return hsp.HealSparseMap.read(str(out)), json.loads(manifest.read_text())


def test_two_exposures_accumulate_to_nexp(tmp_path, monkeypatch):
    """Overlapping CCDs from two exposures give 2 in the overlap, 1 outside.

    The exposures live in DIFFERENT shard directories, which is the campaign
    layout: a glob that assumed one shard would find only one of them.
    """
    products = tmp_path / "products"
    write_footprint(products, "1000001", [box(9.8, 10.2, 19.8, 20.2)])
    write_footprint(products, "2000001", [box(10.0, 10.4, 19.8, 20.2)])

    m, manifest = run(products, tmp_path, monkeypatch)

    values = m[m.valid_pixels]
    assert values.max() == 2
    assert values.min() == 1
    assert (values == 2).sum() > 0
    assert (values == 1).sum() > 0

    assert manifest["n_exposures"] == 2
    assert manifest["exposures"] == ["1000001", "2000001"]
    assert manifest["n_ccds"] == 2
    assert manifest["nside_coverage"] == 32


def test_every_record_on_the_root_is_used(tmp_path, monkeypatch):
    """The map is campaign-cumulative: nothing selects a subset of the records.

    A third exposure appears on the products root with no involvement from any
    caller — the state a reclaimed or out-of-scope exposure is in — and it must
    still be in the map.
    """
    products = tmp_path / "products"
    write_footprint(products, "1000001", [box(9.8, 10.2, 19.8, 20.2)])
    write_footprint(products, "2000001", [box(10.0, 10.4, 19.8, 20.2)])
    write_footprint(products, "3000001", [box(40.0, 40.4, 19.8, 20.2)])

    m, manifest = run(products, tmp_path, monkeypatch)

    assert manifest["exposures"] == ["1000001", "2000001", "3000001"]
    # The far-away exposure is real sky in the map, not just a row in the record.
    assert m.get_values_pos(40.2, 20.0, lonlat=True) == 1


def test_seam_record_is_unwrapped_by_the_builder(tmp_path, monkeypatch):
    """A raw-sky record across RA=0 lands where it belongs, not 360 deg away."""
    products = tmp_path / "products"
    write_footprint(products, "1000001",
                    [([359.8, 0.2, 0.2, 359.8], [19.8, 19.8, 20.2, 20.2])])

    m, _ = run(products, tmp_path, monkeypatch)

    assert m.get_values_pos(0.0, 20.0, lonlat=True) == 1
    assert m.get_values_pos(359.9, 20.0, lonlat=True) == 1
    # Nothing was stamped on the far side of the sky.
    assert m.get_values_pos(180.0, 20.0, lonlat=True) == 0


def test_no_records_is_a_loud_failure(tmp_path, monkeypatch):
    """An empty products root must not write an empty map and exit green."""
    products = tmp_path / "products"
    products.mkdir()
    with pytest.raises(SystemExit) as exc:
        run(products, tmp_path, monkeypatch)
    assert "nothing to build a map from" in str(exc.value)


def test_records_naming_no_ccd_is_a_loud_failure(tmp_path, monkeypatch):
    """Records that name no CCD are the quieter empty map, and equally fatal.

    Every exposure on the root having lost every CCD is a broken PSF stage, not
    a survey with no coverage — but the .hsp it would write is valid and
    plausible, and its consumer would mask everything.
    """
    products = tmp_path / "products"
    write_footprint(products, "1000001", [])
    write_footprint(products, "2000001", [])

    with pytest.raises(SystemExit) as exc:
        run(products, tmp_path, monkeypatch)
    assert "not one names a CCD with a PSF model" in str(exc.value)
