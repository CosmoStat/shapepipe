"""The exp_footprint record, and the CCD-index contract it rests on.

``workflow/scripts/exp_footprint.py`` joins two records that name CCDs in two
different ways and never cross-check each other:

  * ``headers-<exp>.npy`` names a CCD by its POSITION in the array — the same
    ``idx-1`` split_exp used for ``image-<exp>-<idx-1>.fits``;
  * ``exp_persist.json`` names a CCD inside a FILENAME,
    ``validation_psf-<exp>-<ccd>.fits``, written by psfex_interp.

The whole design depends on those being the same integer, and no code asserts
it: split_exp writes the image and the array element in one loop, psfex_interp
inherits the numbering through the file handler, and the coverage map would be
silently WRONG — right pixels, wrong exposure count — if they ever diverged by a
permutation. This module is where that contract is pinned. It inherits the role
of ``tests/module/test_coverage.py::test_get_all_shdus``, which pinned the same
``<exp>-<ccd>`` format against the summary scrape that is being retired.

The fixtures are real astropy WCSs, one per CCD with a DISTINCT centre, so a
permutation of the array shows up as corners on the wrong ``id`` rather than as
a passing test. One CCD straddles RA=0 deliberately: the record is raw sky, and
the seam is the map builder's business (``unwrap_ra``), not this script's.

Deliberately not a module test: nothing here runs shapepipe, only its two
geometry helpers, so it belongs with the fast structural suite.
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io.fits import Header
from astropy.wcs import WCS

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "exp_footprint.py"

EXP = "2605805"
# MegaCam-ish: 2048 x 4612 pixels at 0.187"/pixel, i.e. ~0.106 x 0.240 deg.
NX, NY = 2048, 4612
PIXSCALE = 0.187 / 3600.0


def _load():
    """Import the script by path — ``workflow/scripts`` is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    spec = importlib.util.spec_from_file_location("_exp_footprint", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


exp_footprint = _load()


def ccd_header(ra_centre, dec_centre):
    """A single-CCD header + WCS, as split_exp stores them.

    The header carries plain ``NAXIS1/2`` and no ``ZIMAGE``, which is what
    astropy hands back for a tile-compressed HDU once it is decompressed — the
    case ``_image_shape`` takes its second branch for.
    """
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [NX / 2.0, NY / 2.0]
    w.wcs.crval = [ra_centre, dec_centre]
    w.wcs.cdelt = [-PIXSCALE, PIXSCALE]
    header = w.to_header()
    header["NAXIS"] = 2
    header["NAXIS1"] = NX
    header["NAXIS2"] = NY
    return {"WCS": w, "header": header.tostring()}


def headers_array(centres):
    """``headers-<exp>.npy``'s content: object array, one entry per CCD."""
    arr = np.zeros(len(centres), dtype="O")
    for i, (ra, dec) in enumerate(centres):
        arr[i] = ccd_header(ra, dec)
    return arr


# Six CCDs on a row, each 0.3 deg from the last so no two share a footprint, and
# CCD 3 sitting on the RA=0 seam.
CENTRES = [(10.0, 30.0), (10.3, 30.0), (10.6, 30.0),
           (0.0, 30.0), (11.2, 30.0), (11.5, 30.0)]


@pytest.fixture
def store(tmp_path):
    """A scratch exposure store with a headers npy, and a products root."""
    npy_dir = tmp_path / "exp" / EXP / exp_footprint.HEADERS_DIR
    npy_dir.mkdir(parents=True)
    np.save(npy_dir / f"headers-{EXP}.npy", headers_array(CENTRES))
    return tmp_path


def persist_manifest(path, ccds, patterns=("validation_psf-*.fits",)):
    """An ``exp_persist`` manifest naming exactly ``ccds`` as PSF-bearing.

    Shaped as persist_exp.py writes it, including the decoy member: a keep list
    of several patterns packs files this script must ignore, and reading a CCD
    index out of one of them would be a real bug.
    """
    files = [{"name": f"validation_psf-{EXP}-{c}.fits",
              "pattern": "validation_psf-*.fits", "bytes": 1} for c in ccds]
    files.append({"name": f"{EXP}-0.psf", "pattern": "*.psf", "bytes": 1})
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(
        {"stage": "exp_persist", "unit": EXP, "status": "complete",
         "patterns": list(patterns), "n_files": len(files), "files": files},
        indent=2, sort_keys=True))
    return path


def run(store, persist, out, monkeypatch):
    """Invoke the script's CLI, as the rule's shell does."""
    monkeypatch.setattr(sys, "argv", [
        "exp_footprint.py",
        "--exp-dir", str(store / "exp" / EXP),
        "--exp", EXP,
        "--persist-manifest", str(persist),
        "--manifest", str(out)])
    exp_footprint.main()
    return json.loads(out.read_text())


def test_ccd_index_alignment(store, monkeypatch):
    """npy index ``i`` == ``<exp>-<i>`` == ``validation_psf-<exp>-<i>.fits``.

    The three CCDs with a PSF are non-adjacent on purpose: an off-by-one or a
    "renumber the survivors 0..n" bug passes a contiguous fixture and fails
    this one.
    """
    persist = persist_manifest(store / "prod" / "exp_persist.json", [0, 2, 5])
    record = run(store, persist, store / "prod" / "exp_footprint.json",
                 monkeypatch)

    assert record["n_ccd_headers"] == len(CENTRES)
    assert record["n_valid_psf"] == 3
    assert [c["id"] for c in record["ccds"]] == [
        f"{EXP}-0", f"{EXP}-2", f"{EXP}-5"]
    assert record["ccds_no_psf"] == [f"{EXP}-1", f"{EXP}-3", f"{EXP}-4"]

    # The corners on `<exp>-i` are the corners of ARRAY ELEMENT i, not of the
    # i-th surviving CCD: each fixture CCD has its own centre, so this catches
    # any permutation. Checked against the same two helpers the script imports,
    # which is the point — the contract under test is the INDEXING, not the
    # projection arithmetic (tests/module/test_coverage.py owns that).
    arr = headers_array(CENTRES)
    for entry in record["ccds"]:
        i = int(entry["id"].rsplit("-", 1)[1])
        shape = exp_footprint._image_shape(Header.fromstring(arr[i]["header"]))
        ra, dec = exp_footprint._ccd_corners(arr[i]["WCS"], shape)
        assert entry["ra"] == pytest.approx(list(ra))
        assert entry["dec"] == pytest.approx(list(dec))
        # ... and that centre really is CCD i's, not its neighbour's.
        assert np.mean(entry["dec"]) == pytest.approx(CENTRES[i][1], abs=1e-3)


def test_record_is_raw_sky_across_the_ra_seam(store, monkeypatch):
    """CCD 3 straddles RA=0 and the record says so, unwrapped by nobody.

    The seam is handled once, in ``coverage_map_builder.unwrap_ra``, at the
    moment a polygon is stamped. Unwrapping here as well would put negative RA
    into a durable record that other consumers read, and double-unwrapping is
    not idempotent.
    """
    persist = persist_manifest(store / "prod" / "exp_persist.json", [3])
    record = run(store, persist, store / "prod" / "exp_footprint.json",
                 monkeypatch)

    ra = record["ccds"][0]["ra"]
    assert record["ccds"][0]["id"] == f"{EXP}-3"
    assert min(ra) >= 0.0 and max(ra) < 360.0, "raw sky, not a shifted branch"
    assert max(ra) - min(ra) > 180.0, "the fixture must actually cross RA=0"
    # The four corners land on both sides of the seam, ~0.05 deg out.
    assert sorted(round(r) for r in ra) == [0, 0, 360, 360]


def test_byte_stable_rerun_keeps_the_mtime(store, monkeypatch):
    """A rerun over an unchanged store must not move the manifest's mtime.

    mtime is a rerun trigger, and this manifest is an input of clean_exposure:
    an unconditional rewrite would make every downstream reclamation look out
    of date once per invocation. Same contract as persist_exp.py's.
    """
    persist = persist_manifest(store / "prod" / "exp_persist.json", [0, 2, 5])
    out = store / "prod" / "exp_footprint.json"

    first = run(store, persist, out, monkeypatch)
    raw = out.read_bytes()
    before = out.stat().st_mtime_ns

    second = run(store, persist, out, monkeypatch)
    assert second == first
    assert out.read_bytes() == raw
    assert out.stat().st_mtime_ns == before
    # No stray tmp left behind on the persistent root.
    assert not list(out.parent.glob("*.tmp"))


def test_keep_list_without_the_psf_pattern_is_fatal(store, monkeypatch):
    """A manifest that packs no validation_psf files names no valid-PSF set.

    Writing an empty footprint there would be indistinguishable from an
    exposure that genuinely lost every CCD, and the map would silently lose a
    whole exposure's worth of sky.
    """
    persist = persist_manifest(store / "prod" / "exp_persist.json", [],
                               patterns=("*.psf",))
    with pytest.raises(SystemExit) as exc:
        run(store, persist, store / "prod" / "exp_footprint.json", monkeypatch)
    assert "persist_exp" in str(exc.value)


def test_psf_for_a_ccd_the_split_never_wrote_is_fatal(store, monkeypatch):
    """The one disagreement the index alignment cannot absorb, made loud."""
    persist = persist_manifest(store / "prod" / "exp_persist.json", [0, 99])
    with pytest.raises(SystemExit) as exc:
        run(store, persist, store / "prod" / "exp_footprint.json", monkeypatch)
    assert "99" in str(exc.value)


def test_missing_headers_array_is_fatal(store, monkeypatch):
    """A purged or unbuilt split store fails here, never against VOS.

    This is the other half of the rule declaring only its DURABLE input
    (exposure.smk): a persist manifest outliving its scratch store is a real
    state, and it must cost one error line rather than an exposure rebuild.
    """
    npy = (store / "exp" / EXP / exp_footprint.HEADERS_DIR
           / f"headers-{EXP}.npy")
    npy.unlink()
    persist = persist_manifest(store / "prod" / "exp_persist.json", [0])
    with pytest.raises(SystemExit) as exc:
        run(store, persist, store / "prod" / "exp_footprint.json", monkeypatch)
    assert "no WCS array" in str(exc.value)
