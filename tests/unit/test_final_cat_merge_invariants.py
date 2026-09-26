"""What ``final_cat_merge`` owes the catalogues it merges.

Enforces three contracts on a small campaign of tile catalogues, run through
``workflow/scripts/merge_final_cat.py`` as the rule's shell runs it:

* ``final-cat-param-is-exact-allow-list`` (workflow/CONTRACTS): the merged
  datasets carry the param file's columns, in its order, and nothing else;
* ``read-data-raises-on-missing-column`` (create_final_cat.py): a tile short a
  listed column stops the merge, naming the column;
* ``never-fit-rows-pass-through`` (merge_final_cat.py): objects ngmix never fit
  reach the merged file unchanged.

Plus the property the per-epoch families exist for: an object's slot-n tuple
(``EXP_ID_n``, ``CCD_n``, ``HSM_*_PSF_n``) is the same after the merge as in its
tile catalogue. Objects are matched on ``NUMBER``, so a merge may reorder rows
but not move one field of a row without the others.

Failure modes each test was checked against, by editing the code under test:
a listed column dropped or an unlisted one kept (columns); one family's slots
swapped, or one family's rows reordered (per-epoch); never-fit rows dropped or
their sentinels filled (never-fit); the missing-column raise turned into a skip
(missing column).
"""

import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.lib.recfunctions import repack_fields

h5py = pytest.importorskip("h5py")
fits = pytest.importorskip("astropy.io.fits")

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "merge_final_cat.py"

CAMPAIGN = "fixture-campaign"
TILES = ("210.282", "211.282", "212.283")
N_SLOTS = 3
EPOCH_FAMILIES = ("EXP_ID", "CCD", "HSM_G1_PSF", "HSM_G2_PSF")
SCALARS = (
    ("NUMBER", "i4"), ("XWIN_WORLD", "f8"), ("NGMIX_N_EPOCH", "i4"),
    ("NGMIX_MCAL_FLAGS", "i4"), ("NGMIX_G1_NOSHEAR", "f8"),
    ("NGMIX_T_NOSHEAR", "f8"),
)
EPOCH_DTYPE = {"EXP_ID": "i4", "CCD": "i4",
               "HSM_G1_PSF": "f8", "HSM_G2_PSF": "f8"}
# In each tile but not in the param file: the merge must leave these behind.
UNLISTED = (("MAG_UNLISTED", "f4"), ("EXP_ID_4", "i4"), ("CCD_4", "i4"))

PARAM_LIST = [name for name, _ in SCALARS] + [
    f"{fam}_{n}" for n in range(1, N_SLOTS + 1) for fam in EPOCH_FAMILIES]


def _slot(fam, n):
    return f"{fam}_{n}"


def _tile_catalogue(seed, rows=40, n_unfit=5):
    """One tile's catalogue: listed columns, unlisted ones, never-fit rows."""
    rng = np.random.default_rng(seed)
    columns = list(SCALARS) + [
        (_slot(fam, n), EPOCH_DTYPE[fam])
        for n in range(1, N_SLOTS + 1) for fam in EPOCH_FAMILIES
    ] + list(UNLISTED)
    # The catalogue's own column order is not the param file's.
    columns = [columns[i] for i in rng.permutation(len(columns))]
    cat = np.zeros(rows, dtype=columns)
    cat["NUMBER"] = rng.permutation(rows) + 1
    cat["XWIN_WORLD"] = rng.uniform(0, 360, rows)
    cat["MAG_UNLISTED"] = rng.uniform(18, 25, rows)
    cat["EXP_ID_4"] = cat["CCD_4"] = -1
    n_epoch = rng.integers(1, N_SLOTS + 1, rows)
    unfit = rng.choice(rows, n_unfit, replace=False)
    n_epoch[unfit] = 0
    cat["NGMIX_N_EPOCH"] = n_epoch
    cat["NGMIX_G1_NOSHEAR"] = np.where(n_epoch > 0,
                                       rng.normal(0, 0.3, rows), -10.0)
    cat["NGMIX_T_NOSHEAR"] = np.where(n_epoch > 0,
                                      rng.uniform(0.1, 1, rows), 0.0)
    for n in range(1, N_SLOTS + 1):
        used = n_epoch >= n
        # Distinct values across slots and rows, so any exchange shows.
        cat[_slot("EXP_ID", n)] = np.where(
            used, rng.choice(10**7, rows, replace=False), -1)
        cat[_slot("CCD", n)] = np.where(used, rng.integers(0, 40, rows), -1)
        cat[_slot("HSM_G1_PSF", n)] = np.where(
            used, rng.normal(0, 0.05, rows), -10.0)
        cat[_slot("HSM_G2_PSF", n)] = np.where(
            used, rng.normal(0, 0.05, rows), -10.0)
    return cat


def _campaign(root: Path, drop=None):
    """Lay out a campaign the way the workflow does; return (argv, sources).

    ``drop`` removes one listed column from the last tile.
    """
    products = root / "products"
    sources = {}
    for i, tile in enumerate(TILES):
        cat = _tile_catalogue(seed=1000 + i)
        if drop and tile == TILES[-1]:
            keep = [c for c in cat.dtype.names if c != drop]
            cat = repack_fields(cat[keep])
        path = products / "tiles" / tile[:2] / tile / f"final_cat-{tile}.fits"
        path.parent.mkdir(parents=True)
        fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(cat)]).writeto(path)
        sources[tile] = fits.getdata(path, 1)
    tile_list = root / "tiles.txt"
    tile_list.write_text("\n".join(TILES) + "\n")
    index = root / "index.sqlite"
    con = sqlite3.connect(index)
    con.execute("CREATE TABLE tile_exposures(tile_id TEXT, exp_id TEXT)")
    con.executemany("INSERT INTO tile_exposures VALUES (?, ?)",
                    [(t, "2605805") for t in TILES])
    con.commit()
    con.close()
    param = root / "final_cat.param"
    param.write_text("# fixture schema\nNUMBER\n\n# per-epoch\n"
                     + "\n".join(PARAM_LIST[1:]) + "\n")
    output = products / f"final_cat_{CAMPAIGN}.hdf5"
    argv = [sys.executable, str(SCRIPT), "--products-dir", str(products),
            "--tile-list", str(tile_list), "--index-db", str(index),
            "--output", str(output), "--campaign", CAMPAIGN,
            "--param-file", str(param)]
    return argv, output, sources


@pytest.fixture(scope="module")
def merged(tmp_path_factory):
    argv, output, sources = _campaign(tmp_path_factory.mktemp("campaign"))
    run = subprocess.run(argv, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    with h5py.File(output, "r") as f:
        group = f[f"patches/{CAMPAIGN}"]
        out = {tile: group[tile][()] for tile in group}
    return out, sources


def _by_number(arr):
    return arr[np.argsort(arr["NUMBER"])]


def test_columns_are_exactly_the_param_list(merged):
    out, _ = merged
    assert sorted(out) == sorted(TILES)
    for tile, data in out.items():
        assert list(data.dtype.names) == PARAM_LIST, tile


def test_per_epoch_tuples_survive_the_merge(merged):
    out, sources = merged
    for tile, data in out.items():
        got, want = _by_number(data), _by_number(sources[tile])
        assert np.array_equal(got["NUMBER"], want["NUMBER"]), tile
        for n in range(1, N_SLOTS + 1):
            for fam in EPOCH_FAMILIES:
                col = _slot(fam, n)
                assert np.array_equal(got[col], want[col]), (tile, col)


def test_never_fit_rows_pass_through(merged):
    out, sources = merged
    for tile, data in out.items():
        src = sources[tile]
        assert len(data) == len(src), tile
        got = _by_number(data[data["NGMIX_N_EPOCH"] == 0])
        want = _by_number(src[src["NGMIX_N_EPOCH"] == 0])
        assert len(want) > 0, "fixture lost its never-fit rows"
        for col in ("NUMBER", "NGMIX_G1_NOSHEAR", "NGMIX_T_NOSHEAR",
                    "NGMIX_MCAL_FLAGS"):
            assert np.array_equal(got[col], want[col]), (tile, col)


def test_a_tile_missing_a_listed_column_stops_the_merge(tmp_path):
    missing = _slot("HSM_G1_PSF", N_SLOTS)
    argv, output, _ = _campaign(tmp_path, drop=missing)
    run = subprocess.run(argv, capture_output=True, text=True)
    assert run.returncode != 0, "merge succeeded over a tile short a column"
    assert missing in run.stderr
    assert not output.exists()


def _rewrite_tile(tile_path, retype=None):
    """Rewrite one tile's catalogue, optionally narrowing one column to f4."""
    cat = fits.getdata(tile_path, 1)
    arr = np.array(cat)
    if retype:
        dtype = [(n, "f4" if n == retype else arr.dtype[n])
                 for n in arr.dtype.names]
        arr = arr.astype(dtype)
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(arr)]).writeto(
        tile_path, overwrite=True)


def test_one_column_type_per_campaign(tmp_path):
    """A tile rewritten with the same types refreshes alone (FITS byte order
    is not a type change); a tile whose column changes type beside tiles that
    kept the old one is refused, naming the column and both dtypes, and the
    published catalogue is left as it was."""
    argv, output, _ = _campaign(tmp_path)
    assert subprocess.run(argv, capture_output=True).returncode == 0
    tile = lambda t: (output.parent / "tiles" / t[:2] / t
                      / f"final_cat-{t}.fits")

    _rewrite_tile(tile(TILES[0]))
    run = subprocess.run(argv, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert "0 added, 1 refreshed" in run.stdout, run.stdout

    before = output.read_bytes()
    _rewrite_tile(tile(TILES[-1]), retype="NGMIX_T_NOSHEAR")
    run = subprocess.run(argv, capture_output=True, text=True)
    assert run.returncode != 0, "a campaign with two dtypes for one column"
    assert "NGMIX_T_NOSHEAR" in run.stderr
    assert "float32" in run.stderr and "float64" in run.stderr, run.stderr
    assert output.read_bytes() == before
