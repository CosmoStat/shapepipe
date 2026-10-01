"""UNIT TESTS FOR THE LOCK-FREE SQLITEDICT READER.

``ImmutableSqliteDict`` and ``read_sqlitedict`` must return exactly what
``SqliteDict`` returns for the stores the pipeline writes: same keys, same
insertion order, same decoded values.

"""

import pickle
import subprocess
import sys

import numpy as np
import pytest
from sqlitedict import SqliteDict

from shapepipe.modules.merge_headers_package.merge_headers import merge_headers
from shapepipe.pipeline.sqlite_store import ImmutableSqliteDict, read_sqlitedict

from .test_sextractor_post_process import _write_exposure_headers


def _sqlitedict_items(path):
    with SqliteDict(str(path), flag="r") as db:
        return list(db.items())


def _assert_same_items(got, expected):
    assert [key for key, _ in got] == [key for key, _ in expected]
    for (_, got_value), (_, expected_value) in zip(got, expected):
        assert pickle.dumps(got_value) == pickle.dumps(expected_value)


@pytest.fixture
def wcs_store(tmp_path):
    """Tile header log written by merge_headers from split_exp-style files."""
    header_files = []
    for name in ("2104589", "2104590", "2366971"):
        npy_path = tmp_path / f"headers-{name}.npy"
        _write_exposure_headers(npy_path)
        header_files.append([str(npy_path)])
    merge_headers(header_files, str(tmp_path), tile_number="270.283")
    return tmp_path / "log_exp_headers270.283.sqlite"


@pytest.fixture
def object_store(tmp_path):
    """Per-object store with overwritten keys, sentinels and numpy payloads."""
    path = tmp_path / "objects.sqlite"
    with SqliteDict(str(path)) as db:
        db["1"] = {"2104589-12": {"VIGNET": np.arange(9.0).reshape(3, 3)}}
        db["2"] = "empty"
        db["3"] = {}
        db["1"] = {"2104589-13": {"SHAPES": {"HSM_FLAG_PSF": 0}}}
        db.commit()
    return path


def test_read_sqlitedict_matches_sqlitedict_on_wcs_store(wcs_store):
    """The merge_headers WCS log reads back identically, TILE_ID first."""
    expected = _sqlitedict_items(wcs_store)
    got = read_sqlitedict(wcs_store)

    assert isinstance(got, dict)
    assert next(iter(got)) == "TILE_ID"
    _assert_same_items(list(got.items()), expected)
    assert got["2104590"][1]["WCS"].wcs.compare(
        dict(expected)["2104590"][1]["WCS"].wcs
    )


def test_immutable_mapping_matches_sqlitedict(object_store):
    """Lookups, membership, length and rowid order agree with SqliteDict."""
    expected = _sqlitedict_items(object_store)
    with ImmutableSqliteDict(object_store) as store:
        assert list(store) == ["2", "3", "1"]
        assert len(store) == 3
        _assert_same_items(list(store.items()), expected)
        _assert_same_items([(key, store[key]) for key in store], expected)
        assert "2" in store
        assert "4" not in store
        with pytest.raises(KeyError):
            store["4"]


def test_missing_file_raises_without_creating_it(tmp_path):
    """A missing store is an error, and the read leaves no file behind."""
    path = tmp_path / "absent.sqlite"
    with pytest.raises(FileNotFoundError):
        read_sqlitedict(path)
    assert not path.exists()


def test_read_leaves_no_journal(wcs_store):
    """An immutable read writes nothing next to the store."""
    before = sorted(p.name for p in wcs_store.parent.iterdir())
    read_sqlitedict(wcs_store)
    assert sorted(p.name for p in wcs_store.parent.iterdir()) == before


@pytest.mark.parametrize("suffix", ["-journal", "-wal"])
def test_sidecar_journal_refuses_read(object_store, suffix):
    """A rollback journal or WAL next to the store refuses the read."""
    object_store.with_name(object_store.name + suffix).write_bytes(b"")
    with pytest.raises(RuntimeError, match=suffix):
        ImmutableSqliteDict(object_store)


def test_hot_journal_from_killed_writer_refuses_read(tmp_path):
    """A writer killed mid-transaction leaves a hot journal: refuse, not
    return its uncommitted rows."""
    path = tmp_path / "store.sqlite"
    with SqliteDict(str(path)) as db:
        db["1"] = "committed"
        db.commit()
    script = (
        "import os, sqlitedict\n"
        f"db = sqlitedict.SqliteDict({str(path)!r})\n"
        "db['1'] = 'uncommitted'\n"
        "db['2'] = 'uncommitted'\n"
        "db.conn.select_one('SELECT 1')\n"
        "os._exit(0)\n"
    )
    subprocess.run([sys.executable, "-c", script], check=True)
    if not path.with_name(path.name + "-journal").exists():
        pytest.skip("writer left no journal on this filesystem")
    with pytest.raises(RuntimeError, match="-journal"):
        read_sqlitedict(path)
    with SqliteDict(str(path), flag="r") as db:
        assert dict(db.items()) == {"1": "committed"}


def test_path_with_uri_special_characters(tmp_path):
    """Space, '?', '#' and '%' in the path do not break the file: URI."""
    directory = tmp_path / "a dir?x=1#frag%20"
    directory.mkdir()
    path = directory / "st ore?#%.sqlite"
    with SqliteDict(str(path)) as db:
        db["k"] = {"v": 1}
        db.commit()
    assert read_sqlitedict(path) == {"k": {"v": 1}}
