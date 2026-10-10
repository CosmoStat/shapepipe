"""SQLITE STORE.

Lock-free, read-only access to the SqliteDict stores the pipeline writes
between modules (WCS header logs, vignet catalogues, PSF catalogues).

Reading through :class:`sqlitedict.SqliteDict` costs one sqlite transaction
per key access, and each transaction takes and releases a POSIX lock on the
file. On a networked filesystem those locks are slow (tens to hundreds of
milliseconds, and far worse on a busy server), so a loop of keyed reads can
spend most of its time in lock traffic. :class:`ImmutableSqliteDict` opens the
file with sqlite's ``immutable=1`` URI parameter instead, which skips all
locking and change detection.

@sc [label:custody] immutable-sqlite-store
The writer must close the store before a reader opens it, and no process may
write it while an immutable reader is open. SQLite skips change detection
and journal recovery in this mode, so concurrent writes or an unrecovered
journal can yield invalid data. Opening rejects stores with a rollback
journal or write-ahead log beside them; this guard does not prevent a writer
from starting after the check.

"""

import sqlite3
from collections.abc import Mapping
from pathlib import Path

from sqlitedict import decode


class ImmutableSqliteDict(Mapping):
    """Read-only, lock-free mapping over a SqliteDict file.

    Keys and values are those :class:`sqlitedict.SqliteDict` returns for a
    store written with its default key and value encoding (identity keys,
    pickled values); iteration follows insertion (``rowid``) order, as
    SqliteDict's does. Each lookup is one indexed query; use
    :func:`read_sqlitedict` to load a whole store at once.

    The file must not be written while it is open (see the module
    docstring).

    Parameters
    ----------
    path : str or os.PathLike
        Path to an existing SqliteDict file
    tablename : str, optional
        SqliteDict table name, default ``"unnamed"`` (SqliteDict's default)

    Raises
    ------
    FileNotFoundError
        If ``path`` is not an existing file
    RuntimeError
        If a ``-journal`` or ``-wal`` file sits next to ``path``; the store
        may have an active writer or require recovery before an immutable read

    """

    def __init__(self, path, tablename="unnamed"):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"SqliteDict file not found: '{path}'")
        for suffix in ("-journal", "-wal"):
            sidecar = path.with_name(path.name + suffix)
            if sidecar.exists():
                raise RuntimeError(
                    f"SqliteDict file '{path}' has a '{suffix}' file next to"
                    + " it: it is being written, or a writer died"
                    + " mid-transaction. Open it with SqliteDict to roll the"
                    + " journal back, or rewrite it."
                )
        self.path = path
        self._table = '"' + tablename.replace('"', '""') + '"'
        self._conn = sqlite3.connect(
            f"{path.resolve().as_uri()}?immutable=1",
            uri=True,
            check_same_thread=False,
        )

    def __getitem__(self, key):
        row = self._conn.execute(
            f"SELECT value FROM {self._table} WHERE key = ?", (key,)
        ).fetchone()
        if row is None:
            raise KeyError(key)
        return decode(row[0])

    def __contains__(self, key):
        return (
            self._conn.execute(
                f"SELECT 1 FROM {self._table} WHERE key = ?", (key,)
            ).fetchone()
            is not None
        )

    def __iter__(self):
        for (key,) in self._conn.execute(
            f"SELECT key FROM {self._table} ORDER BY rowid"
        ):
            yield key

    def __len__(self):
        return self._conn.execute(
            f"SELECT COUNT(*) FROM {self._table}"
        ).fetchone()[0]

    def items(self):
        """Iterate over ``(key, value)`` pairs with a single query."""
        for key, value in self._conn.execute(
            f"SELECT key, value FROM {self._table} ORDER BY rowid"
        ):
            yield key, decode(value)

    def close(self):
        """Close the underlying sqlite connection."""
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def read_sqlitedict(path, tablename="unnamed"):
    """Load a whole SqliteDict file into a dict without taking locks.

    The file must not be written during the read (see the module docstring).

    Parameters
    ----------
    path : str or os.PathLike
        Path to an existing SqliteDict file
    tablename : str, optional
        SqliteDict table name, default ``"unnamed"``

    Returns
    -------
    dict
        Every entry of the store, in insertion order

    """
    with ImmutableSqliteDict(path, tablename=tablename) as store:
        return dict(store.items())
