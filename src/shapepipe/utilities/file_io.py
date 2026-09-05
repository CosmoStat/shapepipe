"""FILE I/O UTILITIES.

Small, dependency-light helpers for publishing files the rest of the pipeline
treats as a cache.

:Author: consolidated from workflow/scripts/star_cats.py and
         scripts/python/create_star_cat.py

"""

import os


def write_atomic(table, path):
    """Publish ``table`` at ``path`` all-or-nothing.

    Every caller's only cache test is existence (``Path.exists`` /
    ``os.path.isfile``), so a write killed part-way -- job timeout, OOM, node
    failure -- would otherwise leave a truncated FITS that every later run
    trusts forever, and ``test -s`` passes on partial bytes. Writing to a temp
    and renaming makes the visible file all-or-nothing: ``os.replace`` is atomic
    within a directory.

    The temp keeps the target's suffix, because astropy picks its writer from
    the extension. It is dot-prefixed and PID-tagged so it stays out of the
    ``star_chunk-*`` / ``star_cat*`` globs the rules use, and two concurrent
    writers cannot collide.

    Parameters
    ----------
    table : astropy.table.Table
        Table to write
    path : str or pathlib.Path
        Destination path

    """
    path = os.fspath(path)
    directory, name = os.path.split(path)
    tmp = os.path.join(directory or ".", f".tmp-{os.getpid()}-{name}")
    try:
        table.write(tmp, overwrite=True)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
