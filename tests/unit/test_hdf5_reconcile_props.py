"""Property-based state machine over ``workflow/scripts/hdf5_reconcile.py``.

The module's contract is that an hdf5 catalogue reconciled against a campaign
is a FUNCTION OF ITS INPUT SET — the same units with the same sources give the
same datasets, the same dtypes and the same count attribute, however they got
there. That is a claim about every reachable sequence of appends, refreshes and
removals, not about the three the unit tests happen to walk, so it is tested
here against a model: a random sequence of campaign edits, each followed by a
real plan/apply against a real file on disk, with the model asserted after
every step.

The operations are the four things a campaign can do between invocations —
add a unit, change a unit's source, drop a unit, change the column set — plus
a no-op, which is the one that must leave the file's mtime alone.

Source mtimes are set EXPLICITLY with ``os.utime`` rather than left to the
clock. ``stamp()`` is (size, mtime_ns), so a test that rewrote a file with the
same length inside one filesystem tick would silently exercise "nothing
changed" while believing it exercised a refresh.
"""

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, settings
from hypothesis import strategies as st
from hypothesis.stateful import (
    RuleBasedStateMachine,
    initialize,
    invariant,
    precondition,
    rule,
)

h5py = pytest.importorskip("h5py")

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"


def _load(name):
    path = SCRIPTS / f"{name}.py"
    assert path.exists(), f"{path} not found; the rules call it by path"
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(f"_{name}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


reconcile = _load("hdf5_reconcile")

GROUP = "cat/campaign_a"
COUNT_ATTR = "n_units"
UNITS = ["u0", "u1", "u2", "u3"]
# Two column sets, so a schema change is a real change of dtype and width.
COLUMN_SETS = [("RA", "DEC", "E1"), ("RA", "DEC", "E1", "FWHM")]
# What a rebuild may cost over a from-scratch build of the same campaign: the
# metadata h5py's group copy writes for a moved dataset. Measured at ~1.4 kB
# and constant in the number of rebuilds; the allowance is generous because the
# property being defended is "does not grow with history", not an exact size.
COPY_SLACK = 8192


def _array(columns, rows, seed):
    rng = np.random.default_rng(seed)
    dtype = [(c, "<f8") for c in columns]
    out = np.zeros(rows, dtype=dtype)
    for c in columns:
        out[c] = rng.random(rows)
    return out


def _write_source(path: Path, array, mtime_ns: int) -> None:
    np.save(path, array, allow_pickle=False)
    os.utime(path, ns=(mtime_ns, mtime_ns))


def _read(unit, source):
    return np.load(source, allow_pickle=False)


def _build(output: Path, sources: dict, columns):
    """Plan and apply once, exactly as the two merge rules do; return the plan."""
    units = sorted(sources.items())
    digest = reconcile.schema_digest(columns)
    todo = reconcile.plan(output, GROUP, units, digest)
    if todo.empty():
        return todo
    reconcile.apply(output, GROUP, todo, units, _read, digest, COUNT_ATTR)
    return todo


class ReconcileMachine(RuleBasedStateMachine):
    """A campaign that changes under a catalogue that must keep up with it."""

    @initialize()
    def setup(self):
        self.dir = Path(
            __import__("tempfile").mkdtemp(prefix="reconcile-props-")
        )
        self.output = self.dir / "cat.h5"
        self.columns = COLUMN_SETS[0]
        self.sources = {}          # unit -> source path
        self.expected = {}         # unit -> array as last written
        self.clock = 1_000_000_000_000_000_000
        # Compaction is only claimed of the rebuild path (a plan that removes
        # or refreshes). An add-only plan copies the file and appends, so its
        # layout carries whatever the previous writes left behind.
        self.rebuilt = False

    def teardown(self):
        __import__("shutil").rmtree(self.dir, ignore_errors=True)

    # --- the campaign's moves -------------------------------------------
    def _tick(self):
        self.clock += 1_000_000_000
        return self.clock

    def _step(self, changed):
        before = (self.output.stat().st_mtime_ns
                  if self.output.exists() else None)
        todo = _build(self.output, self.sources, self.columns)
        self.rebuilt = bool(todo.remove or todo.refresh)
        if not changed and before is not None:
            assert self.output.stat().st_mtime_ns == before, (
                "a no-op reconcile rewrote the file; mtime is a rerun trigger"
            )

    @rule(pick=st.integers(0, 2**16), rows=st.integers(1, 5),
          seed=st.integers(0, 2**16))
    @precondition(lambda self: len(self.sources) < len(UNITS))
    def add_unit(self, pick, rows, seed):
        free = sorted(set(UNITS) - set(self.sources))
        unit = free[pick % len(free)]
        path = self.dir / f"{unit}.npy"
        array = _array(self.columns, rows, seed)
        _write_source(path, array, self._tick())
        self.sources[unit] = path
        self.expected[unit] = array
        self._step(changed=True)

    @rule(pick=st.integers(0, 2**16), rows=st.integers(1, 5),
          seed=st.integers(0, 2**16), resize=st.booleans())
    @precondition(lambda self: bool(self.sources))
    def modify_source(self, pick, rows, seed, resize):
        unit = sorted(self.sources)[pick % len(self.sources)]
        old = self.expected[unit]
        rows = rows if resize else len(old)
        array = _array(self.columns, rows, seed)
        _write_source(self.sources[unit], array, self._tick())
        self.expected[unit] = array
        self._step(changed=True)

    @rule(pick=st.integers(0, 2**16))
    @precondition(lambda self: bool(self.sources))
    def remove_unit(self, pick):
        unit = sorted(self.sources)[pick % len(self.sources)]
        self.sources.pop(unit).unlink()
        self.expected.pop(unit)
        self._step(changed=True)

    @rule()
    def change_columns(self):
        """Flip to the other column set — a digest change, so every unit
        refreshes."""
        columns = next(c for c in COLUMN_SETS if c != self.columns)
        self.columns = columns
        # A schema change is a change to how the SOURCES are read, so the
        # sources are rewritten under the new column set as the campaign would.
        for i, (unit, path) in enumerate(sorted(self.sources.items())):
            array = _array(columns, len(self.expected[unit]), 4242 + i)
            _write_source(path, array, self._tick())
            self.expected[unit] = array
        self._step(changed=True)

    @rule()
    def no_op(self):
        self._step(changed=False)

    # --- what must be true after every step ------------------------------
    @invariant()
    def file_matches_campaign(self):
        if not self.expected:
            return
        assert self.output.exists()
        with h5py.File(self.output, "r") as f:
            assert set(f[GROUP]) == set(self.expected), (
                "datasets and campaign units disagree")
            assert f.attrs[COUNT_ATTR] == len(self.expected)
            assert (f.attrs["param_digest"]
                    == reconcile.schema_digest(self.columns))
            dtypes = set()
            for unit, want in self.expected.items():
                got = f[GROUP][unit][...]
                assert got.dtype.names == want.dtype.names
                np.testing.assert_array_equal(got, want)
                dtypes.add(got.dtype)
                stamp = reconcile.stamp(self.sources[unit])
                assert (int(f[GROUP][unit].attrs["src_bytes"]),
                        int(f[GROUP][unit].attrs["src_mtime_ns"])) == stamp
            assert len(dtypes) == 1, (
                "sources share a column list; datasets must share a dtype")

    @invariant()
    def compact(self):
        """A rebuild does not carry the old file's dead space forward.

        HDF5 never reclaims a deleted dataset's space, which is why ``apply``
        builds the tmp FRESH whenever a plan removes or refreshes anything
        instead of copying and editing in place. If that path stopped firing,
        a long-lived campaign would grow by one unit per refresh forever.

        The bound is a from-scratch build of the same campaign plus a fixed
        allowance: moving a dataset across with h5py's group copy costs a
        little more metadata than creating it from an array does, measured at
        ~1.4 kB here and — see the cycle test below — independent of how many
        times the file has been rebuilt. What must never hold is growth that
        tracks the history.
        """
        if not self.rebuilt or not self.output.exists():
            return
        fresh = self.dir / "fresh.h5"
        fresh.unlink(missing_ok=True)
        try:
            _build(fresh, self.sources, self.columns)
            if not fresh.exists():
                return
            assert (self.output.stat().st_size
                    <= fresh.stat().st_size + COPY_SLACK), (
                "a rebuilt file is carrying dead space: "
                f"{self.output.stat().st_size} bytes against "
                f"{fresh.stat().st_size} from scratch")
        finally:
            fresh.unlink(missing_ok=True)


ReconcileMachine.TestCase.settings = settings(
    max_examples=150,
    stateful_step_count=14,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
TestReconcileMachine = ReconcileMachine.TestCase


def test_crash_between_tmp_and_replace_leaves_the_file_untouched():
    """A failed rename must leave the previous catalogue byte-identical.

    Not a hypothesis case: the interesting axis is the crash point, and there
    is one. ``os.replace`` is made to raise where the tmp is moved into place.
    """
    import shutil
    import tempfile

    work = Path(tempfile.mkdtemp(prefix="reconcile-crash-"))
    try:
        output = work / "cat.h5"
        columns = COLUMN_SETS[0]
        sources = {}
        for i, unit in enumerate(UNITS[:2]):
            path = work / f"{unit}.npy"
            _write_source(path, _array(columns, 3, i),
                          1_000_000_000_000_000_000 + i)
            sources[unit] = path
        _build(output, sources, columns)
        before = output.read_bytes()
        before_mtime = output.stat().st_mtime_ns

        # A third unit arrives, and the rename fails.
        path = work / "u2.npy"
        _write_source(path, _array(columns, 3, 99), 1_000_000_000_000_000_099)
        sources["u2"] = path
        units = sorted(sources.items())
        digest = reconcile.schema_digest(columns)
        todo = reconcile.plan(output, GROUP, units, digest)
        assert todo.add == ["u2"]

        real_replace = Path.replace

        def boom(self, target):
            raise OSError("simulated crash between write and rename")

        Path.replace = boom
        try:
            with pytest.raises(OSError):
                reconcile.apply(output, GROUP, todo, units, _read, digest,
                                COUNT_ATTR)
        finally:
            Path.replace = real_replace

        assert output.read_bytes() == before, "the old catalogue was modified"
        assert output.stat().st_mtime_ns == before_mtime
        assert not (work / "cat.h5.tmp").exists(), "tmp outlived the failure"
    finally:
        shutil.rmtree(work, ignore_errors=True)


def test_repeated_refresh_does_not_grow_the_file():
    """The leak the rebuild path exists to prevent, asserted directly.

    Twelve refreshes of one unit in a two-unit campaign. If ``apply`` ever
    copied the file and edited it in place, each would strand the previous
    dataset's bytes and the size would climb monotonically.
    """
    import shutil
    import tempfile

    work = Path(tempfile.mkdtemp(prefix="reconcile-growth-"))
    try:
        output = work / "cat.h5"
        columns = COLUMN_SETS[0]
        sources = {}
        for i, unit in enumerate(("u0", "u1")):
            path = work / f"{unit}.npy"
            _write_source(path, _array(columns, 4, i), 10**18 + i)
            sources[unit] = path
            _build(output, sources, columns)

        sizes = []
        for k in range(12):
            _write_source(sources["u0"], _array(columns, 4, 100 + k),
                          10**18 + 100 + k)
            todo = _build(output, sources, columns)
            assert todo.refresh == ["u0"], todo.describe()
            sizes.append(output.stat().st_size)
        assert len(set(sizes)) == 1, f"file size drifted across refreshes: {sizes}"
    finally:
        shutil.rmtree(work, ignore_errors=True)
