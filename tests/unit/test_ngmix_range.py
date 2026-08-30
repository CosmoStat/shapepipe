"""The ngmix chunk split tiles ``[1, n_obj]`` exactly.

Nothing downstream cross-checks the chunk ranges: merge_sep_cats concatenates
the chunk catalogues, so an overlap silently duplicates objects and a gap
silently drops them. The one invariant that catches both is that the ranges
partition ``1..n_obj`` — that is what this module asserts, over the epoch-weight
distributions the splitter actually has to cope with.

``workflow/scripts/ngmix_range.py`` now runs ONCE per tile (tile_vignets writes
the whole partition to the group's node-local scratch) and is then only read
from, once per chunk. The second thing this module pins is that the trip through
JSON changes nothing: row ``k`` of the written file is exactly what asking the
splitter for chunk ``k`` used to return.

Deliberately container-free: the range function takes a plain sequence of epoch
counts, so nothing here imports astropy, numpy, or shapepipe.
"""

import importlib.util
import json
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "ngmix_range.py"


def _load():
    """Import the script by path — ``workflow/scripts`` is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    spec = importlib.util.spec_from_file_location("_ngmix_range", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ngmix_range = _load()


def assert_tiles(ranges, n_obj, n_chunks):
    """Assert the partition invariant, in the form each failure mode takes."""
    assert len(ranges) == n_chunks
    assert ranges[0][0] == 1
    assert ranges[-1][1] == n_obj
    for (lo, hi), (prev_lo, prev_hi) in zip(ranges[1:], ranges):
        assert lo == prev_hi + 1
        assert hi >= lo - 1, "an empty range is (hi + 1, hi), never narrower"
        assert prev_hi >= prev_lo - 1
    covered = [i for lo, hi in ranges for i in range(lo, hi + 1)]
    assert covered == list(range(1, n_obj + 1))
    # ID_OBJ_MAX <= 0 is ngmix's "unbounded" sentinel: a chunk emitting it
    # re-measures the whole tile instead of nothing.
    assert all(hi >= 1 for _, hi in ranges)


# --------------------------------------------------------------------------- #
# The property
# --------------------------------------------------------------------------- #


@settings(deadline=None)
@given(
    n_obj=st.integers(min_value=1, max_value=400),
    n_chunks=st.integers(min_value=1, max_value=16),
    data=st.data(),
)
def test_ranges_tile_the_catalogue(n_obj, n_chunks, data):
    """Any epoch distribution, any shape: the ranges still partition 1..N."""
    epochs = data.draw(
        st.lists(
            st.integers(min_value=0, max_value=40),
            min_size=n_obj,
            max_size=n_obj,
        )
    )
    assert_tiles(ngmix_range.id_ranges(epochs, n_chunks), n_obj, n_chunks)


@settings(deadline=None)
@given(
    n_obj=st.integers(min_value=1, max_value=60),
    n_chunks=st.integers(min_value=1, max_value=8),
    data=st.data(),
)
def test_ranges_are_deterministic(n_obj, n_chunks, data):
    """Same input, same ranges — the eight processes share nothing else.

    Integer weights make the split a pure function of the epoch list, with no
    float accumulation whose order could matter; this pins that.
    """
    epochs = data.draw(
        st.lists(
            st.integers(min_value=0, max_value=99),
            min_size=n_obj,
            max_size=n_obj,
        )
    )
    first = ngmix_range.id_ranges(epochs, n_chunks)
    assert first == ngmix_range.id_ranges(list(epochs), n_chunks)
    assert all(isinstance(b, int) for lo, hi in first for b in (lo, hi))


# --------------------------------------------------------------------------- #
# Write-once / read-per-chunk is the same partition
# --------------------------------------------------------------------------- #


def _roundtrip(epochs, n_chunks):
    """What a chunk shell sees: the doc as it comes back off disk."""
    return json.loads(json.dumps(ngmix_range.partition(epochs, n_chunks)))


@settings(deadline=None)
@given(
    n_obj=st.integers(min_value=1, max_value=200),
    n_chunks=st.integers(min_value=1, max_value=12),
    data=st.data(),
)
def test_written_rows_reproduce_the_per_chunk_computation(n_obj, n_chunks, data):
    """THE INVARIANCE TEST for materialising the split.

    Writing the whole partition and reading chunk ``k``'s row back must give
    byte-for-byte what the old per-chunk ``id_ranges(...)[k - 1]`` returned. If
    this ever fails, the two mechanisms have drifted and a tile's coverage is
    what pays.
    """
    epochs = data.draw(
        st.lists(
            st.integers(min_value=0, max_value=40),
            min_size=n_obj,
            max_size=n_obj,
        )
    )
    doc = _roundtrip(epochs, n_chunks)
    assert doc["n_obj"] == n_obj and doc["n_chunks"] == n_chunks
    expected = ngmix_range.id_ranges(epochs, n_chunks)
    got = [ngmix_range.chunk_range(doc, k) for k in range(1, n_chunks + 1)]
    assert got == expected
    assert_tiles(got, n_obj, n_chunks)


def test_chunk_index_outside_the_written_partition_is_fatal():
    """A chunk index the file does not hold must not index from the end."""
    doc = _roundtrip([3] * 40, 4)
    for bad in (0, -1, 5):
        with pytest.raises(SystemExit, match="outside 1..4"):
            ngmix_range.chunk_range(doc, bad)


def test_reading_an_absent_ranges_file_is_fatal(tmp_path):
    """NO FALLBACK: a missing file fails, it does not recompute.

    Recomputing per chunk is precisely the mechanism this design removed.
    """
    with pytest.raises(SystemExit, match="no ranges file"):
        ngmix_range.read_ranges(tmp_path / "ngmix_ranges.json")


# --------------------------------------------------------------------------- #
# Degenerate shapes — each pins a decision, not just an absence of crash
# --------------------------------------------------------------------------- #


def test_single_chunk_takes_everything():
    """n_chunks == 1: one range over the whole catalogue."""
    assert ngmix_range.id_ranges([3] * 17, 1) == [(1, 17)]


def test_one_object_per_chunk():
    """n_obj == n_chunks: one object each, however lopsided the weights."""
    assert ngmix_range.id_ranges([0, 9, 1, 40], 4) == [
        (1, 1), (2, 2), (3, 3), (4, 4)
    ]


def test_fewer_objects_than_chunks_pads_with_empty_ranges():
    """DOCUMENTED: the surplus chunks get ``(n_obj + 1, n_obj)``, not (1, 0).

    The old equal-count split gave every surplus chunk ``(1, 0)``, and ngmix
    reads ``ID_OBJ_MAX = 0`` as unbounded — so each of them would have measured
    the entire tile rather than nothing.
    """
    assert ngmix_range.id_ranges([2, 5, 1], 6) == [
        (1, 1), (2, 2), (3, 3), (4, 3), (4, 3), (4, 3)
    ]
    assert_tiles(ngmix_range.id_ranges([2, 5, 1], 6), 3, 6)


def test_zero_objects_is_fatal():
    """No split is meaningful, and (1, 0) is ngmix's unbounded sentinel."""
    with pytest.raises(ValueError, match="zero objects"):
        ngmix_range.id_ranges([], 8)


def test_zero_chunks_is_fatal():
    """There is no zeroth chunk to hand a range to."""
    with pytest.raises(ValueError, match="n_chunks"):
        ngmix_range.id_ranges([1, 2, 3], 0)


def test_equal_weights_reproduce_an_equal_count_split():
    """Uniform epochs: chunk sizes differ by at most one object."""
    ranges = ngmix_range.id_ranges([3] * 1000, 8)
    assert_tiles(ranges, 1000, 8)
    sizes = [hi - lo + 1 for lo, hi in ranges]
    assert max(sizes) - min(sizes) <= 1
    assert sum(sizes) == 1000


def test_zero_epoch_objects_still_weigh_something():
    """ALPHA is why: without it a 0-epoch tail is free and lands in one chunk.

    Ninety-six zero-epoch objects and four 1-epoch ones. Weighing only epochs
    would let a single chunk swallow all ninety-six.
    """
    ranges = ngmix_range.id_ranges([0] * 96 + [1] * 4, 4)
    assert_tiles(ranges, 100, 4)
    assert max(hi - lo + 1 for lo, hi in ranges) < 96


def test_one_enormously_heavy_object_is_isolated():
    """The heavy object gets a chunk to itself; the tail still tiles."""
    epochs = [1] * 20 + [10_000] + [1] * 20
    ranges = ngmix_range.id_ranges(epochs, 4)
    assert_tiles(ranges, 41, 4)
    assert (21, 21) in ranges


# --------------------------------------------------------------------------- #
# The split is the optimal one, not merely a legal one
# --------------------------------------------------------------------------- #


def _brute_force_min_max(weights, n_chunks):
    """Smallest achievable maximum chunk weight, by exhaustive enumeration."""
    n = len(weights)
    best = {}

    def solve(start, chunks):
        if chunks == 1:
            return sum(weights[start:])
        if (start, chunks) not in best:
            best[(start, chunks)] = min(
                max(sum(weights[start:cut]), solve(cut, chunks - 1))
                for cut in range(start + 1, n - chunks + 2)
            )
        return best[(start, chunks)]

    return solve(0, n_chunks)


@settings(deadline=None, max_examples=40)
@given(
    epochs=st.lists(
        st.integers(min_value=0, max_value=12), min_size=4, max_size=14
    ),
    n_chunks=st.integers(min_value=2, max_value=4),
)
def test_slowest_chunk_is_minimal(epochs, n_chunks):
    """The objective is the slowest chunk — check it against brute force."""
    if len(epochs) < n_chunks:
        return
    weights = [
        ngmix_range.MILLI_EPOCH * e + ngmix_range.ALPHA_MILLI_EPOCHS
        for e in epochs
    ]
    ranges = ngmix_range.id_ranges(epochs, n_chunks)
    loads = [sum(weights[lo - 1:hi]) for lo, hi in ranges]
    assert max(loads) == _brute_force_min_max(weights, n_chunks)
