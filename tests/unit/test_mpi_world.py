"""Guard against the silent "rank 0 of N singletons" MPI failure.

When the host MPI launcher and the container's MPI/PMIx stack are
incompatible, every process initialises standalone (``COMM_WORLD`` size 1),
and ShapePipe would otherwise run N uncoordinated copies of the pipeline into
the same output directory -- silently, with exit code 0.
``check_mpi_world`` turns that into an immediate error by comparing the size
that actually wired up against the launcher's intended ``OMPI_COMM_WORLD_SIZE``.
The intended/actual pairs below are the values measured on candide for a
healthy run and for the OpenMPI-4-container / OpenMPI-5-host mismatch.
"""

import pytest

from shapepipe.pipeline.mpi_run import check_mpi_world


class _FakeComm:
    """Minimal stand-in exposing only ``Get_size``."""

    def __init__(self, size):

        self._size = size

    def Get_size(self):

        return self._size


@pytest.mark.parametrize(
    "intended, actual",
    [
        ("4", 4),  # healthy multi-node run
        (None, 1),  # no launcher env -> legitimate single-rank run
        ("1", 1),  # explicit single rank
    ],
    ids=["healthy", "no-launcher-env", "single-rank"],
)
def test_check_mpi_world_passes(monkeypatch, intended, actual):

    if intended is None:
        monkeypatch.delenv("OMPI_COMM_WORLD_SIZE", raising=False)
    else:
        monkeypatch.setenv("OMPI_COMM_WORLD_SIZE", intended)

    check_mpi_world(_FakeComm(actual))


@pytest.mark.parametrize(
    "intended, actual",
    [
        ("4", 1),  # the measured singleton failure
        ("4", 3),  # partial wire-up
    ],
    ids=["singletons", "partial-wireup"],
)
def test_check_mpi_world_aborts_on_mismatch(monkeypatch, intended, actual):

    monkeypatch.setenv("OMPI_COMM_WORLD_SIZE", intended)

    with pytest.raises(RuntimeError, match="MPI world mismatch"):
        check_mpi_world(_FakeComm(actual))
