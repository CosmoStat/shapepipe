"""UNIT TESTS FOR RUN.

This module contains unit tests for the shapepipe.run module, in
particular the MPI-launcher gating of the mpi4py import (#744): a bare
``shapepipe_run`` must never initialize MPI, otherwise the whole process
aborts inside an ``srun``-launched shell whose Open MPI lacks SLURM PMI
support.

:Author: Claude (on behalf of Cail Daley) <cail.daley@cea.fr>

"""

import os
import subprocess
import sys

import pytest

SNIPPET = "import shapepipe.run as r; print(r.import_mpi)"

# Env vars that either mark an MPI launcher (the gate) or make Open MPI
# believe it was direct-launched by srun (the failure mode under test).
_SCRUBBED_PREFIXES = ("OMPI_", "PMI_", "PMIX_", "SLURM_")


def _import_mpi_flag(extra_env):
    """Report shapepipe.run.import_mpi in a subprocess with a clean env."""
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(_SCRUBBED_PREFIXES)
    }
    env.update(extra_env)
    result = subprocess.run(
        [sys.executable, "-c", SNIPPET],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"subprocess failed (exit {result.returncode}): {result.stderr}"
    )
    return result.stdout.strip()


def test_bare_launch_skips_mpi():
    """A bare launch (no MPI launcher env) must not import/init MPI."""
    assert _import_mpi_flag({}) == "False"


def test_mpirun_launch_imports_mpi():
    """An mpirun-style env (OMPI_COMM_WORLD_SIZE) must import MPI."""
    pytest.importorskip("mpi4py")
    assert _import_mpi_flag({"OMPI_COMM_WORLD_SIZE": "1"}) == "True"
