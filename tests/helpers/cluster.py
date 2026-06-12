"""Small SLURM submission helper for candide cluster tests.

Cluster guardrail tests have two phases: *regenerate* (submit the ShapePipe
job chain via SLURM, wait for the outputs) and *evaluate* (run the science
assertion against an outputs directory). This module owns the first phase so
the test bodies stay about the science.

Nothing here imports shapepipe — it is pure ``subprocess`` around ``srun`` /
``sbatch`` / ``squeue`` / ``sacct``, so it is import-safe off-cluster and the
tests that use it are guarded by the ``candide`` marker anyway.

The default star-grid job chain (per Fabian's handoff, 2026-06-12):

    tile_launcher.job
      -> job_per_tile_newversion.job   (sbatch, one per tile)
           -> SP_1z2z_star_grid/job_sp_14.job   (through the container)

A full regeneration is a multi-hour, multi-job affair; this helper is the
seam to drive it on demand, NOT something the test suite runs by default.
"""

import shutil
import subprocess
import time
from pathlib import Path


SP_SIMU_ROOT = Path("/home/hervas/n25/SP_simu_fab")
STAR_GRID_ROOT = SP_SIMU_ROOT / "SP_1z2z_star_grid"


def slurm_available():
    """True when the SLURM client tools are on PATH (i.e. on candide)."""
    return shutil.which("sbatch") is not None and shutil.which("squeue") is not None


def srun(args, *, partition="comp", time_limit="00:30:00", cpus=4, extra=None):
    """Run a command synchronously on a compute node via ``srun``.

    Blocks until the step finishes. Use for short, single-step cluster work
    (e.g. running the R-function evaluation against existing outputs on a
    node). Returns the ``CompletedProcess``; raises on non-zero exit.

    Parameters
    ----------
    args : list of str
        Command + arguments to execute on the node.
    partition : str
        SLURM partition (``comp`` / ``pscomp`` on candide).
    time_limit : str
        Wall-clock limit, ``HH:MM:SS``.
    cpus : int
        ``--cpus-per-task``.
    extra : list of str, optional
        Extra ``srun`` flags (e.g. ``["--exclude", "n09,n17,n36"]``).
    """
    cmd = [
        "srun",
        f"--partition={partition}",
        f"--time={time_limit}",
        f"--cpus-per-task={cpus}",
        *(extra or ["--exclude", "n09,n17,n36"]),
        *args,
    ]
    return subprocess.run(cmd, check=True, text=True, capture_output=True)


def sbatch(job_script, *job_args, dependency=None):
    """Submit a batch job and return its job id (str).

    Parameters
    ----------
    job_script : str or Path
        Path to the ``.job`` / ``.sbatch`` script.
    *job_args
        Positional arguments appended after the script path.
    dependency : str, optional
        A SLURM dependency spec, e.g. ``afterok:12345``.
    """
    cmd = ["sbatch", "--parsable"]
    if dependency is not None:
        cmd.append(f"--dependency={dependency}")
    cmd.extend([str(job_script), *map(str, job_args)])
    out = subprocess.run(cmd, check=True, text=True, capture_output=True).stdout
    # --parsable prints "<jobid>" or "<jobid>;<cluster>"
    return out.strip().split(";")[0]


def submit_star_grid_chain(tiles, *, launcher=None):
    """Submit the star-grid job chain for a list of tile ids.

    Mirrors ``tile_launcher.job``: one ``job_per_tile_newversion.job`` per
    tile id. Returns the list of submitted job ids. Does NOT wait — pair
    with :func:`wait_for_jobs`.

    This is the heavy path (each tile is a multi-hour job). It is wired so a
    cluster test can regenerate outputs on demand; the suite never calls it
    automatically.
    """
    launcher = Path(launcher) if launcher else SP_SIMU_ROOT / "job_per_tile_newversion.job"
    if not launcher.exists():
        raise FileNotFoundError(f"launcher job script not found: {launcher}")
    return [sbatch(launcher, tile) for tile in tiles]


def wait_for_jobs(job_ids, *, poll=60, timeout=4 * 3600):
    """Block until every job id has left the SLURM queue.

    Polls ``squeue`` for the given ids. Returns when none remain pending or
    running, or raises ``TimeoutError`` after ``timeout`` seconds.
    """
    deadline = time.monotonic() + timeout
    ids = set(map(str, job_ids))
    while time.monotonic() < deadline:
        out = subprocess.run(
            ["squeue", "--jobs", ",".join(ids), "--noheader", "--format=%i"],
            check=False,
            text=True,
            capture_output=True,
        ).stdout
        still_running = {line.strip() for line in out.splitlines() if line.strip()}
        if not (still_running & ids):
            return
        time.sleep(poll)
    raise TimeoutError(f"jobs {ids} did not finish within {timeout}s")
