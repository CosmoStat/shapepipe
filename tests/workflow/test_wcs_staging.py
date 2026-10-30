"""Concurrent containers stage WCS without sharing PID-derived temp names."""

import os
import subprocess


def test_wcs_staging_survives_shared_container_pid(psfex_dag, tmp_path):
    """Two PID namespaces can give their shells the same $$.

    A barrier makes both copies finish before either rename, exposing reuse
    of a PID-only name instead of relying on a probabilistic race.
    """
    tile = next(iter(psfex_dag.campaign.ready))
    prologue = psfex_dag.namespace["tile_local"](tile)
    store = tmp_path / "store"
    store.mkdir()
    prologue = prologue.replace("/local/scratch", str(store)).replace("$$", "14")
    run = tmp_path / "run"
    source_dir = run / "output/run_sp_tile_Mh_exp/merge_headers_runner/output"
    source_dir.mkdir(parents=True)
    source = source_dir / f"log_exp_headers-{tile.replace('.', '-')}.sqlite"
    source.write_bytes(b"complete WCS store\n" * 1000)

    tools = tmp_path / "bin"
    tools.mkdir()
    barrier = tmp_path / "barrier"
    barrier.mkdir()
    cp = tools / "cp"
    cp.write_text(
        "#!/bin/bash\nset -euo pipefail\n"
        '/bin/cp "$@"\n'
        'touch "$BARRIER/$BASHPID"\n'
        'for i in $(seq 1 500); do\n'
        '  [ "$(find "$BARRIER" -type f | wc -l)" -ge 2 ] && exit 0\n'
        '  sleep 0.01\n'
        'done\nexit 1\n'
    )
    cp.chmod(0o755)
    env = dict(os.environ, SP_RUN=str(run), BARRIER=str(barrier),
               PATH=f"{tools}:{os.environ['PATH']}")
    children = [subprocess.Popen(
        ["bash", "-euc", prologue], env=env, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True,
    ) for _ in range(2)]
    outputs = [child.communicate(timeout=30) for child in children]
    assert [child.returncode for child in children] == [0, 0], outputs
    staged = list(store.glob("*/wcs/log_exp_headers-*.sqlite"))
    assert len(staged) == 1
    assert staged[0].read_bytes() == source.read_bytes()
    assert not list(store.glob("*/wcs/.log_exp_headers.*"))
