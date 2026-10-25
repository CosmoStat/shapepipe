"""Preserve durable products across publication and hydration failures."""

import fcntl
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "allocation_products", SCRIPTS / "allocation_products.py"
)
products = importlib.util.module_from_spec(spec)
spec.loader.exec_module(products)


def seed(destination, files):
    """Publish fixture files as one owned generation."""
    source = destination.parent / (destination.name + "-seed")
    source.mkdir(parents=True)
    for name, value in files.items():
        target = source / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(value)
    products.publish(source, destination)
    products.shutil.rmtree(source)


@pytest.mark.parametrize("existing", [False, True])
def test_failed_copy_preserves_published_generation(
    tmp_path, monkeypatch, existing
):
    """Keep the public root unchanged when a generation copy fails."""
    local, destination = tmp_path / "local", tmp_path / "published"
    local.mkdir()
    (local / "catalogue").write_text("replacement")
    previous = None
    if existing:
        seed(destination, {"catalogue": "complete original"})
        previous = destination.resolve()

    def interrupted(source, target):
        Path(target).write_text("partial")
        raise OSError("copy interrupted")

    monkeypatch.setattr(products.shutil, "copy2", interrupted)
    with pytest.raises(OSError, match="interrupted"):
        products.publish(local, destination)
    if existing:
        assert destination.resolve() == previous
        assert (destination / "catalogue").read_text() == "complete original"
        assert (destination / products.COMPLETE).exists()
    else:
        assert not destination.exists()
        assert not destination.is_symlink()


def test_hydration_restores_nanoseconds_on_coarse_filesystem(tmp_path):
    """Restore original timestamps after coarse publication storage."""
    local, destination, restored = (
        tmp_path / p for p in ("local", "published", "restored")
    )
    local.mkdir()
    (local / "fragment").write_bytes(b"immutable map fragment")
    stamp = 1700000000123456789
    os.utime(local / "fragment", ns=(stamp, stamp))
    products.publish(local, destination)
    os.utime(destination / "fragment", ns=(1700000000000000000,) * 2)
    products.hydrate(destination, restored)
    assert (restored / "fragment").stat().st_mtime_ns == stamp
    assert (restored / "fragment").read_bytes() == b"immutable map fragment"


def test_interruption_before_pointer_swap_keeps_old_generation(
    tmp_path, monkeypatch
):
    """Keep the old generation selected until promotion completes."""
    local, destination = tmp_path / "local", tmp_path / "published"
    local.mkdir()
    for name in ("a", "b"):
        (local / name).write_text("new")
    seed(destination, {"a": "old", "b": "old"})
    previous = destination.resolve()

    def interrupted(generation, root):
        assert (generation / products.COMPLETE).exists()
        assert (generation / "a").read_text() == "new"
        assert (generation / "b").read_text() == "new"
        raise OSError("before pointer swap")

    monkeypatch.setattr(products, "promote", interrupted)
    with pytest.raises(OSError, match="pointer swap"):
        products.publish(local, destination)
    assert destination.resolve() == previous
    products.hydrate(destination, tmp_path / "restored")
    assert (tmp_path / "restored/a").read_text() == "old"
    assert (tmp_path / "restored/b").read_text() == "old"


def test_publish_refuses_ordinary_directory_and_unowned_root_symlink(tmp_path):
    """Reject roots that the publication helper does not own."""
    local = tmp_path / "local"
    local.mkdir()
    ordinary = tmp_path / "ordinary"
    ordinary.mkdir()
    with pytest.raises(ValueError, match="fresh products"):
        products.publish(local, ordinary)
    link = tmp_path / "link"
    link.symlink_to(ordinary)
    with pytest.raises(ValueError, match="Unowned"):
        products.publish(local, link)


def test_interrupted_hydration_retries_complete_file(tmp_path, monkeypatch):
    """Retry hydration without accepting a truncated target."""
    destination, local = tmp_path / "published", tmp_path / "local"
    seed(destination, {"catalogue": "complete catalogue"})
    original = products.shutil.copy2

    def interrupted(source, target):
        Path(target).write_text("partial")
        raise OSError("hydrate interrupted")

    monkeypatch.setattr(products.shutil, "copy2", interrupted)
    with pytest.raises(OSError, match="hydrate interrupted"):
        products.hydrate(destination, local)
    assert not (local / "catalogue").exists()
    monkeypatch.setattr(products.shutil, "copy2", original)
    products.hydrate(destination, local)
    assert (local / "catalogue").read_text() == "complete catalogue"


def test_hydration_refuses_corrupt_published_bytes(tmp_path):
    """Reject copied bytes that fail their published checksum."""
    destination, local = tmp_path / "published", tmp_path / "local"
    seed(destination, {"catalogue": "complete catalogue"})
    (destination / "catalogue").write_text("corrupt")
    with pytest.raises(ValueError, match="checksum mismatch"):
        products.hydrate(destination, local)
    assert not (local / "catalogue").exists()


@pytest.mark.parametrize("name", ["../outside", "/absolute"])
def test_manifest_cannot_escape_root(tmp_path, name):
    """Reject manifest paths outside the published root."""
    published = tmp_path / "published"
    published.mkdir()
    (published / products.MANIFEST).write_text(
        json.dumps({name: {"mtime_ns": 1, "sha256": "x"}})
    )
    with pytest.raises(ValueError, match="path"):
        products.hydrate(published, tmp_path / "local")


def test_product_symlinks_are_refused_before_publication(tmp_path):
    """Reject internal product links before publishing."""
    local = tmp_path / "local"
    local.mkdir()
    (local / "link").symlink_to(tmp_path / "outside")
    with pytest.raises(ValueError, match="symlinks"):
        products.publish(local, tmp_path / "published")


def test_prepare_keeps_lexical_paths_across_allocations(tmp_path):
    """Keep DAG product paths stable when the physical store changes."""
    state = tmp_path / "state"
    state.mkdir()
    config = tmp_path / "config.yaml"
    destination = tmp_path / "published"
    seed(destination, {"fragment": "durable"})
    (destination / products.COMPLETE).write_text("previous success")
    config.write_text(
        yaml.safe_dump(
            {
                "run": "campaign",
                "outputs": {
                    "run_dir": str(tmp_path / "run"),
                    "products_dir": str(destination),
                },
            }
        )
    )
    first, _, _ = products.prepare(config, None, state, tmp_path / "job1")
    first_path = products.run_config.load(str(config), str(first))["outputs"][
        "products_dir"
    ]
    second, local, _ = products.prepare(
        config, None, state, tmp_path / "job2", dry_run=True
    )
    second_path = products.run_config.load(str(config), str(second))[
        "outputs"
    ]["products_dir"]
    assert first_path == second_path == str(state / "products")
    assert str(local) not in second_path
    assert (state / "products/fragment").read_text() == "durable"
    assert state.joinpath("products").resolve() == local


def test_dry_prepare_preserves_success_marker(tmp_path):
    """Leave published completion evidence unchanged during dry runs."""
    state = tmp_path / "state"
    state.mkdir()
    destination = tmp_path / "published"
    seed(destination, {})
    (destination / products.COMPLETE).write_text("previous success")
    config = tmp_path / "config.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "run": "campaign",
                "outputs": {
                    "run_dir": str(tmp_path / "run"),
                    "products_dir": str(destination),
                },
            }
        )
    )
    products.prepare(config, None, state, tmp_path / "job", dry_run=True)
    assert (destination / products.COMPLETE).read_text() == "previous success"


@pytest.mark.parametrize("failure", [False, True])
def test_launcher_routes_compute_and_publication(tmp_path, failure):
    """Exercise publication, failed compute, resume and launcher locks."""
    tools = tmp_path / "bin"
    tools.mkdir()
    fake = tools / "snakemake"
    fake.write_text("""#!/usr/bin/env python
import os, pathlib, sys, yaml
if os.environ.get("EXPECT_NONALLOC"):
    assert "SP_IN_ALLOCATION" not in os.environ
    sys.exit(0)
assert os.environ.get("SP_IN_ALLOCATION") == "1"
config = yaml.safe_load(pathlib.Path(os.environ["SP_RUN_CONFIG"]).read_text())
assert config["clean"] is False and config["clean_tiles"] is False
with open(os.environ["CALLS"], "a") as stream:
    stream.write(os.environ["SP_PHASE"] + "\\n")
if os.environ["SP_PHASE"] == "compute":
    helper = pathlib.Path(os.environ["SNAPSHOT_PUBLISHER"])
    source = helper.read_text()
    helper.write_text(
        "import os, pathlib, sys\\n"
        "if len(sys.argv) > 1 and sys.argv[1] == 'publish':\\n"
        "    pathlib.Path(os.environ['PUBLISH_MARKER']).touch()\\n"
        + source
    )
    root = pathlib.Path(config["outputs"]["products_dir"])
    (root / "catalogue").write_text("new complete catalogue")
    sys.exit(int(os.environ["FAIL_COMPUTE"]))
""")
    fake.chmod(0o755)
    apptainer = tools / "apptainer"
    apptainer.write_text("#!/bin/sh\nexit 0\n")
    apptainer.chmod(0o755)
    destination = tmp_path / "published"
    seed(destination, {"catalogue": "previous catalogue"})
    (destination / products.COMPLETE).write_text("previous success")
    tiles = tmp_path / "tiles.txt"
    tiles.write_text("123.456\n")
    run = tmp_path / "run.yaml"
    run.write_text(
        yaml.safe_dump(
            {
                "run": "fixture",
                "tile_list": str(tiles),
                "outputs": {
                    "run_dir": str(tmp_path / "run"),
                    "index_db": str(tmp_path / "index.sqlite"),
                    "products_dir": str(destination),
                },
            }
        )
    )
    environment = dict(
        os.environ,
        PATH=str(tools) + os.pathsep + os.environ["PATH"],
        SP_SNAKEMAKE_ENV=str(tmp_path / "no-venv"),
        SP_PROFILE="nibi",
        SLURM_JOB_ID="12345",
        SLURM_TMPDIR=str(tmp_path / "node"),
        SLURM_MEM_PER_NODE="64000",
        CALLS=str(tmp_path / "calls"),
        FAIL_COMPUTE=str(int(failure)),
        SNAPSHOT_PUBLISHER=str(
            tmp_path / "run-state/code/workflow/scripts/allocation_products.py"
        ),
        PUBLISH_MARKER=str(tmp_path / "snapshot-publisher-called"),
    )
    environment.pop("SP_RUN_CONFIG", None)
    environment.pop("SP_STATE_DIR", None)
    result = subprocess.run(
        [
            str(SCRIPTS.parent / "bin/sp"),
            "run",
            "--in-allocation",
            "-c",
            str(run),
        ],
        env=environment,
        capture_output=True,
        text=True,
    )
    assert result.returncode == int(failure), result.stdout + result.stderr
    assert (tmp_path / "calls").read_text().splitlines() == [
        "prepare",
        "compute",
    ]
    marker = tmp_path / "snapshot-publisher-called"
    assert marker.exists() is (not failure)
    assert (destination / "catalogue").read_text() == (
        "previous catalogue" if failure else "new complete catalogue"
    )
    assert (destination / products.COMPLETE).exists()
    assert (
        tmp_path / "run-state/products/catalogue"
    ).read_text() == "new complete catalogue"
    if failure:
        # Model lost allocation storage with durable manifests remaining.
        tile_manifest = (
            tmp_path / "run/tiles/12/123.456/manifests/tile_vignets.json"
        )
        tile_manifest.parent.mkdir(parents=True)
        tile_manifest.write_text("{}")
        exposure = tmp_path / "run/exp/12/12345"
        (exposure / "manifests").mkdir(parents=True)
        (exposure / "manifests/exp_psf.json").write_text("{}")
        (exposure / "output").symlink_to(
            tmp_path / "node/exp-store/fixture/12345"
        )
        products.shutil.rmtree(tmp_path / "node")
        environment.update(FAIL_COMPUTE="0", SLURM_JOB_ID="12346")
        resumed = subprocess.run(
            [
                str(SCRIPTS.parent / "bin/sp"),
                "run",
                "--in-allocation",
                "-c",
                str(run),
            ],
            env=environment,
            capture_output=True,
            text=True,
        )
        assert resumed.returncode == 0, resumed.stdout + resumed.stderr
        assert marker.exists()
        assert not tile_manifest.exists()
        assert not (exposure / "output").is_symlink()
        assert not (exposure / "manifests").exists()
        assert (
            destination / "catalogue"
        ).read_text() == "new complete catalogue"
        assert (destination / products.COMPLETE).exists()
    alias_target = (tmp_path / "run-state/products").readlink()
    with (tmp_path / "run-state/launcher.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        blocked = subprocess.run(
            [
                str(SCRIPTS.parent / "bin/sp"),
                "run",
                "--in-allocation",
                "-c",
                str(run),
            ],
            env=environment,
            capture_output=True,
            text=True,
        )
    assert blocked.returncode == 2
    assert "another launcher" in blocked.stderr
    assert (tmp_path / "run-state/products").readlink() == alias_target
    assert (destination / products.COMPLETE).exists()
    environment.update(EXPECT_NONALLOC="1", SP_IN_ALLOCATION="1")
    ordinary = subprocess.run(
        [str(SCRIPTS.parent / "bin/sp"), "run", "-n", "-c", str(run)],
        env=environment,
        capture_output=True,
        text=True,
    )
    assert ordinary.returncode == 0, ordinary.stdout + ordinary.stderr


def test_memory_budget_uses_step_and_parent_limits(tmp_path):
    """Use the tightest memory cap with scheduler headroom."""
    cgroup = tmp_path / "cgroup"
    step = cgroup / "job/step"
    step.mkdir(parents=True)
    (cgroup / "memory.max").write_text("max")
    (cgroup / "job/memory.max").write_text(str(64000 * 2**20))
    (step / "memory.max").write_text(str(16000 * 2**20))
    proc = tmp_path / "proc"
    proc.write_text("0::/job/step\n")
    assert products.memory_budget(cgroup, proc, {}) == 12000


def test_profile_bind_sources_exist_before_container_launch(tmp_path):
    """Create physical bind sources before the container starts."""
    profile = tmp_path / "profile.yaml"
    profile.write_text('apptainer-args: "--cleanenv"\n')
    local = tmp_path / "node/products/campaign"
    exposure = tmp_path / "node/exp-store/campaign"
    products.bind_stores(profile, local, exposure)
    assert local.is_dir() and exposure.is_dir()
    assert f"--bind {exposure}:{exposure}" in profile.read_text()


def test_cpu_budget_ignores_openmp_and_caps_to_step(monkeypatch):
    """Use CPU affinity and SLURM limits independently of OpenMP."""
    monkeypatch.setattr(
        products.os, "sched_getaffinity", lambda pid: set(range(16))
    )
    assert products.cpu_budget({"OMP_NUM_THREADS": "1"}) == 16
    assert (
        products.cpu_budget(
            {"OMP_NUM_THREADS": "1", "SLURM_CPUS_PER_TASK": "4"}
        )
        == 4
    )


def test_same_label_distinct_campaigns_never_share_local_products(tmp_path):
    """Separate local products for campaigns sharing a run label."""
    local_roots = []
    for name in ("a", "b"):
        state = tmp_path / name / "state"
        state.mkdir(parents=True)
        destination = tmp_path / name / "published"
        seed(destination, {"catalogue": name})
        config = tmp_path / name / "config.yaml"
        config.write_text(
            yaml.safe_dump(
                {
                    "run": "same-label",
                    "outputs": {
                        "run_dir": str(tmp_path / name / "run"),
                        "products_dir": str(destination),
                    },
                }
            )
        )
        _, local, _ = products.prepare(
            config, None, state, tmp_path / "shared-job"
        )
        local_roots.append(local)
        assert (local / "catalogue").read_text() == name
    assert local_roots[0] != local_roots[1]
    assert (local_roots[0] / "catalogue").read_text() == "a"
