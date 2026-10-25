"""Stage DAG products locally and atomically publish complete generations."""

import argparse
import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

try:
    from . import run_config
except ImportError:
    import run_config
import yaml

MANIFEST = ".allocation-products.json"
COMPLETE = "allocation-complete.json"
OWNER = ".allocation-owner.json"


def store_key(run, state):
    """Separate campaigns sharing a label but not their state."""
    identity = hashlib.sha256(str(Path(state).resolve()).encode()).hexdigest()[
        :16
    ]
    return f"{Path(str(run)).name}-{identity}"


def cpu_budget(environ=None):
    """Read the process cpuset independently of OpenMP thread caps."""
    environ = os.environ if environ is None else environ
    affinity = len(os.sched_getaffinity(0))
    declared = int(environ.get("SLURM_CPUS_PER_TASK") or affinity)
    if min(affinity, declared) <= 0:
        raise ValueError("Cannot derive a positive allocation CPU budget")
    return min(affinity, declared)


def memory_budget(
    cgroup=Path("/sys/fs/cgroup"), proc=Path("/proc/self/cgroup"), environ=None
):
    """Use cgroup memory caps with a SLURM fallback and headroom."""
    environ = os.environ if environ is None else environ
    limits = []
    try:
        relative = next(
            line.split(":", 2)[2].strip()
            for line in proc.read_text().splitlines()
            if line.startswith("0::")
        )
        current = cgroup / relative.lstrip("/")
        while current == cgroup or cgroup in current.parents:
            cap = current / "memory.max"
            if cap.exists() and cap.read_text().strip() not in ("", "max"):
                limits.append(int(cap.read_text()) // 2**20)
            if current == cgroup:
                break
            current = current.parent
    except (OSError, StopIteration, ValueError):
        pass
    if not limits and environ.get("SLURM_MEM_PER_NODE"):
        limits.append(int(environ["SLURM_MEM_PER_NODE"]))
    total = min(limits) if limits else 0
    budget = total - max(4000, total // 20)
    if budget <= 0:
        raise ValueError("Cannot derive a positive allocation memory budget")
    return budget


def bind_stores(profile, local, exposure):
    """Pass ephemeral stores through cleanenv and bind their paths."""
    for path in (local, exposure):
        path.mkdir(parents=True, exist_ok=True)
    arguments = f" --env SP_EXP_STORE={exposure}"
    arguments += "".join(
        f" --bind {path}:{path}" for path in (local, exposure)
    )
    text, count = re.subn(
        r'(apptainer-args:\s*"[^"]*)"',
        lambda match: match.group(1) + arguments + '"',
        profile.read_text(),
        count=1,
    )
    if count != 1:
        raise ValueError(
            "Profile needs quoted apptainer-args for local stores"
        )
    profile.write_text(text)


def entries(root):
    """Yield regular products, refusing internal symlinks."""
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Product symlinks are unsupported: {path}")
        if (
            path.is_file()
            and not path.name.startswith((".sp-publish-", ".sp-record-"))
            and path.relative_to(root).as_posix()
            not in (MANIFEST, COMPLETE, OWNER)
        ):
            yield path


def digest(path):
    """Return the SHA256 digest of one product."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def manifest_entries(root):
    """Read and validate published product paths and metadata."""
    path = root / MANIFEST
    record = json.loads(path.read_text()) if path.exists() else {}
    for name, metadata in record.items():
        relative = Path(name)
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or not relative.parts
        ):
            raise ValueError(f"Invalid published product path: {name!r}")
        if (
            not isinstance(metadata, dict)
            or not isinstance(metadata.get("mtime_ns"), int)
            or not isinstance(metadata.get("sha256"), str)
        ):
            raise ValueError(f"Invalid product metadata: {name!r}")
    return record


def hydrate(source, local):
    """Restore published files and their original timestamp precision."""
    record = manifest_entries(source)
    files = list(entries(source)) if source.exists() else []
    if {path.relative_to(source).as_posix() for path in files} != set(record):
        raise ValueError(
            "Published product tree does not match its integrity manifest"
        )
    local.mkdir(parents=True, exist_ok=True)
    for path in files:
        name = path.relative_to(source).as_posix()
        target = local / name
        if target.exists():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        atomic_copy(path, target)
        if digest(target) != record[name]["sha256"]:
            target.unlink()
            raise ValueError(f"Published product checksum mismatch: {name}")
        stamp = record[name]["mtime_ns"]
        os.utime(target, ns=(stamp, stamp))


def atomic_copy(source, target):
    """Replace one destination file only after its copy is complete."""
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=".sp-publish-", dir=target.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_json(value, target):
    """Atomically replace one JSON metadata file."""
    fd, name = tempfile.mkstemp(prefix=".sp-record-", dir=target.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(value, stream, indent=2)
            stream.write("\n")
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def published_generation(destination):
    """Accept only complete generations owned by this publication root."""
    if not destination.is_symlink():
        if destination.exists():
            raise ValueError(
                "Use a fresh products root; ordinary directories are not "
                f"migrated: {destination}"
            )
        return None
    generation = destination.resolve(strict=True)
    prefix = f".{destination.name}.sp-generation-"
    if (
        generation.parent != destination.parent.resolve()
        or not generation.name.startswith(prefix)
    ):
        raise ValueError(f"Unowned publication symlink: {destination}")
    owner = generation / OWNER
    if not owner.is_file() or json.loads(owner.read_text()) != {
        "destination": str(destination.absolute())
    }:
        raise ValueError(f"Unowned product generation: {generation}")
    if not (generation / COMPLETE).is_file():
        raise ValueError(f"Incomplete published generation: {generation}")
    return generation


def promote(generation, destination):
    """Atomically switch a public root to a complete sibling generation.

    @sc [label:safety] complete-generation-promotion
    Only a fully copied generation with integrity and completion records is
    promoted. tests/unit/test_allocation_products.py injects pre-swap failure.
    """
    fd, name = tempfile.mkstemp(
        prefix=f".{destination.name}.sp-pointer-", dir=destination.parent
    )
    os.close(fd)
    temporary = Path(name)
    temporary.unlink()
    try:
        temporary.symlink_to(generation.name, target_is_directory=True)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def publish(local, destination):
    """Publish a complete generation through atomic symlink replacement."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    published_generation(destination)
    files = list(entries(local))
    generation = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.sp-generation-",
            dir=destination.parent,
        )
    )
    record = {}
    for source in files:
        name = source.relative_to(local).as_posix()
        target = generation / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        record[name] = {
            "mtime_ns": source.stat().st_mtime_ns,
            "sha256": digest(source),
        }
    atomic_json(record, generation / MANIFEST)
    atomic_json(
        {"destination": str(destination.absolute())}, generation / OWNER
    )
    atomic_json(
        {"slurm_job_id": os.environ.get("SLURM_JOB_ID"), "files": len(record)},
        generation / COMPLETE,
    )
    promote(generation, destination)


def prepare(config, run, state, temporary, dry_run=False):
    """Hydrate before planning with one lexical product root across jobs.

    @sc [label:correctness] stable-allocation-product-path
    DAG params use the state/products lexical alias across allocations, never
    the physical local root. The lexical-path fixture guards this boundary.

    @sc [label:safety] no-unpublished-product-reclamation
    Persistent tombstones cannot certify products lost with a timed-out job.
    The launcher fixture exercises dead-disk resume with both cleaners off.
    """
    resolved = run_config.load(str(config), str(run) if run else "")
    destination = Path(
        resolved["outputs"].get("products_dir")
        or resolved["outputs"]["run_dir"]
    )
    if destination.resolve() == Path(resolved["outputs"]["run_dir"]).resolve():
        raise ValueError(
            "In-allocation products need a separate outputs.products_dir"
        )
    alias = state / "products"
    local = temporary / "products" / store_key(resolved["run"], state)
    if alias.exists() and not alias.is_symlink():
        raise ValueError(f"Product alias must be a symlink: {alias}")
    generation = published_generation(destination)
    if generation is None:
        local.mkdir(parents=True, exist_ok=True)
    else:
        hydrate(generation, local)
    pending = state / "products.next"
    pending.unlink(missing_ok=True)
    pending.symlink_to(local, target_is_directory=True)
    os.replace(pending, alias)
    resolved["outputs"]["products_dir"] = str(alias)
    resolved["clean"] = False
    resolved["clean_tiles"] = False
    effective = state / "allocation-config.yaml"
    effective.write_text(yaml.safe_dump(resolved))
    return effective, local, destination


def main():
    """Dispatch the allocation helper command."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("memory")
    sub.add_parser("cpus")
    key = sub.add_parser("key")
    key.add_argument("run")
    key.add_argument("state", type=Path)
    bind = sub.add_parser("bind")
    for name in ("profile", "local", "exposure"):
        bind.add_argument(name, type=Path)
    stage = sub.add_parser("prepare")
    for name in ("config", "run", "state", "temporary"):
        stage.add_argument(f"--{name}", required=name != "run", type=Path)
    stage.add_argument("--dry-run", action="store_true")
    export = sub.add_parser("publish")
    export.add_argument("local", type=Path)
    export.add_argument("destination", type=Path)
    args = parser.parse_args()
    if args.command == "key":
        print(store_key(args.run, args.state))
    elif args.command == "cpus":
        print(cpu_budget())
    elif args.command == "memory":
        print(memory_budget())
    elif args.command == "bind":
        bind_stores(args.profile, args.local, args.exposure)
    elif args.command == "prepare":
        for path in prepare(
            args.config, args.run, args.state, args.temporary, args.dry_run
        ):
            print(path)
    else:
        publish(args.local, args.destination)


if __name__ == "__main__":
    main()
