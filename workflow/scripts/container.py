#!/usr/bin/env python3
"""Manage this user's copy of the ShapePipe container image.

Two layers, and the second only exists when you ask for one:

* the **SIF** (``~/.cache/shapepipe/shapepipe.sif``) -- a pristine, read-only
  copy of the published image, pulled into your own cache. Per-user by
  construction: one file, one owner, nobody else's refresh moves the ground
  under a running job.
* an optional **sandbox** (``~/.cache/shapepipe/sandbox/``) -- the same image
  unpacked into a writable directory, so a ``pip install`` inside it sticks.
  The escape hatch for work that needs a package the image does not carry yet.

Resolution order, shared by this CLI and by the workflow: **sandbox if it
exists, else the cached SIF if it exists, else the ``container:`` path in
workflow/config.yaml**. That last one is the current shared /project image, so
a checkout with an empty cache behaves exactly as it did before this verb
existed, and a package installed into your sandbox is there for your workflow
jobs too.

Subcommands, exposed as ``sp container <verb>``::

    sp container pull                     # fetch the tag into the cache
    sp container status                   # what is here, and how current it is
    sp container sandbox                  # unpack the SIF into a writable dir
    sp container exec <cmd...>            # run something inside it
    sp container exec --writable <cmd...> # ... with writes that persist

``pull`` needs the network. Compute nodes on Alliance clusters generally have
none, so run it on a login node or inside an ``salloc`` allocation -- never
from a batch job.

Deliberately **stdlib-only**: it runs on the bare host, outside the container,
where the science stack is not installed, and so must import without it.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# The published image. CI pushes one tag per branch, sanitized; `-runtime` is
# the slim variant the workflow runs.
CONTAINER_URI = "docker://ghcr.io/cosmostat/shapepipe:develop-runtime"

# The single source of truth for the fallback image: the workflow's own
# `container:` key, which is also what the Snakefile reads. Written down once,
# here, so the CLI and the workflow cannot disagree about the default.
CONFIG_FILE = Path(__file__).resolve().parents[1] / "config.yaml"
CONFIG_KEY = "container"

# ~/.cache/shapepipe by default; SP_CACHE_DIR moves the whole cache (e.g. onto
# a filesystem with room), XDG_CACHE_HOME moves it with everything else.
CACHE_DIR = Path(
    os.environ.get("SP_CACHE_DIR")
    or Path(os.environ.get("XDG_CACHE_HOME", "~/.cache")) / "shapepipe"
).expanduser()

# This user's pristine image. Override with ``SP_CONTAINER`` (absolute path).
DEFAULT_SIF = CACHE_DIR / "shapepipe.sif"

# The optional writable unpacking of it. Override with ``SP_SANDBOX``.
DEFAULT_SANDBOX = CACHE_DIR / "sandbox"

# Bind mounts for `exec`, matching the nibi profile's apptainer-args (the two
# cluster filesystems this workflow reads and writes, plus the home that holds
# ~/.ssl/cadcproxy.pem). Override wholesale with ``SP_APPTAINER_BINDS``.
DEFAULT_BINDS = "/project,/scratch,/home"


def configured_default():
    """Return the ``container:`` path from workflow/config.yaml, or ``None``.

    A deliberately minimal scalar read (the same one bin/sp does in sed): this
    module is stdlib-only, so there is no yaml to import.
    """
    try:
        text = CONFIG_FILE.read_text()
    except OSError:
        return None
    match = re.search(rf"^{CONFIG_KEY}:[ \t]*(\S+)", text, re.MULTILINE)
    return match.group(1) if match else None


def local_sif():
    """Return this user's cached image path (may not exist yet)."""
    override = os.environ.get("SP_CONTAINER")
    return (Path(override) if override else DEFAULT_SIF).expanduser()


def local_sandbox():
    """Return this user's writable sandbox directory (may not exist)."""
    override = os.environ.get("SP_SANDBOX")
    return (Path(override) if override else DEFAULT_SANDBOX).expanduser()


def resolve_image():
    """Return ``(path, kind)`` for the image everything should run.

    ``kind`` is ``"sandbox"``, ``"sif"``, ``"configured"`` (the shared
    /project image named in config.yaml -- the default when the cache is
    empty) or ``"none"``.
    """
    sandbox = local_sandbox()
    if sandbox.is_dir():
        return str(sandbox), "sandbox"
    sif = local_sif()
    if sif.exists():
        return str(sif), "sif"
    default = configured_default()
    if default:
        return default, "configured"
    return "", "none"


def image_labels(image):
    """Return the image's OCI labels as a dict, or ``{}`` if unreadable.

    Never raises: a missing file, a missing ``apptainer`` or a corrupt image
    all mean "we don't know", which every caller treats as non-fatal.
    """
    path = Path(image)
    if not path.exists() or shutil.which("apptainer") is None:
        return {}
    try:
        out = subprocess.run(
            ["apptainer", "inspect", "--labels", str(path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return {}
    if out.returncode != 0:
        return {}
    labels = {}
    for line in out.stdout.splitlines():
        key, sep, value = line.partition(":")
        if sep:
            labels[key.strip()] = value.strip()
    return labels


def image_revision(image):
    """Return the commit the image was built from, or ``None``."""
    return image_labels(image).get("org.opencontainers.image.revision")


def _require_apptainer():
    """Exit unless ``apptainer`` is on PATH (bin/sp loads the module)."""
    if shutil.which("apptainer") is None:
        sys.exit("apptainer is not on PATH (bin/sp loads apptainer/1.4.5)")


def _git(*args, cwd=None):
    """Run git, returning stripped stdout or ``None`` on any failure."""
    try:
        out = subprocess.run(
            ["git", *args], capture_output=True, text=True, cwd=cwd, timeout=30
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def compare_revision(revision, repo=None):
    """Place an image revision relative to this checkout's HEAD.

    One of ``"in-sync"``, ``"behind"`` (the image predates HEAD), ``"ahead"``,
    ``"diverged"``, or ``"unknown"`` (no label, no git, or a commit this clone
    has never fetched).
    """
    if not revision:
        return "unknown"
    repo = repo or Path(__file__).resolve().parents[2]
    head = _git("rev-parse", "HEAD", cwd=repo)
    if head is None:
        return "unknown"
    if head == revision:
        return "in-sync"
    if _git("cat-file", "-e", f"{revision}^{{commit}}", cwd=repo) is None:
        return "unknown"
    if _git("merge-base", "--is-ancestor", revision, head, cwd=repo) is not None:
        return "behind"
    if _git("merge-base", "--is-ancestor", head, revision, cwd=repo) is not None:
        return "ahead"
    return "diverged"


def cmd_pull(args):
    """Pull ``--tag`` into the cache, atomically."""
    _require_apptainer()
    sif = local_sif()
    sif.parent.mkdir(parents=True, exist_ok=True)
    # Pull to a sibling temp name and rename: an atomic rename within one
    # directory, so an in-flight job sees either the whole old image or the
    # whole new one. Pulling in place leaves the file half-written for the many
    # minutes the pull takes. Jobs already running hold the old inode open.
    tmp = sif.with_name(sif.name + f".pull.{os.getpid()}")
    print(f"pulling {args.tag}\n     -> {sif}")
    try:
        subprocess.run(
            ["apptainer", "pull", "--force", "--name", str(tmp), args.tag], check=True
        )
        os.replace(tmp, sif)
    except subprocess.CalledProcessError as exc:
        tmp.unlink(missing_ok=True)
        sys.exit(f"pull failed ({exc.returncode}); {sif} is unchanged")
    except KeyboardInterrupt:
        tmp.unlink(missing_ok=True)
        raise
    labels = image_labels(sif)
    print(f"revision: {labels.get('org.opencontainers.image.revision', 'unknown')}")
    print(f"version:  {labels.get('org.opencontainers.image.version', 'unknown')}")
    return 0


def cmd_sandbox(args):
    """Unpack the image into a writable directory -- the opt-in escape hatch."""
    _require_apptainer()
    sandbox = local_sandbox()
    if sandbox.exists() and not args.force:
        sys.exit(
            f"sandbox already exists at {sandbox}\n"
            "pass --force to discard it and rebuild from a clean image"
        )
    source = args.source
    if not source:
        image, kind = resolve_image()
        if kind == "sandbox":
            # Rebuilding from the sandbox itself would just re-copy the drift.
            image = str(local_sif()) if local_sif().exists() else (
                configured_default() or CONTAINER_URI
            )
        source = image or CONTAINER_URI
    sandbox.parent.mkdir(parents=True, exist_ok=True)
    print(f"building sandbox from {source}\n     -> {sandbox}")
    # Build beside the target and swap it in, as `pull` does -- and for a
    # sharper reason. A half-written .sif fails loudly, but a half-unpacked
    # sandbox *directory* is still a directory, so resolve_image() would elect
    # it and every job would silently run a broken tree. Staging also means a
    # --force rebuild that fails leaves the sandbox you already had intact.
    #
    # `--fix-perms` so the tree can be deleted again later. No `--fakeroot`: an
    # unprivileged build from an existing image goes through user namespaces,
    # which is what the Alliance clusters provide.
    staging = sandbox.with_name(f"{sandbox.name}.build.{os.getpid()}")
    shutil.rmtree(staging, ignore_errors=True)
    try:
        subprocess.run(
            ["apptainer", "build", "--sandbox", "--fix-perms", str(staging), source],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        shutil.rmtree(staging, ignore_errors=True)
        sys.exit(f"sandbox build failed ({exc.returncode}); {sandbox} is unchanged")
    except (KeyboardInterrupt, OSError):
        shutil.rmtree(staging, ignore_errors=True)
        raise

    if sandbox.exists():
        print(f"replacing {sandbox}")
        shutil.rmtree(sandbox, ignore_errors=True)
        if sandbox.exists():
            shutil.rmtree(staging, ignore_errors=True)
            sys.exit(f"could not remove {sandbox}; remove it by hand and retry")
    os.replace(staging, sandbox)
    print(
        "\nthis sandbox now takes precedence over the SIF everywhere, including "
        "workflow jobs.\ninstall into it with: sp container exec --writable pip "
        "install <pkg>\nreset to a clean image with: sp container pull && "
        "sp container sandbox --force"
    )
    return 0


def cmd_status(args):
    """Report which image layer is live, its revision, and how current it is."""
    sif = local_sif()
    sandbox = local_sandbox()
    active, kind = resolve_image()

    if sif.exists():
        print(f"SIF:       {sif} ({sif.stat().st_size / 1e9:.1f} GB)")
    else:
        print(f"SIF:       absent ({sif})")
    if sandbox.is_dir():
        print(f"sandbox:   {sandbox} (writable; may carry local modifications)")
    else:
        print("sandbox:   none")
    print(f"configured: {configured_default() or 'unset'} ({CONFIG_FILE})")

    if kind == "none":
        print("\nactive:    NONE -- no cached image and no container: in config.yaml")
        print(f"run: sp container pull --tag {CONTAINER_URI}")
        return 1

    print(f"\nactive:    {active} ({kind})")
    if kind == "configured" and not Path(active).exists():
        print("           WARNING: that path does not exist on this host")
        return 1

    labels = image_labels(active)
    revision = labels.get("org.opencontainers.image.revision")
    source = ""
    if revision is None and kind == "sandbox" and sif.exists():
        # Some sandbox trees do not carry the original labels through. The SIF
        # beside it is the best remaining evidence of what it was built from --
        # a guess, so it is labelled as one rather than printed as fact.
        revision = image_revision(sif)
        if revision:
            source = " (inferred from the SIF beside it, not read from the sandbox)"
    print(f"revision:  {revision or 'unknown'}{source}")
    print(f"version:   {labels.get('org.opencontainers.image.version', 'unknown')}")
    if kind == "sandbox":
        print(
            "           (that revision is what the sandbox was built from; "
            "anything\n            installed into it since is in no label)"
        )
    verdict = compare_revision(revision)
    explain = {
        "in-sync": "matches this checkout's HEAD",
        "behind": "older than this checkout's HEAD -- pull to refresh",
        "ahead": "newer than this checkout's HEAD",
        "diverged": "on a different branch from this checkout",
        "unknown": "cannot compare (no label, or a commit this clone lacks)",
    }[verdict]
    print(f"checkout:  {verdict} ({explain})")
    return 0


def cmd_exec(args):
    """Run a command inside the resolved image -- the one-off path."""
    _require_apptainer()
    command = [a for a in args.command if a != "--"]
    if not command:
        sys.exit("nothing to run; pass a command after `exec`")
    binds = args.bind or os.environ.get("SP_APPTAINER_BINDS", DEFAULT_BINDS)

    if args.writable:
        # A SIF is a read-only filesystem, so `--writable` against one fails
        # obscurely; only a sandbox takes writes.
        sandbox = local_sandbox()
        if not sandbox.is_dir():
            sys.exit(
                f"--writable needs a sandbox, and there is none at {sandbox}\n"
                "build one with: sp container sandbox"
            )
        image, extra = str(sandbox), ["--writable"]
    else:
        image, kind = resolve_image()
        if kind == "none":
            sys.exit("no image resolved; run: sp container pull")
        extra = []

    cmd = ["apptainer", "exec", *extra, "--cleanenv", "--bind", binds, image, *command]
    return subprocess.run(cmd).returncode


def cmd_resolve(args):
    """Print just the resolved image path -- what the Snakefile consumes."""
    image, kind = resolve_image()
    if kind == "none":
        sys.exit("no image resolved; run: sp container pull")
    print(image)
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        prog="sp container",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="subcommand", required=True)

    p_pull = sub.add_parser(
        "pull",
        help="fetch the image into the cache (login node or salloc: compute "
        "nodes have no network)",
    )
    p_pull.add_argument(
        "--tag", default=CONTAINER_URI, help=f"image to pull (default: {CONTAINER_URI})"
    )
    p_pull.set_defaults(func=cmd_pull)

    p_status = sub.add_parser(
        "status", help="report which image layer is live and how current it is"
    )
    p_status.set_defaults(func=cmd_status)

    p_sandbox = sub.add_parser(
        "sandbox", help="unpack the image into a writable directory (opt-in)"
    )
    p_sandbox.add_argument(
        "--source", help="image to unpack (default: the cached SIF, or the config path)"
    )
    p_sandbox.add_argument(
        "--force",
        action="store_true",
        help="discard an existing sandbox and rebuild from a clean image",
    )
    p_sandbox.set_defaults(func=cmd_sandbox)

    p_exec = sub.add_parser("exec", help="run a command inside the resolved image")
    p_exec.add_argument("--bind", help=f"bind mounts (default: {DEFAULT_BINDS})")
    p_exec.add_argument(
        "--writable",
        action="store_true",
        help="run against the sandbox so writes (e.g. pip install) persist",
    )
    p_exec.add_argument("command", nargs=argparse.REMAINDER)
    p_exec.set_defaults(func=cmd_exec)

    p_resolve = sub.add_parser(
        "resolve", help="print the resolved image path (what the workflow runs)"
    )
    p_resolve.set_defaults(func=cmd_resolve)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
