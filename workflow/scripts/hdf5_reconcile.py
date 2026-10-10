#!/usr/bin/env python3
"""Reconcile an HDF5 catalogue with a campaign, one dataset per unit.

``merge_final_cat.py`` uses one dataset per tile; ``merge_star_cat.py`` uses one
per exposure. Reconciliation reads only added or changed units, avoiding about
800 GB of source reads at DR6 scale to append a few tens of MB.

``plan`` adds missing datasets, removes units outside the campaign, and refreshes
datasets whose source size/mtime or root schema digest differs. The digest
covers column names; types are checked when units are read. A conflict with
kept datasets triggers a full refresh; conflicting source types then fail.

@sc [label:schema] hdf5-reconcile-uniform-column-types
All units must agree on each column's kind, itemsize and shape, independent of
byte order. Refuse conflicting sources rather than silently promote values.
A name digest alone cannot enforce this constraint.

Publication copies or rebuilds the whole HDF5 file beside the target before
atomic replacement, requiring temporary disk space even for an append. An
add-only plan copies the existing file; refreshes and removals rebuild it to
avoid retaining space occupied by deleted datasets. ``check_free_space`` applies
a conservative margin before rewriting an existing catalogue.

The same units and sources yield the same datasets, columns and count attribute,
regardless of arrival order. HDF5 byte layout can differ with insertion order,
so no-op detection uses the plan rather than a byte comparison.

@sc [label:operations] hdf5-reconcile-no-op-mtime
Plan against a read-only open and have callers skip ``apply`` for an empty plan.
This preserves mtime and avoids triggering downstream reruns. ``apply`` itself
does not short-circuit an empty plan.
"""

import hashlib
import json
import shutil
import sys
from pathlib import Path

import h5py


# Every unit's dataset is written lzf-compressed: about half the bytes of an
# uncompressed catalogue, at no measurable write cost. lzf ships with h5py, so
# any h5py reader decodes it transparently; plain libhdf5 tools (h5dump, C
# readers) do not carry the filter. Kept datasets move across with the
# library's group copy, which preserves their filter as written.
COMPRESSION = "lzf"


def schema_digest(columns) -> str:
    """Fingerprint the supplied column names in their supplied order."""
    return hashlib.md5("\n".join(columns).encode()).hexdigest()[:16]


def code_provenance(snapshot_json) -> dict:
    """The launch code's identity, to stamp onto the merged file's root.

    ``snapshot_json`` is ``sp run``'s code snapshot (``bin/sp``'s
    ``$STATE_DIR/code/snapshot.json``), passed through by the calling rule. A
    workflow driven outside ``sp run`` has no such file — the merge still
    succeeds, and the caller writes ``code_head = "unknown"`` rather than
    failing an otherwise-good build.
    """
    if not snapshot_json or not Path(snapshot_json).exists():
        return {"head": "unknown"}
    data = json.loads(Path(snapshot_json).read_text())
    out = {k: data[k] for k in ("head", "branch", "dirty", "taken_at")
           if k in data}
    if data.get("dirty") and data.get("dirty_files"):
        out["dirty_files"] = data["dirty_files"]
    return out


def stamp(path: Path) -> tuple:
    """A source's identity, as recorded on the dataset built from it.

    @sc [label:coupling] hdf5-reconcile-source-stamp
    Source changes must alter size or nanosecond mtime. This stamp is not a
    checksum: rewriting a source with identical size and mtime leaves its
    dataset indistinguishable from an unchanged one.
    """
    st = Path(path).stat()
    return st.st_size, st.st_mtime_ns


def column_types(dtype) -> dict:
    """``{column: (kind, itemsize, shape)}``: a dataset's schema, as compared
    across units. Byte order is left out: a source may store either, and
    neither changes a value."""
    return {n: (dtype[n].base.kind, dtype[n].base.itemsize, dtype[n].shape)
            for n in dtype.names}


def type_conflict(a_unit, a_dtype, b_unit, b_dtype) -> str:
    """Name the first column whose type differs between two units' dtypes."""
    a, b = column_types(a_dtype), column_types(b_dtype)
    for col in a:
        if a[col] != b.get(col):
            got = b_dtype[col].base.name if col in b else "absent"
            return (f"column {col} is {a_dtype[col].base.name} in {a_unit} "
                    f"but {got} in {b_unit}")
    return f"{b_unit} carries columns {a_unit} does not"


class _Retyped(Exception):
    """A read unit's types differ from the datasets `apply` would keep."""


class Plan:
    """What reconciling requires: three unit lists.

    ``add`` and ``refresh`` are both "read the source and write the dataset";
    they are separate only so the log can say which happened, because a refresh
    means a finished unit's source moved under us and that is worth seeing.
    """

    def __init__(self, add, refresh, remove):
        self.add, self.refresh, self.remove = add, refresh, remove

    def empty(self):
        return not (self.add or self.refresh or self.remove)

    def describe(self):
        return (f"{len(self.add)} added, {len(self.refresh)} refreshed, "
                f"{len(self.remove)} removed")


def plan(output: Path, group_path: str, units: list, digest: str) -> Plan:
    """Compare the file on disk with the campaign without writing anything."""
    if not output.exists():
        return Plan([u for u, _ in units], [], [])

    want = {unit for unit, _ in units}
    add, refresh = [], []
    with h5py.File(output, "r") as f:
        stale_schema = f.attrs.get("param_digest") != digest
        have = dict(f[group_path].items()) if group_path in f else {}
        present = set(have)
        for unit, source in units:
            if unit not in present:
                add.append(unit)
            elif stale_schema:
                refresh.append(unit)
            else:
                attrs = have[unit].attrs
                if (int(attrs.get("src_bytes", -1)),
                        int(attrs.get("src_mtime_ns", -1))) != stamp(source):
                    refresh.append(unit)
    return Plan(add, refresh, sorted(present - want))


# Twice the file, plus a tenth of it again: the copy and the original coexist,
# and hdf5 is not a format to run to the last byte of a filesystem on.
FREE_SPACE_MARGIN = 2.1


def check_free_space(output: Path) -> None:
    """Refuse to start a rewrite the filesystem cannot hold.

    A merge that fills /project does not just fail: it fails everything else
    writing there at the same time, and it can leave a truncated tmp beside a
    catalogue people trust. Cheaper to say so first.
    """
    if not output.exists():
        return
    size = output.stat().st_size
    free = shutil.disk_usage(output.parent).free
    if free < size * FREE_SPACE_MARGIN:
        sys.exit(
            f"hdf5_reconcile: {output.parent} has {free / 1e9:.1f} GB free and "
            f"this merge needs about {size * FREE_SPACE_MARGIN / 1e9:.1f} GB — "
            f"it rewrites {output.name} ({size / 1e9:.1f} GB) through a tmp "
            f"copy beside it. Free space or move products_dir; the existing "
            f"catalogue is untouched.")


def check_sole_group(output: Path, group_path: str) -> None:
    """One file, one campaign — refuse to half-update a file holding two.

    @sc [label:schema] hdf5-reconcile-single-campaign-group
    No sibling campaign group may exist under the requested parent. Reconciling
    only one group would leave the other stale while the root count describes
    just the updated group. Use a separate output path for each campaign.
    """
    if not output.exists() or "/" not in group_path:
        return
    parent, leaf = group_path.rsplit("/", 1)
    with h5py.File(output, "r") as f:
        if parent not in f:
            return
        others = sorted(k for k in f[parent] if k != leaf)
    if others:
        sys.exit(
            f"hdf5_reconcile: {output} already holds {parent}/"
            f"{', '.join(others)} beside {group_path}. One file is one "
            f"campaign: reconciling would freeze the other group and count "
            f"only this one. Point `run:` back, or write to a new path.")


def apply(output: Path, group_path: str, todo: Plan, units: list, read,
          digest: str, count_attr: str, provenance: dict | None = None) -> None:
    """Carry the plan out; re-read every unit if the column types changed.

    See ``_apply``. When a unit read under the plan has other column types
    than the datasets the plan would keep, the plan widens to refresh every
    kept unit, so the file keeps one type per column (module docstring).
    """
    try:
        _apply(output, group_path, todo, units, read, digest, count_attr,
               provenance)
    except _Retyped as exc:
        every = Plan(todo.add, [u for u, _ in units if u not in todo.add],
                     todo.remove)
        print(f"[hdf5_reconcile] {exc}; re-reading all {len(units)} unit(s)")
        _apply(output, group_path, every, units, read, digest, count_attr,
               provenance)


def _apply(output: Path, group_path: str, todo: Plan, units: list, read,
           digest: str, count_attr: str, provenance: dict | None) -> None:
    """Carry the plan out on a tmp file, then move it into place.

    ``read(unit, source)`` returns the structured array for one unit; it is
    called only for the units the plan names, which is what makes an append
    cheap.

    ``provenance`` (``code_provenance()``'s return) is stamped onto the file's
    root as ``code_head``/``code_branch``/``code_dirty``/``code_snapshot_at``,
    plus ``code_dirty_files`` (newline-joined) when the snapshot was dirty. It
    is written here, alongside ``count_attr`` and ``param_digest``, rather than
    on every no-op invocation: reconciling is planned against a read-only open,
    and an empty plan must leave the file's mtime alone (see the module
    docstring), so a run that changes no data never touches the file even if
    the code that would have produced it has moved on.

    See the module docstring for copy-versus-rebuild selection. Kept datasets
    move through h5py's group copy without loading rows into NumPy.

    @sc [label:custody] hdf5-reconcile-atomic-publication
    Write beside the output and replace it only after the merge succeeds.
    A crash before replacement leaves the published catalogue intact. SIGKILL
    can leave the tmp behind; the next writer removes it before checking space.
    One writer per output is required.
    """
    sources = dict(units)
    rewrite = bool(todo.remove or todo.refresh)
    tmp = output.with_name(output.name + ".tmp")
    # A tmp left by a killed run is this run's to delete (one writer per
    # output), and deleting it first keeps it from counting against the space
    # this run needs.
    tmp.unlink(missing_ok=True)
    check_free_space(output)
    check_sole_group(output, group_path)
    written = set(todo.add) | set(todo.refresh)
    keep = [u for u, _ in units if u not in written]
    try:
        if output.exists() and not rewrite:
            shutil.copy2(output, tmp)
        with h5py.File(tmp, "a") as f:
            group = (f[group_path] if group_path in f
                     else f.create_group(group_path))
            if rewrite and output.exists():
                with h5py.File(output, "r") as src:
                    for unit in keep:
                        # Copy kept datasets through the owning HDF5 file.
                        src.copy(f"{group_path}/{unit}", group, name=unit)
            # The kept datasets' types are the reference a read unit must
            # match; with nothing kept, the first unit read is the reference.
            # Conflicts among kept datasets also trigger a full refresh.
            kept = [(u, group[u].dtype) for u in keep if u in group]
            ref = kept[0] if kept else None
            for unit, dtype in kept[1:]:
                if column_types(dtype) != column_types(ref[1]):
                    raise _Retyped(type_conflict(*ref, unit, dtype))
            for unit in todo.add + todo.refresh:
                source = sources[unit]
                data = read(unit, source)
                if ref is None:
                    ref = (unit, data.dtype)
                elif column_types(data.dtype) != column_types(ref[1]):
                    conflict = type_conflict(*ref, unit, data.dtype)
                    if keep:
                        raise _Retyped(conflict)
                    sys.exit(f"hdf5_reconcile: {conflict}. One catalogue "
                             f"holds one type per column; remake the units "
                             f"whose sources are stale. {output} is "
                             f"untouched.")
                dset = group.create_dataset(unit, data=data, dtype=data.dtype,
                                            compression=COMPRESSION)
                # The dataset's own record of what it was read from; this is
                # what lets a later invocation leave it alone.
                dset.attrs["src_bytes"], dset.attrs["src_mtime_ns"] = \
                    stamp(source)
            f.attrs[count_attr] = len(group)
            f.attrs["param_digest"] = digest
            if provenance:
                # One record at a time: an add-only merge copied the last
                # one's attributes, and a snapshot-less record must not keep
                # its branch, dirty flag, dirty files or time.
                for attr in [a for a in f.attrs if a.startswith("code_")]:
                    del f.attrs[attr]
                f.attrs["code_head"] = provenance.get("head", "unknown")
                for key, attr in (("branch", "code_branch"),
                                  ("dirty", "code_dirty"),
                                  ("taken_at", "code_snapshot_at")):
                    if key in provenance:
                        f.attrs[attr] = provenance[key]
                if provenance.get("dirty_files"):
                    f.attrs["code_dirty_files"] = \
                        "\n".join(provenance["dirty_files"])
        tmp.replace(output)              # atomic: same filesystem
    finally:
        tmp.unlink(missing_ok=True)
