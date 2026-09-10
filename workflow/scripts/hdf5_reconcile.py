#!/usr/bin/env python3
"""Bring an hdf5 catalogue into agreement with a campaign, one dataset per unit.

Shared by the two campaign-level merges — ``merge_final_cat.py`` (one dataset
per tile) and ``merge_star_cat.py`` (one per exposure) — because they want
exactly the same thing of their output and disagreeing about it would be a bug
waiting to happen rather than a difference worth having.

WHY RECONCILE RATHER THAN REBUILD. The output must be a function of the input
set — that is what makes the rules' fingerprints mean anything — but reading
every unit to add one is ~800 GB of IO at DR6 scale for a few tens of MB of new
data. So the file is brought INTO AGREEMENT with the campaign instead:

  * a unit with no dataset is read and added;
  * a dataset whose unit has left the campaign is deleted;
  * a dataset whose SOURCE has changed is re-read. Each records its source's
    size and mtime as attributes, and a mismatch is what changed means. This is
    the only reason a finished unit is read twice, and it is why the file
    cannot drift from its inputs the way an append-only tool does;
  * a dataset whose column set was written under a DIFFERENT SCHEMA is re-read.
    The column set is the one input nothing else can see: it is not a source
    file, so no stamp moves when it changes. It travels as a digest on the
    file's root.
  * a dataset that agrees with its source and its schema is left alone, unread.

An append therefore READS exactly the appended units. It still WRITES the whole
file: the existing one is copied so the result can be moved into place
atomically, which costs one pass over it and, briefly, twice its size on disk.
That is the cheap half by orders of magnitude — copying a 1 GB hdf5 against
re-reading 800 GB of catalogues — but it is not free, and `apply` refuses rather
than filling the filesystem when the free space is not there.

WHAT IS AND IS NOT A FUNCTION OF THE INPUT SET. The file's CONTENT is: the same
units with the same sources give the same datasets, the same columns and the
same count attribute, whether they arrived at once or one batch at a time. Its
BYTE LAYOUT is not, because hdf5 lays a group out in the order things were
added. That is the trade for not re-reading the campaign, and it is why the
no-op case compares ACTIONS rather than bytes.

UNTOUCHED ON A NO-OP, which is stronger than byte-stable and cheaper to
establish. Reconciling is PLANNED against a read-only open; an empty plan never
opens the file for writing, so its mtime cannot move — and mtime is a rerun
trigger, so an unconditional rewrite would make every invocation look like a
change.
"""

import hashlib
import shutil
import sys
from pathlib import Path

import h5py


def schema_digest(columns) -> str:
    """A fingerprint of the COLUMN SET the datasets were written with."""
    return hashlib.md5("\n".join(columns).encode()).hexdigest()[:16]


def stamp(path: Path) -> tuple:
    """A source's identity, as recorded on the dataset built from it.

    Size and mtime, not a checksum: the question is "did this change since we
    read it", which mtime answers for a pipeline that writes a file once. A
    campaign that rewrote a source in place with identical size and mtime would
    defeat it, and nothing does.
    """
    st = Path(path).stat()
    return st.st_size, st.st_mtime_ns


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
    """Compare the file on disk with the campaign, WITHOUT writing anything."""
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

    Renaming `campaign:` mid-flight points the rule at a NEW group inside the
    SAME file (the path carries the campaign only on the tile side, where the
    group does). Reconciling would then add a second group beside the first,
    leave the first frozen and stale, and set a count attribute describing only
    one of them. Nothing downstream reads such a file correctly, and no rule
    here means to produce one. Say what is there and stop.
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
            f"only this one. Point `campaign:` back, or write to a new path.")


def apply(output: Path, group_path: str, todo: Plan, units: list, read,
          digest: str, count_attr: str) -> None:
    """Carry the plan out on a tmp file, then move it into place.

    ``read(unit, source)`` returns the structured array for one unit; it is
    called only for the units the plan names, which is what makes an append
    cheap.

    TWO WAYS TO BUILD THE TMP, and which one is used is about SPACE, not speed.
    HDF5 never reclaims the space a deleted dataset occupied, so a file that is
    copied and then edited in place grows for the life of the campaign — every
    refresh of a unit leaks that unit. So:

      * a plan that only ADDS copies the existing file and appends to it. There
        is nothing to reclaim, and copying beats rewriting. It is still a pass
        over the whole file — an append is cheap in READS, not in writes.
      * a plan that removes or refreshes anything builds the tmp FRESH, moving
        the datasets it keeps across with h5py's own group copy — a
        dataset-level copy inside the library that never reads a row into numpy
        — and writing only the units that actually changed. The result is
        compact.

    Either way the tmp is moved into place at the end, so a crash mid-merge
    leaves the old catalogue intact rather than a half-written one. A SIGKILL
    between writing the tmp and renaming it leaves the tmp behind — one file,
    beside the catalogue, overwritten by the next run; the rename itself is
    atomic, which is the property that matters.
    """
    sources = dict(units)
    rewrite = bool(todo.remove or todo.refresh)
    check_free_space(output)
    check_sole_group(output, group_path)
    written = set(todo.add) | set(todo.refresh)
    keep = [u for u, _ in units if u not in written]
    tmp = output.with_name(output.name + ".tmp")
    try:
        tmp.unlink(missing_ok=True)
        if output.exists() and not rewrite:
            shutil.copy2(output, tmp)
        with h5py.File(tmp, "a") as f:
            group = (f[group_path] if group_path in f
                     else f.create_group(group_path))
            if rewrite and output.exists():
                with h5py.File(output, "r") as src:
                    for unit in keep:
                        # File.copy, not Dataset.copy — the latter does not
                        # exist, and the difference only shows when a plan both
                        # rewrites and keeps something.
                        src.copy(f"{group_path}/{unit}", group, name=unit)
            for unit in todo.add + todo.refresh:
                source = sources[unit]
                data = read(unit, source)
                dset = group.create_dataset(unit, data=data, dtype=data.dtype)
                # The dataset's own record of what it was read from; this is
                # what lets a later invocation leave it alone.
                dset.attrs["src_bytes"], dset.attrs["src_mtime_ns"] = \
                    stamp(source)
            f.attrs[count_attr] = len(group)
            f.attrs["param_digest"] = digest
        tmp.replace(output)              # atomic: same filesystem
    finally:
        tmp.unlink(missing_ok=True)
