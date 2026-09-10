"""Property-based state machine over ``workflow/scripts/persist_exp.py``.

``exp_persist`` packs one exposure's keepable PSF products into a tar on
/project and writes a manifest describing it, and its central promise is that
RETENTION IS ADDITIVE: an existing tar is a floor, so shrinking the campaign's
keep list can never delete a product from the backed-up filesystem. That is a
claim about every sequence of keep lists and store states the campaign can
walk through, so it is tested here against a model — random keep lists over a
random set of present products, packed repeatedly, with the tar and the
manifest asserted after every pack.

The keep lists mix product NAMES (``psf_model``), RAW GLOBS (``*.fits``) and
overlapping combinations of the two, because overlap is the case that once
failed every exposure in a campaign: two patterns matching one file is one
file, not a name collision. A genuine collision — two DIFFERENT source paths
landing on one flat member name — must still be fatal, and has its own test.

The script is driven through ``main()`` with a patched ``sys.argv`` rather than
a subprocess: the rule invokes it as a script, but a subprocess per hypothesis
step would put this file out of reach of a login node's time budget.
"""

import fnmatch
import hashlib
import importlib.util
import json
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from hypothesis.stateful import (
    RuleBasedStateMachine,
    initialize,
    precondition,
    rule,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"


def _load(name):
    path = SCRIPTS / f"{name}.py"
    assert path.exists(), f"{path} not found; the rule calls it by path"
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(f"_{name}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


persist = _load("persist_exp")
ALWAYS = persist.ALWAYS

# One concrete file name per catalogued product, in the module output dir the
# real chain writes it to. The names are shaped like the campaign's (module
# tag, exposure, CCD) so the catalogue's globs match them for the same reason
# they match the real thing.
LAYOUT = {
    "star_selection": ("setools", "mask", "star_selection-2079614-5.fits"),
    "star_train": ("setools", "rand_split",
                   "star_split_ratio_80-2079614-5.fits"),
    "star_test": ("setools", "rand_split",
                  "star_split_ratio_20-2079614-5.fits"),
    "star_stats": ("setools", "stat", "star_stat-2079614-5.txt"),
    "psf_model": ("psfex", "", "star_split_ratio_80-2079614-5.psf"),
    "psfex_cat": ("psfex", "", "psfex_cat-2079614-5.cat"),
    "psf_validation": ("psfex_interp", "", "validation_psf-2079614-5.fits"),
}
OPTIONAL = sorted(set(LAYOUT) - {ALWAYS})
# What a campaign can write in `persist_exp:` — names, raw globs, and one name
# the catalogue does not know, which must be refused before any work happens.
ENTRIES = OPTIONAL + ["*.fits", "*.psf", "star_*", "validation_psf-*.fits"]
UNKNOWN = "psf_residuals"

EXP = "2079614"


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def _members(tar: Path) -> list:
    with tarfile.open(tar) as tf:
        return [ti.name for ti in tf.getmembers() if ti.isfile()]


class Store:
    """One exposure's scratch store, its destination, and how to pack it."""

    def __init__(self):
        self.root = Path(tempfile.mkdtemp(prefix="persist-exp-props-"))
        self.exp_dir = self.root / "exp" / EXP
        self.dest = self.root / "products" / "psf"
        self.manifest = self.root / "products" / "manifests" / f"{EXP}.json"
        self.tar = self.dest / f"{EXP}.tar"

    def close(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def path_of(self, product: str) -> Path:
        module, sub, name = LAYOUT[product]
        base = (self.exp_dir / "output" / persist.RUN_NAME
                / f"run_sp_{module}" / "output")
        return (base / sub / name) if sub else (base / name)

    def write(self, product: str, payload: bytes) -> None:
        path = self.path_of(product)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)

    def drop(self, product: str) -> None:
        self.path_of(product).unlink(missing_ok=True)

    def pack(self, keep: list) -> int:
        """Run the script's ``main`` as the rule does. 0 on success."""
        argv = ["persist_exp.py", "--exp-dir", str(self.exp_dir),
                "--exp", EXP, "--dest", str(self.dest),
                "--manifest", str(self.manifest)]
        for entry in keep:
            argv += ["--pattern", entry]
        old = sys.argv
        sys.argv = argv
        try:
            persist.main()
            return 0
        except SystemExit as exc:
            return 1 if exc.code not in (0, None) else 0
        finally:
            sys.argv = old


class PersistExpMachine(RuleBasedStateMachine):
    """A store that gains and loses products under a keep list that changes."""

    @initialize()
    def setup(self):
        self.store = Store()
        self.present = set()
        self.prior_members = []       # members of the last tar written
        self.prior_labels = {}        # member -> product recorded for it

    def teardown(self):
        self.store.close()

    # --- the store and the config move ----------------------------------
    @rule(product=st.sampled_from(sorted(LAYOUT)), size=st.integers(1, 64))
    def add_product(self, product, size):
        self.store.write(product, bytes([len(product) % 251]) * size)
        self.present.add(product)

    @rule(product=st.sampled_from(sorted(LAYOUT)))
    def drop_product(self, product):
        self.store.drop(product)
        self.present.discard(product)

    @rule(keep=st.lists(st.sampled_from(ENTRIES), max_size=4, unique=True))
    def pack(self, keep):
        self._pack_and_check(keep)

    @rule(keep=st.lists(st.sampled_from(ENTRIES), max_size=3, unique=True))
    def pack_with_unknown_product(self, keep):
        """An unknown product name is refused before anything is written."""
        before = (_md5(self.store.tar) if self.store.tar.exists() else None)
        code = self.store.pack(keep + [UNKNOWN])
        assert code != 0, "an unknown product name was accepted"
        after = (_md5(self.store.tar) if self.store.tar.exists() else None)
        assert after == before, "a refused keep list still touched the tar"

    @rule()
    @precondition(lambda self: ALWAYS in self.present)
    def pack_twice_unchanged(self):
        """A rerun over an unchanged store must not move a single byte."""
        keep = sorted(OPTIONAL)[:2]
        self._pack_and_check(keep)
        tar_md5, man_md5 = _md5(self.store.tar), _md5(self.store.manifest)
        tar_mtime = self.store.tar.stat().st_mtime_ns
        man_mtime = self.store.manifest.stat().st_mtime_ns
        assert self.store.pack(keep) == 0
        assert _md5(self.store.tar) == tar_md5, "the tar is not byte-stable"
        assert _md5(self.store.manifest) == man_md5, "the manifest is not byte-stable"
        assert self.store.tar.stat().st_mtime_ns == tar_mtime, (
            "an unchanged rerun rewrote the tar; mtime is a rerun trigger")
        assert self.store.manifest.stat().st_mtime_ns == man_mtime, (
            "an unchanged rerun rewrote the manifest")

    # --- what a pack must leave behind -----------------------------------
    def _pack_and_check(self, keep):
        had_tar = self.store.tar.exists()
        tar_before = _md5(self.store.tar) if had_tar else None
        man_before = (_md5(self.store.manifest)
                      if self.store.manifest.exists() else None)
        code = self.store.pack(keep)

        if ALWAYS not in self.present:
            # The star catalogue's input is not optional: the job fails and
            # nothing downstream may be told the store is safe to reclaim.
            assert code != 0, (
                f"{ALWAYS} is missing and the pack still succeeded")
            assert (_md5(self.store.tar) if self.store.tar.exists()
                    else None) == tar_before, "a failed pack touched the tar"
            assert (_md5(self.store.manifest)
                    if self.store.manifest.exists()
                    else None) == man_before, (
                "a failed pack wrote a manifest; clean_exposure would take "
                "that as permission to delete the store")
            return

        assert code == 0, f"pack failed with {ALWAYS} present and keep={keep}"
        assert self.store.tar.exists() and self.store.manifest.exists()
        members = _members(self.store.tar)
        assert len(members) == len(set(members)), (
            f"duplicate member names in the tar: {members}")

        # ADDITIVE: an existing tar is a floor.
        assert set(members) >= set(self.prior_members), (
            "members vanished from the tar: "
            f"{sorted(set(self.prior_members) - set(members))}")

        body = json.loads(self.store.manifest.read_text())
        listed = {f["name"] for f in body["files"]}
        assert listed == set(members), (
            "manifest and tar disagree about what was packed: "
            f"{sorted(listed ^ set(members))}")
        assert body["n_files"] == len(members)
        assert body["unit"] == EXP and body["status"] == "complete"

        entries = [ALWAYS] + [e for e in keep if e != ALWAYS]
        assert body["products"] == entries
        for f in body["files"]:
            if f["src"] is None:            # carried from the previous tar
                assert f["name"] in self.prior_members
                assert f["product"] == self.prior_labels.get(f["name"], "?")
                continue
            assert f["product"] in entries, (
                f"{f['name']} labelled {f['product']!r}, not in the keep list")
            assert fnmatch.fnmatch(f["name"], persist.resolve(f["product"])), (
                f"{f['name']} does not match {f['product']!r}'s glob")
            assert Path(f["src"]).exists()
            assert f["bytes"] == Path(f["src"]).stat().st_size

        # Every present product the keep list asks for is in there.
        for entry in entries:
            glob = persist.resolve(entry)
            for product in self.present:
                if fnmatch.fnmatch(LAYOUT[product][2], glob):
                    assert LAYOUT[product][2] in listed, (
                        f"{product} matched {entry!r} but was not packed")

        self.prior_members = members
        self.prior_labels = {f["name"]: f["product"] for f in body["files"]}


PersistExpMachine.TestCase.settings = settings(
    max_examples=120,
    stateful_step_count=12,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
TestPersistExpMachine = PersistExpMachine.TestCase


@pytest.fixture()
def store():
    s = Store()
    yield s
    s.close()


def _seed(store, products=(ALWAYS,)):
    for i, product in enumerate(products):
        store.write(product, bytes([i + 1]) * (16 + i))


def test_corrupt_existing_tar_is_refused_and_left_alone(store):
    """A tar that cannot be read may still hold the only copy of something."""
    _seed(store, (ALWAYS, "psf_model"))
    assert store.pack(["psf_model"]) == 0
    store.tar.write_bytes(b"not a tar at all, not even close" * 8)
    corrupt = store.tar.read_bytes()
    man_before = _md5(store.manifest)

    assert store.pack(["psf_model"]) != 0, "a corrupt tar was overwritten"
    assert store.tar.read_bytes() == corrupt, "the corrupt tar was modified"
    assert _md5(store.manifest) == man_before, (
        "a manifest was written over a tar that could not be read")
    assert not store.tar.with_name(store.tar.name + ".tmp").exists()


def test_two_sources_with_one_member_name_is_fatal(store):
    """Members are flat, so a real name clash would silently overwrite."""
    _seed(store, (ALWAYS,))
    # The same file name under a second module output dir.
    clash = (store.exp_dir / "output" / persist.RUN_NAME / "run_sp_setools"
             / "output" / "new_cat" / LAYOUT[ALWAYS][2])
    clash.parent.mkdir(parents=True, exist_ok=True)
    clash.write_bytes(b"a different file with the same name")

    assert store.pack([]) != 0, "two different sources shared a member name"
    assert not store.manifest.exists()
    assert not store.tar.exists()


def test_shrinking_the_keep_list_cannot_delete_a_product(store):
    """The property the additive rule exists for, stated end to end."""
    _seed(store, (ALWAYS, "psf_model", "star_train"))
    assert store.pack(["psf_model", "star_train"]) == 0
    wide = set(_members(store.tar))
    assert LAYOUT["psf_model"][2] in wide

    # The campaign changes its mind, and the scratch store is gone.
    for product in ("psf_model", "star_train"):
        store.drop(product)
    assert store.pack([]) == 0
    assert set(_members(store.tar)) == wide, (
        "shrinking persist_exp: deleted products from the backed-up tar")
    body = json.loads(store.manifest.read_text())
    carried = {f["name"] for f in body["files"] if f["src"] is None}
    assert LAYOUT["psf_model"][2] in carried
    assert {f["name"]: f["product"] for f in body["files"]}[
        LAYOUT["psf_model"][2]] == "psf_model", (
        "a carried member lost the product label the old manifest had")


@settings(max_examples=80, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(st.sampled_from(
    [ALWAYS, "*.fits", "validation_psf-*.fits", "psf_validation",
     "star_*", "*.psf", "psf_model", "star_train"]),
    min_size=1, max_size=5))
def test_overlapping_patterns_never_fail(keep):
    """Two patterns matching one file is one file, not a name collision.

    Overlap is ordinary — ``validation_psf-*.fits`` beside ``*.fits`` is a
    perfectly reasonable way to say "the validation catalogues, and everything
    else FITS while we are here" — and treating the second match as a clash
    once failed every exposure in a campaign.

    A fresh store per example, because the additive rule makes packing
    stateful and this property is about ONE pack.
    """
    s = Store()
    try:
        _seed(s, tuple(LAYOUT))
        assert s.pack(keep) == 0, f"overlapping keep list failed: {keep}"
        members = _members(s.tar)
        assert len(members) == len(set(members)), members
        # Every product present matched something, so all seven are packed.
        assert set(members) == {name for _, _, name in LAYOUT.values()} & set(
            members)
    finally:
        s.close()
