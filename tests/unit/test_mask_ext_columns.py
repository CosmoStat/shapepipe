"""The MASK_* columns final_cat.param asks for are the ones make_cat writes.

make_cat writes one ``MASK_<label>`` column per ``label:path`` pair in
``MASK_EXT_PATHS`` (config_tile_Mc.ini); the post-processing merge reads
``final_cat.param`` and fails every tile on a name make_cat did not write.
Each label is ``<flag value>_<name>``, and the flag value is the one in the
map's file name: a label on the wrong map would publish one mask under
another's name. Enforces contract ``mask-ext-ladder-columns``.
"""

import configparser
import os
import re
from pathlib import Path

from shapepipe.modules.make_cat_package import make_cat


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "workflow" / "config" / "cfis"


def _mask_ext_paths(monkeypatch):
    monkeypatch.setenv("SP_INPUT_MASKS", "/dummy/masks")
    config = configparser.ConfigParser(interpolation=None)
    config.read(CONFIG_DIR / "config_tile_Mc.ini")
    raw = config.get("MAKE_CAT_RUNNER", "MASK_EXT_PATHS")
    band_paths = make_cat.parse_mask_ext_paths(os.path.expandvars(raw))
    for path in band_paths.values():
        assert path.startswith("/dummy/masks/"), path
    return band_paths


def _written_mask_columns(monkeypatch):
    return [f"MASK_{label}" for label in _mask_ext_paths(monkeypatch)]


def _requested_mask_columns():
    lines = (CONFIG_DIR / "final_cat.param").read_text().splitlines()
    return [line.strip() for line in lines if re.match(r"MASK_", line.strip())]


def test_final_cat_mask_columns_match_mask_ext_paths(monkeypatch):
    """final_cat.param's MASK_* lines are MASK_EXT_PATHS's labels, in order."""
    written = _written_mask_columns(monkeypatch)
    requested = _requested_mask_columns()
    assert len(written) == 11, written
    assert requested == written, (
        "mask-ext-ladder-columns: final_cat.param's MASK_* lines "
        f"{requested} != the columns MASK_EXT_PATHS makes make_cat write "
        f"{written}"
    )


def test_mask_label_flag_value_is_the_maps(monkeypatch):
    """Each label's flag value is the one in its map's file name."""
    for label, path in _mask_ext_paths(monkeypatch).items():
        value, _, name = label.partition("_")
        assert name, f"label {label!r} is not <flag value>_<name>"
        assert path.endswith(f"_n{value}.hsp"), (
            f"mask-ext-ladder-columns: MASK_{label} is read from {path}"
        )
