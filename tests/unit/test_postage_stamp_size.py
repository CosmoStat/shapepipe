"""Configuration invariant for tile and multi-epoch vignette sizes."""

import configparser
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG = REPO_ROOT / "workflow" / "config"
CFIS_CONFIG = CONFIG / "cfis"


def _ini_int(path, section, option):
    parser = configparser.ConfigParser(interpolation=None)
    assert parser.read(path) == [str(path)]
    return parser.getint(section, option)


def _multi_epoch_stamp_size(path):
    """STAMP_SIZE of the one vignetmaker section run in MULTI-EPOCH mode."""
    parser = configparser.ConfigParser(interpolation=None)
    assert parser.read(path) == [str(path)]
    sections = [
        section
        for section in parser.sections()
        if section.startswith("VIGNETMAKER_RUNNER")
        and parser.get(section, "MODE", fallback="") == "MULTI-EPOCH"
    ]
    assert len(sections) == 1, f"{path.name}: multi-epoch sections {sections}"
    return sections[0], parser.getint(sections[0], "STAMP_SIZE")


def _param_vignet_size(path):
    matches = []
    for line in path.read_text(encoding="utf-8").splitlines():
        setting = line.split("#", 1)[0].strip()
        match = re.match(r"^VIGNET\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)", setting)
        if match:
            width, height = map(int, match.groups())
            assert width == height, f"{path.name}: VIGNET must be square"
            matches.append(width)

    assert len(matches) == 1, f"{path.name}: expected one VIGNET setting"
    return matches[0]


@pytest.mark.decision("postage_stamp_size")
def test_tile_and_epoch_vignet_sizes_match():
    """The tile markers align pixel-for-pixel with each multi-epoch stamp."""
    param_file = CFIS_CONFIG / "default_noimaflags.param"
    sizes = {
        "config_tile_Uc.ini#READ_EXT_SEXCAT_RUNNER.VIGNET_SIZE": _ini_int(
            CFIS_CONFIG / "config_tile_Uc.ini",
            "READ_EXT_SEXCAT_RUNNER",
            "VIGNET_SIZE",
        ),
        f"{param_file.name}#VIGNET": _param_vignet_size(param_file),
    }

    vignet_configs = {
        path.resolve()
        for flavour in ("cfis", "cfis_image_sims")
        for path in (CONFIG / flavour).glob("config_tile_PiViVi_*.ini")
    }
    assert vignet_configs, "No config_tile_PiViVi_*.ini files found"
    for path in sorted(vignet_configs):
        section, size = _multi_epoch_stamp_size(path)
        sizes[f"{path.parent.name}/{path.name}#{section}.STAMP_SIZE"] = size

    assert len(set(sizes.values())) == 1, f"Vignette sizes differ: {sizes}"
