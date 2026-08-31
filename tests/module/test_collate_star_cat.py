"""UNIT TESTS FOR STAR-CATALOGUE COLLATION PATHS.

Pin the patch vs patch-less (v2.0) path and filename convention of
``scripts/python/collate_star_cat.py``. Runs up to v1.6 carry a ``P<patch>``
token in both the input run directory and the output filename; v2.0 is
patch-less (``patch is None``) and drops that token, reading from a single
``<base>/output`` root and writing ``validation_psf_conv-<idx>.fits`` — the
name still matched by the downstream ``validation_psf_conv-*`` glob.
"""

import importlib.util
from pathlib import Path

import pytest

# The collation script lives under scripts/python (not an importable package),
# so load it by path.
_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "python"
    / "collate_star_cat.py"
)
_spec = importlib.util.spec_from_file_location("collate_star_cat", _SCRIPT)
collate_star_cat = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(collate_star_cat)


@pytest.mark.parametrize(
    "patch, exp_input, exp_output",
    [
        ("3", "in/P3/output/", "out/P3"),
        (None, "in/output/", "out"),
    ],
)
def test_collate_paths(patch, exp_input, exp_output):
    """v1.x carries the P<patch> token; v2.0 (patch None) drops it."""
    assert collate_star_cat.collate_paths("in", "out", patch) == (
        exp_input,
        exp_output,
    )


@pytest.mark.parametrize(
    "patch, expected",
    [
        ("3", "validation_psf_conv-3-0.fits"),
        (None, "validation_psf_conv-0.fits"),
    ],
)
def test_output_filename(patch, expected):
    """The patch token is present for v1.x and absent for v2.0."""
    assert collate_star_cat.output_filename("validation_psf", patch, 0) == expected


def test_output_filename_matches_downstream_glob():
    """Both layouts stay under the downstream ``validation_psf_conv-*`` glob."""
    for patch in ("1", None):
        assert collate_star_cat.output_filename(
            "validation_psf", patch, 5
        ).startswith("validation_psf_conv-")


@pytest.mark.parametrize("bad", ["v2", "2.0", "v1.7", ""])
def test_invalid_version_raises(bad):
    """A mistyped -V is rejected rather than falling through to v1.x."""
    obj = collate_star_cat.Convert()
    obj._params["version_cat"] = bad
    with pytest.raises(ValueError):
        obj.run()
