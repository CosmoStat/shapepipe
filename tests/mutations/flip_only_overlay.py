"""Verify that registration regressions kill a flip-only ngmix mutant.

Run inside the test environment from any directory::

    python tests/mutations/flip_only_overlay.py

The mutant lives in a temporary source tree; the checkout is not modified.
This command succeeds only when all three selected regressions fail.
"""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
TARGET = Path("shapepipe/modules/ngmix_package/ngmix.py")
NATIVE = (
    "        tile_seg, neighbour, off_tile = "
    "(np.asarray(epoch[name]) for name in required)\n"
)
FLIP_ONLY = """        tile_seg = Ngmix.MegaCamFlip(np.copy(tile_cat.seg[i_tile]), int(ccd_n))
        tile_vign = Ngmix.MegaCamFlip(np.copy(tile_cat.vign[i_tile]), int(ccd_n))
        neighbour, off_tile = split_tile_markers(tile_vign, gal_vign.shape, tile_seg, obj_id)
"""


def main():
    """Run the independent-grid, off-tile-cut and retained-footprint guards."""
    with tempfile.TemporaryDirectory(prefix="flip-only-mutant-") as temporary:
        source = Path(temporary) / "src"
        shutil.copytree(ROOT / "src", source)
        target = source / TARGET
        text = target.read_text()
        assert text.count(NATIVE) == 1, "Native-map consumption site changed"
        target.write_text(text.replace(NATIVE, FLIP_ONLY))
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            (str(source), str(ROOT), env.get("PYTHONPATH", ""))
        )
        env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        result = subprocess.run(
            [sys.executable, "-m", "pytest",
             str(ROOT / "tests/module/test_categorical_registration.py"),
             "-k", "ngmix_consumes or off_tile_reaches or retained_188307",
             "-o", "addopts=", "-q", "--tb=short"],
            cwd=ROOT, env=env, capture_output=True, text=True,
        )
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        assert result.returncode == 1 and "3 failed" in result.stdout, (
            "The flip-only overlay mutant must fail all three regressions"
        )
        print("Flip-only mutant killed: all three regressions failed.")


if __name__ == "__main__":
    main()
