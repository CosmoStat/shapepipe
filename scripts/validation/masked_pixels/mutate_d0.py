"""Restore a 7/13-px epoch rule in a private source copy; require its death.

Run inside the campaign SIF on a compute allocation, from the repository root::

    python scripts/validation/masked_pixels/mutate_d0.py --out /path/mutant.log

The checkout is untouched. A killed mutant must fail assertions, not collection.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


RADIUS_MUTANT = '''def footprint_defect_vetoes(defect, seg, object_number):
    interpolated = interpolated_defects(
        defect, np.zeros_like(defect), "uberseg", "interpolate"
    )
    return bool(np.any(defect & (centre_distance(defect.shape)
                                < np.where(interpolated, 7, 13))))


'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    with tempfile.TemporaryDirectory(prefix="d0-mutant-") as tmp:
        source = Path(tmp) / "src"
        shutil.copytree(root / "src", source)
        path = source / "shapepipe/modules/ngmix_package/ngmix.py"
        text = path.read_text()
        start = text.index("def footprint_defect_vetoes(")
        stop = text.index("def centre_distance(", start)
        path.write_text(text[:start] + RADIUS_MUTANT + text[stop:])
        env = dict(os.environ, PYTHONPATH=f"{source}:{root}")
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-o", "addopts=",
             "tests/unit/test_footprint_defect_veto.py",
             "tests/module/test_ngmix_defect_fill.py", "-k", "footprint"],
            cwd=root, env=env, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    Path(args.out).write_text(result.stdout)
    killed = (result.returncode == 1 and "AssertionError" in result.stdout
              and "ERROR collecting" not in result.stdout)
    print(f"radius mutant {'KILLED' if killed else 'SURVIVED/INVALID'}; log: {args.out}")
    return 0 if killed else 1


if __name__ == "__main__":
    raise SystemExit(main())
