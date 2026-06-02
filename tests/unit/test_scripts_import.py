"""Smoke-import the standalone validation scripts shipped under ``scripts/``.

``test_imports.py`` walks the installed ``shapepipe`` package, but the
standalone scripts under ``scripts/`` are not part of the package, so the
package walk never reaches them. This test extends import-smoke coverage to the
validation scripts added alongside the ngmix v2.0 work: each is expected to at
least import cleanly (its ``if __name__ == "__main__"`` guard keeps ``main``
from running). A broken top-level import — e.g. importing a helper that no
longer exists in the module — fails here instead of only at run time.
"""

import importlib.util
from pathlib import Path

import pytest
import shapepipe


REPO_ROOT = Path(shapepipe.__file__).resolve().parents[2]

# Scripts added by the ngmix v2.0 work that are meant to run against the current
# shapepipe module and carry a __main__ guard (so importing them is side-effect
# free). The jupyter export scripts/jupyter/test_centroid_shift.py is excluded:
# it has no __main__ guard and executes on import.
VALIDATION_SCRIPTS = [
    "scripts/validation/centroid/centroid_bias.py",
    "scripts/validation/centroid/centroid_bias_v2.py",
    "scripts/python/fitting.py",
]


@pytest.mark.parametrize("relpath", VALIDATION_SCRIPTS)
def test_validation_script_imports_cleanly(relpath):

    path = REPO_ROOT / relpath
    if not path.exists():
        pytest.skip(f"{relpath} not found at {path}")

    spec = importlib.util.spec_from_file_location(
        f"_smoke_{path.stem}", path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
