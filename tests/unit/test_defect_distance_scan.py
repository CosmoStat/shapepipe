"""The no-veto scan must measure the fill it labels, without epoch cuts."""

import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[2] / "scripts/validation/masked_pixels/bias_vs_distance.py"
spec = importlib.util.spec_from_file_location("distance_scan", PATH)
scan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan)


def test_napoli_grid_reaches_the_centre_and_both_old_boundaries():
    jobs = scan.napoli_jobs()
    assert len(jobs) == 696
    for hlr, psf in scan.GALAXIES:
        for psf_shear in ((0.0, 0.0), (0.05, 0.02)):
            for kind in scan.KINDS:
                for fill in (("production", "noise") if kind in scan.INTERPOLATED
                             else ("production",)):
                    distances = {j[2] for j in jobs if j[0] == kind and j[1] == fill
                                 and j[3:5] == (hlr, psf) and j[5] == psf_shear}
                    assert distances == set(range(8 if fill == "production"
                                                  and kind in scan.INTERPOLATED else 14))


@pytest.mark.parametrize("fill,expected", [("noise", "noise"), ("production", "interpolate")])
def test_scan_passes_fill_explicitly_without_patching_shared_module(monkeypatch, fill, expected):
    from shapepipe.modules.ngmix_package import ngmix as ngm
    original = ngm.interpolated_defects
    calls = []

    def response(bad, **kwargs):
        calls.append(kwargs)
        assert bad[25, 25] and bad.sum() == 51
        return {key: [0.0, 0.0] for key in ("m", "m_err", "c", "c_err")}

    monkeypatch.setitem(sys.modules, "tests.helpers.defect_response",
                        SimpleNamespace(defect_response=response))
    row = scan.one(("column", fill, 0, 0.3, 0.7, (0.05, 0.02)))
    assert "error" not in row, row
    assert calls[0]["options"] == {"defect_fill": expected}
    assert calls[0]["psf_shear"] == (0.05, 0.02)
    assert row["psf_shear"] == (0.05, 0.02)
    assert ngm.interpolated_defects is original
