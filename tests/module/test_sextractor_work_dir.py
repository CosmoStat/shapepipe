"""UNIT TESTS FOR SEXTRACTOR_RUNNER'S WORK_DIR.

With ``WORK_DIR`` set, SExtractor writes its catalogue and check images in a
fresh directory under it, the SEG_VIGNET cut and the join rewrite the
catalogue there, and the files are moved to the run's output directory once
complete. The published files must be the bytes the in-place run writes, the
work directory must not outlive the run, and a failed run must publish
nothing.

SExtractor is replaced by the CROWD scene of test_sextractor_seg_vignet,
written to the paths the command line names.
"""

import logging
import re

import numpy as np
import pytest
from astropy.io import fits

from shapepipe.modules import sextractor_runner as runner_module
from shapepipe.pipeline.config import CustomParser
from tests.module.test_sextractor_seg_vignet import (
    EXT_HEADER,
    _write_crowd,
    _write_crowd_external,
)

NUM = "-001-001"
FILES = [f"sexcat{NUM}.fits", f"background{NUM}.fits",
         f"segmentation{NUM}.fits"]


def _run(root, monkeypatch, work_dir=None, ext=None):
    """sextractor_runner under blend_handling: uberseg with the DR6 join, in
    place or through ``work_dir``; returns the SExtractor command line."""
    out = root / "output"
    tmp = root / "tmp"
    out.mkdir(parents=True)
    tmp.mkdir()
    calls = []

    def fake_execute(command_line):
        calls.append(command_line)
        cat = re.search(r"-CATALOG_NAME (\S+)", command_line).group(1)
        background, seg = re.search(
            r"-CHECKIMAGE_NAME (\S+)", command_line).group(1).split(",")
        _write_crowd(cat, seg)
        fits.PrimaryHDU(np.zeros((30, 30), np.float32)).writeto(background)
        return "", "All done"

    monkeypatch.setattr(runner_module, "execute", fake_execute)
    dot_param = root / "default.param"
    dot_param.write_text("NUMBER\nX_IMAGE\nY_IMAGE\nVIGNET(9,9)\n")
    if ext is None:
        ext = _write_crowd_external(root / f"CFIS_cat{NUM}.cat")
    section = {
        "EXEC_PATH": "source-extractor",
        "DOT_SEX_FILE": "d.sex", "DOT_PARAM_FILE": str(dot_param),
        "DOT_CONV_FILE": "d.conv", "WEIGHT_IMAGE": "True",
        "FLAG_IMAGE": "False", "PSF_FILE": "False",
        "DETECTION_IMAGE": "False", "DETECTION_WEIGHT": "False",
        "ZP_FROM_HEADER": "False", "BKG_FROM_HEADER": "False",
        "CHECKIMAGE": "BACKGROUND, SEGMENTATION",
        "SEG_VIGNET": "True", "MATCH_CATALOGUE": str(ext),
        "MATCH_RADIUS": "1.0", "MATCH_MIN_FRACTION": "0.98",
        "MATCH_TOLERATED_UNPAIRED": "20", "MAKE_POST_PROCESS": "False",
    }
    if work_dir is not None:
        section["WORK_DIR"] = str(work_dir)
    config = CustomParser()
    config["SEXTRACTOR_RUNNER"] = section
    runner_module.sextractor_runner(
        [str(root / f"image{NUM}.fits"), str(root / f"weight{NUM}.fits")],
        {"output": str(out), "tmp": str(tmp)}, NUM, config,
        "SEXTRACTOR_RUNNER", logging.getLogger("test"),
    )
    return calls[0]


def test_work_dir_publishes_the_in_place_bytes(tmp_path, monkeypatch):
    """The output directory receives the same files, byte for byte, and
    SExtractor wrote them under WORK_DIR, which is left empty."""
    _run(tmp_path / "in_place", monkeypatch)
    work = tmp_path / "local"
    work.mkdir()
    command_line = _run(tmp_path / "staged", monkeypatch, work_dir=work)

    assert f"-CATALOG_NAME {work}/sp-detect{NUM}." in command_line
    assert str(tmp_path / "staged" / "output") not in command_line
    assert list(work.iterdir()) == []
    staged = tmp_path / "staged" / "output"
    assert sorted(p.name for p in staged.iterdir()) == sorted(FILES)
    for name in FILES:
        assert ((staged / name).read_bytes()
                == (tmp_path / "in_place" / "output" / name).read_bytes())


def test_failed_run_publishes_nothing(tmp_path, monkeypatch):
    """A join that refuses the catalogue leaves neither a partial sexcat in
    the output directory nor the work files on the node."""
    ext = tmp_path / "dup.cat"
    ext.write_text(EXT_HEADER + "".join(
        f"{n:10d} {x:11.4f} {y:11.4f}\n"
        for n, x, y in [(5, 10.0, 10.0), (5, 10.0, 13.0), (6, 22.0, 22.0)]))
    work = tmp_path / "local"
    work.mkdir()
    with pytest.raises(ValueError, match="NUMBER"):
        _run(tmp_path / "staged", monkeypatch, work_dir=work, ext=ext)
    assert list(work.iterdir()) == []
    assert list((tmp_path / "staged" / "output").iterdir()) == []
