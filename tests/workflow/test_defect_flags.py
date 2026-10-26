"""The committed detection config supplies the adaptive-flag footprint."""

from pathlib import Path
import re

import numpy as np
import pytest

from tests.workflow.harness import Campaign
from tests.workflow.test_dag import _committed_value


@pytest.mark.parametrize("input_type", ["data", "image_sims"])
@pytest.mark.parametrize("blend", [None, "noisefill", "uberseg"])
def test_committed_detection_supplies_defect_flag_footprint(
        tmp_path, resolve_dag, monkeypatch, input_type, blend):
    """Every handling enables SEG_VIGNET and Tile_cat consumes that column.

    This guards the shipped path, rather than assigning tile.seg by hand:
    removing the workflow export, disabling the committed ini option, or
    dropping the loader's SEG_VIGNET read must not silently turn an ordinary
    clean object into an unknown-footprint (-1) object. The SExtractor runner
    column-writing behavior is covered by test_sextractor_seg_vignet.
    """
    from shapepipe.modules.ngmix_package.ngmix import Postage_stamp, Tile_cat
    from tests.module.test_ngmix_uberseg import _write_tile_cat

    campaign = Campaign(tmp_path / "campaign", input_type, "psfex")
    if blend is not None:
        campaign.config["blend_handling"] = blend
    campaign.write_config()
    with resolve_dag(campaign) as dag:
        job = next(job for job in dag.jobs if job.rule.name == "tile_detect")
        exports = dict(re.findall(r"^export (SP_SEG_VIGNET)='([^']*)'$",
                                  job.params.pre, re.MULTILINE))
        config_dir = (Path(__file__).parents[2] / "workflow" / "config"
                      / ("cfis" if input_type == "data" else "cfis_image_sims"))
        enabled = _committed_value(job.shellcmd, config_dir, "SEG_VIGNET",
                                   exports, monkeypatch)
    assert enabled == "True", "committed detection disables SEG_VIGNET"
    seg = np.zeros((3, 5, 5), dtype=np.int32)
    seg[:, 2, 2] = [4, 5, 9]
    path = _write_tile_cat(tmp_path / "tile.fits", seg)
    tile = Tile_cat(str(path))
    stamp = Postage_stamp()
    stamp.record_footprint_defects(
        np.zeros((5, 5), dtype=np.int16), np.zeros((5, 5), dtype=bool),
        None if tile.seg is None else tile.seg[0], tile.obj_id[0],
    )
    assert stamp.imaflags_iso == stamp.n_epoch_interp_iso == 0
