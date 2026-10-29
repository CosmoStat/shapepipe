"""A zero-row catalogue uses the resolved DAG's real partition commands.

tile_vignets' prelude runs ngmix_range before the modules, so planning alone
does not exercise the empty-tile path (Fir tile 177.257).
Execute that writer and every scheduled reader from the resolved job params,
then feed their exported bounds to the real module-to-publication regression.
"""

import json
import os
from pathlib import Path
import subprocess

import h5py

from tests.module import test_zero_object_tile as modules
from tests.workflow.harness import Campaign, resolve


def test_zero_object_tile_dag_partition_through_publication(tmp_path, monkeypatch):
    campaign = Campaign(tmp_path / "campaign", "data", "psfex")
    campaign.config["ngmix_chunks"] = 16  # Fir's configured static fan-out
    campaign.write_config()
    tile = "124.456"  # one exposure in the campaign fixture
    run = campaign.run_dir / "tiles" / tile[:2] / tile
    local = tmp_path / "node-local"
    local.mkdir()
    selected = lambda jobs: [j for j in jobs if j.wildcards.tile == tile]
    seen = []
    real_runner = modules.ngmix_runner

    with resolve(campaign, monkeypatch) as dag:
        vignets, = selected(dag.jobs_for("tile_vignets"))
        chunks = selected(dag.jobs_for("tile_ngmix"))
        merged, = selected(dag.jobs_for("tile_merge_cats"))
        made, = selected(dag.jobs_for("tile_make_cat"))
        assert len(chunks) == 16
        assert set(map(str, merged.input.manifests)) == {
            str(j.output.manifest) for j in chunks}
        assert str(merged.output.manifest) in map(str, made.input)
        assert all(str(vignets.output.manifest) in map(str, j.input)
                   for j in chunks)
        assert str(campaign.final_cat(tile)) in map(str, made.output)
        writer, = [line for line in vignets.params.pre.splitlines()
                   if "ngmix_range.py" in line and "--write" in line]
        readers = {}
        for job in chunks:
            lookup, = [line for line in job.params.pre.splitlines()
                       if "ngmix_range.py" in line and "--read" in line]
            assert 'eval "$ngmix_range_out"' in job.params.pre
            readers[int(job.wildcards.chunk)] = lookup

        def scheduled_runner(inputs, dirs, number, config, section, log):
            # The fixture has now produced its real zero-row LDAC/EPOCH tables.
            # Run the exact shell command Snakemake resolved for tile_vignets,
            # not partition() directly or a test copy of its CLI arguments.
            env = dict(os.environ, SP_RUN=str(run), SP_LOCAL=str(local))
            if not seen:
                written = subprocess.run(["bash", "-ec", writer], env=env,
                                         capture_output=True, text=True)
                assert written.returncode == 0, written.stdout + written.stderr
                doc = json.loads((local / "ngmix_ranges.json").read_text())
                assert doc["n_obj"] == 0 and doc["n_chunks"] == 16
                assert len(doc["chunks"]) == len(chunks)
            k = len(seen) + 1
            command = (readers[k] + '\neval "$ngmix_range_out"\n'
                       + 'printf "%s %s" "$NGMIX_ROW_MIN" "$NGMIX_ROW_MAX"')
            read = subprocess.run(["bash", "-ec", command], env=env,
                                  capture_output=True, text=True)
            assert read.returncode == 0, read.stdout + read.stderr
            lo, hi = map(int, read.stdout.split())
            assert (lo, hi) == (2, 1)  # bounded empty interval, not unbounded 0
            config.set(section, "ID_OBJ_MIN", str(lo))
            config.set(section, "ID_OBJ_MAX", str(hi))
            seen.append(k)
            return real_runner(inputs, dirs, number, config, section, log)

        monkeypatch.setattr(modules, "ngmix_runner", scheduled_runner)
        modules.test_zero_object_tile_modules_and_publication(
            run.parent, monkeypatch, 12, tile=tile, n_chunks=len(chunks),
            products_dir=campaign.products_dir, exposure=campaign.ready[tile][0])
        assert seen == list(range(1, 17))
        for job in [vignets, *chunks, merged, made]:
            manifest = Path(job.output.manifest)
            assert json.loads(manifest.read_text())["status"] == "complete"
        with h5py.File(campaign.final_cat(tile)) as cat:
            assert all(column.shape[0] == 0 for column in cat.values())
        stats = json.loads(Path(made.output.stats).read_text())
        assert stats["catalogue_match"]["n_sextractor"] == 0
        assert stats["epoch_cuts"]["n_chunks"] == 16
