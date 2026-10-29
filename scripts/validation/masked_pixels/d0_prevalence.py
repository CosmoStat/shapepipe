"""Measure D0 through prepare_postage_stamps on retained survey inputs.

Run in an allocation and the campaign SIF with PYTHONPATH=src:.
Epoch stores use immutable read-only connections; WCS headers open
read-only. No fits are performed.
The inputs JSON is the survey's provenance record with an ``inputs`` mapping.
"""

import argparse
from collections import Counter
import json
from multiprocessing import get_context
from pathlib import Path
import pickle
import sqlite3
import subprocess
from types import SimpleNamespace

import numpy as np
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package import ngmix as ngm

TILE = PATHS = VIGNET = None


class ReadOnlyStore:
    """Decode one sqlitedict value without opening a writer connection."""

    def __init__(self, path):
        self.db = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)

    def __getitem__(self, key):
        row = self.db.execute(
            "SELECT value FROM unnamed WHERE key=?", (str(key),)
        ).fetchone()
        if row is None:
            raise KeyError(key)
        return pickle.loads(bytes(row[0]))


def initialize():
    global VIGNET
    stores = {key: ReadOnlyStore(PATHS[key]) for key in
              ("image", "weight", "flag", "rms", "psf", "background")}
    with SqliteDict(PATHS["headers"], flag="r") as db:
        headers = dict(db.items())
    VIGNET = SimpleNamespace(
        gal_vign_cat=stores["image"], weight_vign_cat=stores["weight"],
        flag_vign_cat=stores["flag"], bkg_rms_vign_cat=stores["rms"],
        psf_vign_cat=stores["psf"], bkg_vign_cat=stores["background"],
        f_wcs_file=headers,
    )


def measure_object(i):
    number = int(TILE.obj_id[i])
    psf = VIGNET.psf_vign_cat[str(number)]
    gal = VIGNET.gal_vign_cat[str(number)]
    if not isinstance(psf, dict) or not isinstance(gal, dict):
        return i, {}, {}, []
    # Observe the actual predicate calls, after the production fraction cut.
    # Retain the real predicate and run it on every epoch, without a fallback.
    predicate = ngm.footprint_defect_vetoes
    edge = Counter()
    decisions = []

    def record(defect, seg, obj_id):
        foot = seg == obj_id
        empty = not foot.any()
        border = bool(foot[[0, -1], :].any() or foot[:, [0, -1]].any())
        reject = predicate(defect, seg, obj_id)
        edge["eligible"] += 1
        edge["empty"] += empty
        edge["border"] += border
        edge["border_rejected"] += border and reject
        decisions.append(reject)
        return reject

    ngm.footprint_defect_vetoes = record
    try:
        stamp = ngm.prepare_postage_stamps(
            VIGNET, number, i, TILE, psf_obj=psf, gal_obj=gal,
            blend_handling="uberseg", defect_fill="interpolate",
        )
    finally:
        ngm.footprint_defect_vetoes = predicate
    return i, dict(stamp.epoch_cuts), dict(edge), decisions


def main():
    global TILE, PATHS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--survey", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--nproc", type=int, default=4)
    args = parser.parse_args()
    PATHS = json.loads(Path(args.inputs).read_text())["inputs"]
    TILE = ngm.Tile_cat(PATHS["sexcat"])
    with np.load(args.survey) as data:
        reference = data["epochs"]
        assert np.array_equal(TILE.obj_id, data["objects"]["NUMBER"])
    own = np.stack([TILE.seg[i] == number for i, number in enumerate(TILE.obj_id)])
    empty = ~own.any(axis=(1, 2))
    border = own[:, [0, -1], :].any(axis=(1, 2)) | own[:, :, [0, -1]].any(axis=(1, 2))
    cuts, edge = Counter(), Counter()
    per_object = np.zeros(len(TILE.obj_id), dtype=np.int64)
    decisions = []
    with get_context("fork").Pool(args.nproc, initializer=initialize) as pool:
        for i, counts, frequencies, rejected in pool.imap(measure_object, range(len(own)), chunksize=64):
            cuts.update(counts)
            edge.update(frequencies)
            per_object[i] = counts.get("footprint_defect", 0)
            decisions.extend(rejected)
            if i % 2000 == 0:
                print(f"objects {i}/{len(own)} cuts {dict(cuts)}", flush=True)
    eligible = reference[reference["eligible"]]
    actual = np.asarray(decisions, dtype=bool)
    assert len(actual) == len(eligible) == edge["eligible"]
    assert cuts["considered"] == len(reference)
    assert cuts["masked_fraction"] == int((~reference["eligible"]).sum())
    raw_d0 = eligible["hit0"] | empty[eligible["oi"]]
    fixed = (eligible["d_interp"] < 7) | (eligible["d_noise"] < 13)
    survey_d0 = np.where(eligible["known"], eligible["hit0"], fixed)
    report = dict(
        command=" ".join(__import__("sys").argv), inputs=PATHS,
        git_sha=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        objects=len(own), objects_empty=int(empty.sum()),
        objects_border=int(border.sum()), epoch_cuts=dict(cuts),
        eligible_edge_cases=dict(edge),
        considered_empty=int(empty[reference["oi"]].sum()),
        considered_border=int(border[reference["oi"]].sum()),
        eligible_epochs=len(eligible), kept_epochs=int((~actual).sum()),
        epoch_loss_fraction=float(actual.mean()),
        survey_d0_rejected=int(survey_d0.sum()),
        survey_d0_loss_fraction=float(survey_d0.mean()),
        survey_policy_differences=int((actual != survey_d0).sum()),
        raw_overlap_mismatches=int((actual != raw_d0).sum()),
        objects_emptied=int(((np.bincount(eligible["oi"], minlength=len(own))
                             - per_object) == 0).sum()),
    )
    Path(args.out).write_text(json.dumps(report, indent=2) + "\n")
    np.savez_compressed(Path(args.out).with_suffix(".npz"), actual=actual,
                        survey_d0=survey_d0, raw_d0=raw_d0,
                        empty=empty, border=border, eligible=eligible)
    print(json.dumps(report, indent=2), flush=True)
    assert report["raw_overlap_mismatches"] == 0


if __name__ == "__main__":
    main()
