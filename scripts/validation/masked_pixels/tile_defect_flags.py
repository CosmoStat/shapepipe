"""Measure footprint defect flags without fitting, using retained tile stores.

Run inside the ShapePipe container on a compute allocation, with PYTHONPATH=src::

    python scripts/validation/masked_pixels/tile_defect_flags.py \
        --sexcat SEXCAT --image IMAGE_SQLITE \
        --flag FLAG_SQLITE --weight WEIGHT_SQLITE --psf GALAXY_PSF_SQLITE \
        --headers HEADERS_SQLITE --background BACKGROUND_SQLITE \
        --rms BACKGROUND_RMS_SQLITE --blend-handling uberseg \
        --nflagged NFLAGGED.hsp --out /scratch/cdaley/defectflag-check/TILE \
        --plot scripts/validation/masked_pixels/tile_defect_flags-TILE.png

All stores open read-only; nothing is written beneath an input directory.
SEXCAT must carry the committed workflow's SEG_VIGNET column. Omit
--background only for already background-subtracted images; omit --rms when the run
used stamp-estimated RMS. Pass the production blend/fill choices explicitly.
The output contains every tile row, including unavailable PSFs, zero-kept-
epoch objects and unknown footprints, distinguished by status and -1.
The Poisson expectation uses 0.0028 as the exposure DATASEC pixel flag rate
from #887. It is an unconditioned reference, not a calibration: spatially
correlated defects and the footprint cut violate its independence assumption.
After measuring, --plot-only --out PREFIX --nflagged MAP --plot PNG regenerates
just the figure and spatial comparison from PREFIX.npz and PREFIX.json.
"""

import argparse
from collections import Counter
from contextlib import ExitStack
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package.ngmix import (
    Ngmix, Tile_cat,
    prepare_postage_stamps,
)


def measure(args):
    """Read production stamps and invoke the pipeline's epoch preparation."""
    tile = Tile_cat(args.sexcat)
    if tile.seg is None:
        raise ValueError("A real-tile validation requires segmentation stamps")
    n = len(tile.obj_id)
    data = dict(number=tile.obj_id, ra=tile.ra, dec=tile.dec)
    for key in ("imaflags_iso", "n_epoch_interp_iso", "n_epoch_interp",
                "n_kept", "area", "n_candidate",
                "precut_imaflags_iso"):
        data[key] = np.full(n, -1, dtype=np.int64)
    data["footprint_clipped"] = np.zeros(n, dtype=bool)
    status, dtypes, cuts = Counter(), Counter(), Counter()
    stamp_pixels = stamp_flagged = footprint_pixels = footprint_flagged = 0
    with ExitStack() as stack:
        def store(path):
            if path is None:
                return None
            if not Path(path).is_file():
                raise FileNotFoundError(path)
            db = SqliteDict(path, flag="r")
            if getattr(args, "cache_stores", False):
                print(f"Caching {Path(path).name}", flush=True)
                cached = dict(db.items())
                db.close()
                print(f"Cached {len(cached)} rows", flush=True)
                return cached
            stack.callback(db.close)
            return db

        vignet = SimpleNamespace(
            gal_vign_cat=store(args.image), flag_vign_cat=store(args.flag),
            weight_vign_cat=store(args.weight), psf_vign_cat=store(args.psf),
            f_wcs_file=store(args.headers), bkg_vign_cat=store(args.background),
            bkg_rms_vign_cat=store(args.rms),
        )
        # Decode exposure-level WCS objects once rather than for every object.
        vignet.f_wcs_file = dict(vignet.f_wcs_file.items())
        for i, obj_id in enumerate(tile.obj_id):
            if i % 1000 == 0:
                print(f"{i}/{n}", flush=True)
            footprint = tile.seg[i] == obj_id
            data["area"][i] = footprint.sum()
            data["footprint_clipped"][i] = (
                footprint[[0, -1], :].any() or footprint[:, [0, -1]].any()
            )
            key = str(obj_id)
            if key not in vignet.psf_vign_cat or key not in vignet.gal_vign_cat:
                status["missing_store_row"] += 1
                continue
            psf, gal = vignet.psf_vign_cat[key], vignet.gal_vign_cat[key]
            if isinstance(psf, str) or isinstance(gal, str):
                status["empty_store_row"] += 1
                continue
            flags = vignet.flag_vign_cat[key]
            data["n_candidate"][i] = len(psf)
            if footprint.any() and not data["footprint_clipped"][i]:
                bits = 0
                for name in psf:
                    own = Ngmix.MegaCamFlip(footprint, int(name.split("-")[1]))
                    values = np.asarray(flags[name]["VIGNET"])[own]
                    bits |= int(np.bitwise_or.reduce(values.astype(np.int64)))
                data["precut_imaflags_iso"][i] = bits
            # Reuse the fetched flag row when the production preparation reads
            # it; no second SQLite fetch is needed for the pre-cut comparison.
            per_object = SimpleNamespace(**(vars(vignet) | {
                "flag_vign_cat": {key: flags},
            }))
            stamp = prepare_postage_stamps(
                per_object, obj_id, i, tile, bkg_sub=args.background is not None,
                psf_obj=psf, gal_obj=gal, blend_handling=args.blend_handling,
                defect_fill=args.defect_fill,
            )
            data["n_kept"][i] = len(stamp.flags)
            cuts.update(stamp.epoch_cuts)
            for flag, epoch_seg in zip(stamp.flags, stamp.segs):
                dtypes[str(flag.dtype)] += 1
                flagged_pixels = flag != 0
                stamp_pixels += flag.size
                stamp_flagged += int(flagged_pixels.sum())
                if stamp.imaflags_iso >= 0:
                    own = epoch_seg == obj_id
                    footprint_pixels += int(own.sum())
                    footprint_flagged += int((flagged_pixels & own).sum())
            if not stamp.flags:
                status["no_kept_epoch"] += 1
                continue
            for name, value in stamp.defect_diagnostics().items():
                if name in data:
                    data[name][i] = value
            status["unknown_footprint" if stamp.imaflags_iso < 0 else "known"] += 1
    valid = data["imaflags_iso"] >= 0
    kept = data["n_kept"] > 0
    interp = data["n_epoch_interp_iso"] > 0
    flagged = data["imaflags_iso"] > 0
    wide = data["n_epoch_interp"] > 0
    summary = dict(
        inputs=vars(args), n_objects=n, status=dict(status),
        kept_flag_vignet_dtypes=dict(dtypes), epoch_cuts=dict(cuts),
        kept_stamp_flagged_pixel_fraction=(
            stamp_flagged / stamp_pixels if stamp_pixels else None
        ),
        kept_known_footprint_flagged_pixel_fraction=(
            footprint_flagged / footprint_pixels if footprint_pixels else None
        ),
        footprints=dict(
            n_edge_touching=int(data["footprint_clipped"].sum()),
            fraction_edge_touching=float(data["footprint_clipped"].mean()),
            n_empty=int((data["area"] == 0).sum()),
            n_edge_touching_with_kept_epoch=int(
                (data["footprint_clipped"] & kept).sum()
            ),
            n_unknown_with_kept_epoch=int((~valid & kept).sum()),
        ),
        n_interp_iso=int(interp.sum()), n_imaflags_iso=int(flagged.sum()),
        n_stamp_interp=int(wide.sum()),
        fractions_all_rows=dict(interp_iso=float(interp.mean()),
                                imaflags_iso=float(flagged.mean()),
                                stamp_interp=float(wide.mean())),
        bits={str(1 << bit): int(((data["imaflags_iso"][valid] & (1 << bit)) != 0).sum())
              for bit in range(63)
              if ((data["imaflags_iso"][valid] & (1 << bit)) != 0).any()},
    )
    if valid.any():
        summary["fractions_known_footprint"] = dict(
            interp_iso=float(interp[valid].mean()),
            imaflags_iso=float(flagged[valid].mean()),
        )
        summary["fraction_precut_instrument_flags_known_kept"] = float(
            (data["precut_imaflags_iso"][valid] > 0).mean()
        )
        summary["poisson_before_cuts_known_kept"] = float(np.mean(
            -np.expm1(-0.0028 * data["area"][valid] * data["n_candidate"][valid])
        ))
        for area in ("area",):
            summary[f"poisson_{area}"] = float(np.mean(
                -np.expm1(-0.0028 * data[area][valid] * data["n_kept"][valid])
            ))
        summary["mean_area_known"] = float(data["area"][valid].mean())
        summary["median_area_known"] = float(np.median(data["area"][valid]))
        summary["mean_kept_epochs_known"] = float(data["n_kept"][valid].mean())
    if kept.any():
        summary["fraction_stamp_interp_kept"] = float(wide[kept].mean())
        summary["fraction_unknown_footprint_kept"] = float((~valid)[kept].mean())
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(str(out) + ".npz", **data)
    Path(str(out) + ".json").write_text(json.dumps(summary, indent=2) + "\n")
    if args.plot:
        plot(args, data, summary)
    print(json.dumps(summary, indent=2))


def plot(args, data, summary):
    """Overlay interpolated-footprint objects on the retained nflagged map."""
    import healsparse
    import hpgeom
    from matplotlib.colors import LogNorm
    from figstyle import INK2, ORANGE, PAPER

    if not args.nflagged:
        raise ValueError("--plot requires the exposure --nflagged map")
    sky = healsparse.HealSparseMap.read(args.nflagged)
    ra, dec = data["ra"], data["dec"]
    # Unwrap a tile crossing RA=0 before laying out its sky grid.
    centre = float(np.rad2deg(np.angle(np.mean(np.exp(1j * np.deg2rad(ra))))))
    x = (ra - centre + 180) % 360 - 180
    extent = [x.min(), x.max(), dec.min(), dec.max()]
    xx, yy = np.meshgrid(np.linspace(*extent[:2], 700),
                         np.linspace(*extent[2:], 700))
    image = sky.get_values_pos((xx + centre) % 360, yy, lonlat=True).astype(float)
    image[image < 0] = np.nan
    at_objects = sky.get_values_pos(ra, dec, lonlat=True)
    selected = data["n_epoch_interp_iso"] > 0
    summary["map_nflagged_positive_at_interp_iso"] = int(
        ((at_objects > 0) & selected).sum()
    )
    summary["map_grid_fraction_positive"] = float((image > 0).mean())
    summary["map_fraction_positive_at_other_centroids"] = float(
        (at_objects[~selected] > 0).mean()
    )
    # A flagged footprint pixel need not share the centroid's 1.6-arcsec
    # HealSparse cell. Eight arcsec covers a 51-px stamp's half-diagonal
    # (6.6 arcsec) plus one map-cell half-diagonal.
    def nearby(rows):
        return sum(bool((sky.get_values_pix(hpgeom.query_circle(
            sky.nside_sparse, ra[i], dec[i], 8 / 3600, inclusive=True,
        )) > 0).any()) for i in rows)

    selected_rows = np.flatnonzero(selected)
    controls = np.flatnonzero(~selected & (data["imaflags_iso"] >= 0))
    rng = np.random.default_rng(0)
    controls = rng.choice(controls, size=min(len(controls), 2000), replace=False)
    summary["map_positive_within_8_arcsec"] = dict(
        flagged_n=len(selected_rows), flagged_matches=nearby(selected_rows),
        control_n=len(controls), control_matches=nearby(controls),
    )
    fig, ax = plt.subplots(figsize=(9, 8), facecolor=PAPER)
    ax.set_facecolor(PAPER)
    im = ax.imshow(image, origin="lower", extent=extent, cmap="Blues",
                   norm=LogNorm(vmin=1, vmax=max(2, np.nanmax(image))),
                   interpolation="nearest", aspect="auto")
    ax.scatter(x[~selected], dec[~selected], s=0.6, c=INK2, alpha=.12,
               linewidths=0, label="Other tile detections")
    ax.scatter(x[selected], dec[selected], s=13, c=ORANGE, linewidths=.3,
               edgecolors="white", label=f"Interpolated defect on footprint ({selected.sum():,})")
    ax.set(xlabel=f"RA offset from {centre % 360:.3f}° [deg]",
           ylabel="Declination [deg]",
           title=f"{Path(args.out).name}: kept-epoch defects on object footprints")
    ax.invert_xaxis()
    ax.legend(loc="upper right", markerscale=2)
    fig.colorbar(im, ax=ax, label="Flagged CCD pixels per sky pixel (nflagged)")
    fig.tight_layout()
    Path(args.plot).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, dpi=180)
    plt.close(fig)
    Path(args.out + ".json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sexcat", "image", "flag", "weight", "psf", "headers", "out"):
        parser.add_argument("--" + name, required=(name == "out"))
    for name in ("background", "rms", "nflagged", "plot"):
        parser.add_argument("--" + name)
    parser.add_argument("--blend-handling", choices=("noisefill", "uberseg"))
    parser.add_argument("--defect-fill", choices=("interpolate", "noise"), default="interpolate")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--cache-stores", action="store_true",
                        help="Cache read-only stores in RAM before measuring")
    args = parser.parse_args()
    if args.plot_only:
        if not args.plot or not args.nflagged:
            parser.error("--plot-only requires --plot and --nflagged")
        summary = json.loads(Path(args.out + ".json").read_text())
        with np.load(args.out + ".npz") as data:
            plot(args, data, summary)
    else:
        for name in ("sexcat", "image", "flag", "weight", "psf", "headers",
                     "blend_handling"):
            if getattr(args, name) is None:
                parser.error(f"--{name.replace('_', '-')} is required for measurement")
        measure(args)


if __name__ == "__main__":
    main()
