"""Measure footprint defect flags without fitting, using retained tile stores.

Run inside the ShapePipe container on a compute allocation, with PYTHONPATH=src::

    python scripts/validation/masked_pixels/tile_defect_flags.py \
        --sexcat SEXCAT --segcat SEG_VIGNET_CAT --image IMAGE_SQLITE \
        --flag FLAG_SQLITE --weight WEIGHT_SQLITE --psf GALAXY_PSF_SQLITE \
        --headers HEADERS_SQLITE --background BACKGROUND_SQLITE \
        --rms BACKGROUND_RMS_SQLITE --blend-handling uberseg \
        --nflagged NFLAGGED.hsp --out /scratch/cdaley/defectflag-check/TILE \
        --plot scripts/validation/masked_pixels/tile_defect_flags-TILE.png

All stores open read-only; nothing is written beneath an input directory.
Omit --segcat only if SEXCAT itself carries SEG_VIGNET. Omit --background
only for already background-subtracted images; omit --rms only when the run
used stamp-estimated RMS. Pass the production blend/fill choices explicitly.
The output contains every tile row, including unavailable PSFs, zero-kept-
epoch objects and unknown footprints, distinguished by status and -1.
The Poisson expectation uses 0.0028 as the exposure DATASEC pixel flag rate
from #887. It is an unconditioned reference, not a calibration: spatially
correlated defects and the central veto violate its independence assumption.
"""

import argparse
from collections import Counter
from contextlib import ExitStack
import json
from pathlib import Path
from types import SimpleNamespace

import fitsio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package.ngmix import (
    EPOCH_INTERPOLATED_DEFECT_RADIUS, Tile_cat, centre_distance,
    prepare_postage_stamps,
)


def measure(args):
    """Read production stamps and invoke the pipeline's epoch preparation."""
    tile = Tile_cat(args.sexcat, seg_cat_path=args.segcat)
    if tile.seg is None:
        with fitsio.FITS(args.sexcat) as cat:
            objects = next(h for h in cat if h.get_extname() == "LDAC_OBJECTS")
            if "SEG_VIGNET" in objects.get_colnames():
                tile.seg = objects.read_column("SEG_VIGNET")
    if tile.seg is None:
        raise ValueError("A real-tile validation requires segmentation stamps")
    n = len(tile.obj_id)
    data = dict(number=tile.obj_id, ra=tile.ra, dec=tile.dec)
    for key in ("imaflags_iso", "n_epoch_interp_iso", "n_epoch_interp",
                "n_kept", "area", "area_outside_veto"):
        data[key] = np.full(n, -1, dtype=np.int64)
    status, dtypes = Counter(), Counter()
    with ExitStack() as stack:
        def store(path):
            if path is None:
                return None
            if not Path(path).is_file():
                raise FileNotFoundError(path)
            db = SqliteDict(path, flag="r")
            stack.callback(db.close)
            return db

        vignet = SimpleNamespace(
            gal_vign_cat=store(args.image), flag_vign_cat=store(args.flag),
            weight_vign_cat=store(args.weight), psf_vign_cat=store(args.psf),
            f_wcs_file=store(args.headers), bkg_vign_cat=store(args.background),
            bkg_rms_vign_cat=store(args.rms),
        )
        radius = centre_distance(tile.seg[0].shape)
        for i, obj_id in enumerate(tile.obj_id):
            footprint = tile.seg[i] == obj_id
            data["area"][i] = footprint.sum()
            data["area_outside_veto"][i] = (
                footprint & (radius >= EPOCH_INTERPOLATED_DEFECT_RADIUS)
            ).sum()
            key = str(obj_id)
            if key not in vignet.psf_vign_cat or key not in vignet.gal_vign_cat:
                status["missing_store_row"] += 1
                continue
            psf, gal = vignet.psf_vign_cat[key], vignet.gal_vign_cat[key]
            if isinstance(psf, str) or isinstance(gal, str):
                status["empty_store_row"] += 1
                continue
            stamp = prepare_postage_stamps(
                vignet, obj_id, i, tile, bkg_sub=args.background is not None,
                psf_obj=psf, gal_obj=gal, blend_handling=args.blend_handling,
                defect_fill=args.defect_fill,
            )
            data["n_kept"][i] = len(stamp.flags)
            for flag in stamp.flags:
                dtypes[str(flag.dtype)] += 1
            if not stamp.flags:
                status["no_kept_epoch"] += 1
                continue
            for name, value in stamp.defect_diagnostics().items():
                if name in data:
                    data[name][i] = value
            status["unknown_footprint" if stamp.imaflags_iso < 0 else "known"] += 1
            if i % 1000 == 0:
                print(f"{i}/{n}", flush=True)
    valid = data["imaflags_iso"] >= 0
    kept = data["n_kept"] > 0
    interp = data["n_epoch_interp_iso"] > 0
    flagged = data["imaflags_iso"] > 0
    wide = data["n_epoch_interp"] > 0
    summary = dict(
        inputs=vars(args), n_objects=n, status=dict(status),
        kept_flag_vignet_dtypes=dict(dtypes),
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
        for area in ("area", "area_outside_veto"):
            summary[f"poisson_{area}"] = float(np.mean(
                -np.expm1(-0.0028 * data[area][valid] * data["n_kept"][valid])
            ))
        summary["mean_area_known"] = float(data["area"][valid].mean())
        summary["mean_kept_epochs_known"] = float(data["n_kept"][valid].mean())
    if kept.any():
        summary["fraction_stamp_interp_kept"] = float(wide[kept].mean())
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
    fig, ax = plt.subplots(figsize=(9, 8), facecolor=PAPER)
    ax.set_facecolor(PAPER)
    im = ax.imshow(image, origin="lower", extent=extent, cmap="Blues",
                   interpolation="nearest", aspect="auto")
    ax.scatter(x[~selected], dec[~selected], s=0.6, c=INK2, alpha=.12,
               linewidths=0, label="Other tile detections")
    ax.scatter(x[selected], dec[selected], s=13, c=ORANGE, linewidths=.3,
               edgecolors="white", label=f"Interpolated defect on footprint ({selected.sum():,})")
    ax.set(xlabel=f"RA offset from {centre % 360:.3f}° [deg]",
           ylabel="Declination [deg]", title="Kept-epoch defects on tile-object footprints")
    ax.invert_xaxis()
    ax.legend(loc="upper right", markerscale=2)
    fig.colorbar(im, ax=ax, label="Flagged exposures per sky pixel (nflagged)")
    fig.tight_layout()
    Path(args.plot).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, dpi=180)
    plt.close(fig)
    Path(args.out + ".json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sexcat", "image", "flag", "weight", "psf", "headers", "out"):
        parser.add_argument("--" + name, required=True)
    for name in ("segcat", "background", "rms", "nflagged", "plot"):
        parser.add_argument("--" + name)
    parser.add_argument("--blend-handling", choices=("noisefill", "uberseg"), required=True)
    parser.add_argument("--defect-fill", choices=("interpolate", "noise"), default="interpolate")
    measure(parser.parse_args())


if __name__ == "__main__":
    main()
