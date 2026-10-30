"""Resolve fixed pixel source lists against durable success manifests.

An operator can explicitly abandon a source in run_dir/abandoned.txt. Missing
manifests otherwise block publication, including failed or killed jobs.
"""

import csv
from pathlib import Path


def check_exposure_pixels(source_list, exposure, pixels):
    """Reject an exposure's rows or map coverage outside its planned pixels."""
    with Path(source_list).open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not {"pixel", "exp_id"} <= set(reader.fieldnames):
            raise ValueError(f"Invalid exposure list: {source_list}")
        allowed = {int(row['pixel']) for row in reader if row['exp_id'] == exposure}
    missing = set(map(int, pixels)) - allowed
    if missing:
        raise ValueError(f"Exposure {exposure} has data in unlisted pixel {min(missing)}")


# @sc [decision:survey_publication.pixel_completeness]
def finished_sources(pixel, kind, source_list, products_dir, run_dir):
    """Return succeeded units, refusing publication while any unit is pending.

    Tiles use the durable tile_make_cat verdict; star catalogues use
    exp_persist and maps use exp_maps. Abandonment takes precedence over a
    missing manifest, never over a success manifest.
    """
    columns = {"tiles": "tile_id", "stars": "exp_id", "maps": "exp_id"}
    if kind not in columns:
        raise ValueError(f"Unknown source kind: {kind}")
    with Path(source_list).open(newline="") as stream:
        reader = csv.DictReader(stream)
        column = columns[kind]
        if reader.fieldnames is None or not {"pixel", column} <= set(reader.fieldnames):
            raise ValueError(f"Invalid source list: {source_list}")
        units = [row[column] for row in reader if int(row["pixel"]) == pixel]
    if len(units) != len(set(units)):
        raise ValueError(f"Duplicate source in pixel {pixel}")
    abandoned_path = Path(run_dir) / "abandoned.txt"
    abandoned = set()
    if abandoned_path.exists():
        abandoned = {line.split("#", 1)[0].strip()
                     for line in abandoned_path.read_text().splitlines()}
    products = Path(products_dir)
    succeeded, blocked = [], []
    for unit in units:
        if kind == "tiles":
            manifest = products / "tiles" / unit[:2] / unit / "tile_make_cat.json"
        else:
            stage = "exp_persist" if kind == "stars" else "exp_maps"
            manifest = products / "exp" / unit[:2] / unit / "manifests" / f"{stage}.json"
        if manifest.is_file():
            succeeded.append(unit)
        elif unit not in abandoned:
            blocked.append(unit)
    if blocked:
        raise RuntimeError(f"Pixel {pixel} blocked by {kind}: {', '.join(blocked)}")
    return succeeded
