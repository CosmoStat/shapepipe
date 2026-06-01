---
name: v2.0 ShapePipe run — workflow wishlist
tags:
    - shapepipe
    - planning
created-at: 2026-05-31T22:09:28.586302759+02:00
outcome: Forward-looking design wishlist for a future 'v2.0' production run, rescued from the obsolete docs/source/work_flow_v2.0.md before deletion (it was an unrealized plan sitting in published docs). Not yet implemented.
---

## The wishlist (from work_flow_v2.0.md, 2020-era)

A future "v2.0" run was sketched with these goals — none yet implemented; kept
here so the intent isn't lost when the stale docs page is removed:

- **Run per tile** — pipeline runs a single tile end-to-end, or to any
  user-specified stopping point.
- **Central exposure store** — used exposures stored once in a central place and
  linked from each tile, rather than the current mix of per-patch central storage
  and per-tile copies.
- **Drop the patch organisation** — runs no longer organised into patches.
- **Use `/scratch`** for processing speed.
- **Smarter file layout** — ~20000 tiles and a few ×10⁴ single-HDU exposure files
  need subdir bucketing, e.g. `tiles/123/123.456`, `exp/12/…`.
- **Metadata over symlinks** — replace the symlink-based "used exposures per tile"
  with a YAML (or DB) manifest.
- **Update ngmix** (can come later).
- **External detection catalogue** — replace SExtractor detection on tiles with an
  external ASCII catalogue; still run SExtractor on exposures to identify stars.
- **Remove spread-model (SM)** — note: `spread_model` still exists in the code as
  of 2026-05, so this remains unrealized.
