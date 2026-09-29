# Exposure-level maps

A workflow campaign builds two HealSparse maps from its exposures, beside the
merged catalogues on the products root. Both are per-exposure records combined
by one campaign job, at the mask ladder's resolution (`nside` 131072 over
`nside_coverage` 128, set once under `exposure_maps:` in
`workflow/config.yaml`), and both are produced automatically whenever the input
can support them — no configuration is needed.

| map | chain | product | built for |
| --- | --- | --- | --- |
| defect | `exp_split` → `exp_defect_map` → `defect_map_merge` | `defect_map/defect_map_<run>.hsp` | data runs |
| exposure count | `exp_persist` → `exp_footprint` → `nexp_map` | `nexp_map/nexp_map_<run>.hsp` | fitted-PSF runs |

**The defect map** is boolean, `True` = masked: every healpix pixel touched by
a flagged CCD pixel (bad columns, saturated pixels, bleed trails) of any
exposure. It carries the instrument flags, which otherwise never leave the pixel
domain, into the same form as the sky masks. Image simulations' flag images are
blank, so sims build none.

**The exposure-count map** counts, per pixel, the exposures with a valid PSF
model covering it — the count behind sp_validation's `npoint >= 3` cut.
`exp_footprint` records the sky corners of each CCD whose PSF fit succeeded;
`nexp_map` stamps every record on the products root into the map, so it grows as
tiles are appended. `psf_model: fake` fits no PSF, so sims build none.
`exposure_maps.nexp.enabled: false` opts a campaign out; the map is rebuilt
whole, so a campaign appended in many small batches may prefer to build it once
at the end.

Nothing in the workflow reads either map back. Plot the count map by hand with
`plot_coverage_map -i <products_dir>/nexp_map/nexp_map_<run>.hsp ...`, using the
sky windows under `exposure_maps.nexp.plot`. The design and the measurements
are in `workflow/README.md`.
