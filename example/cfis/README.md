# CFIS PSF-validation and post-processing configs

These configuration files cover PSF validation (`Ms`, `Pl`, `MsPl`,
`valjoint`) and assorted post-processing / statistics steps that are **not
(yet) covered by the Snakemake workflow**.  They are kept here for reference
and manual use with `shapepipe_run -c <file>`.

The per-module pipeline configurations driven by the Snakemake workflow live
in [`workflow/config/cfis/`](../../workflow/config/cfis) instead; see
[`workflow/README.md`](../../workflow/README.md).

The shared SExtractor / PSFEx / catalogue-parameter files (`default.sex`
variants, `default.param`, `default.psfex`, `default.conv`, `final_cat.param`,
`star_selection.setools`) also live in `workflow/config/cfis/`, which is what
`$SP_CONFIG` points at.  None of the configs kept here reference them by path,
but if you add one that does, point it at `workflow/config/cfis/` rather than
copying the file back into this directory.
