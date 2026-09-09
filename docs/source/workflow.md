# Snakemake Workflow

The end-to-end CFIS production run is driven by a
[Snakemake](https://snakemake.readthedocs.io) workflow that lives in the
`workflow/` directory of the repository. It replaces the older bit-coded bash
job chain (`job_sp`) that this documentation previously described.

The workflow owns:

- the per-module ShapePipe configuration files, in `workflow/config/cfis/`,
  which is what `$SP_CONFIG` points at;
- the rule graph and job dispatch, in `workflow/Snakefile` and
  `workflow/rules/`;
- the run-level settings, in `workflow/config.yaml`.

Full usage — dependencies, profiles, how to launch a run, and how the rules map
onto the ShapePipe modules — is documented in `workflow/README.md` in the
repository.

```{note}
The PSF-validation and post-processing configurations that the workflow does
not yet cover are kept separately in `example/cfis/`; see the README there and
[Post-processing](post_processing.md).
```
