# Workflow DAG checks

Run these checks inside the development container:

```bash
python -m pytest tests/workflow -o addopts='' -q -p no:cacheprovider
```

`testpaths = ["tests"]` includes this directory in the full suite and the image-build CI run.
These tests check planning, not job execution, container validity, or SLURM group execution.
They need no survey files, cluster access, or nested Apptainer process.

## Campaign and API

`harness.Campaign` writes run YAML, a tile list, exposure lists, and empty prepared-stage manifests under `tmp_path`.
Two ready tiles use two exposures, one of them shared; the declared list also contains a duplicate and an unready tile.
An out-of-scope consumer has finished its vignets, and another consumer is explicitly ignored.
`build_index.build()` seeds those two consumers in the SQLite history; the Snakefile's compute parse builds the current batch's index from the exposure lists.
A prepare parse does not build the index, so no prepare invocation is necessary.

The run configuration exercises `machine: candide` and `$base_dir`/`$run` expansion for both input types.
Scratch and persistent roots differ, and the products directory's basename deliberately differs from `run:`.
The three modes are `data+psfex`, `data+mccd`, and `image_sims+fake`.
`final_cat_merge` is present in all three: it merges galaxy catalogues regardless of the PSF source.
Only `exp_persist` and `star_cat_merge` disappear for fake PSFs.

`resolve()` uses the same `SP_PHASE`, `SP_PROFILE`, `SP_RUN_CONFIG`, image selection, and state-directory conventions as `workflow/bin/sp`.
It isolates the environment, source cache, bare script imports, and Snakemake's shared global namespace.
An empty image sentinel satisfies parse-time image resolution without deploying software.
The candide profile supplies the resource defaults, resource overrides, Apptainer arguments, and rerun triggers; no executor starts.

The API sequence is `SnakemakeApi()` → `api.workflow(...)` → `workflow_api.dag(DAGSettings(targets={"all"}, ...))` → `dag_api.printdag()`.
`print_dag_as="dot"` matches Snakemake 9's renderer, which compares the CLI string rather than its enum default.
The resolved jobs come from `workflow_api._workflow.dag.jobs`; this last access is private because the public API exposes graph printing but not Job objects.
`ResolvedDAG` exposes the active rule set, declared rule set, namespace, and `jobs_for(rule)` for input/output, wildcard, params, and rendered-shell inspection.
Jobs retain the rules' thread counts rather than a local `--cores` cap.
The API context remains open while tests inspect jobs and closes before the fixture restores process state.

## Campaign-boundary pin

`params_pin.json` pins SHA-256 digests of:

- `unit_pre()` rendered for every stage;
- every rule's shell template;
- every resolved job's `params.pre` and formatted shell, including all ngmix chunks.

Rules without `params.pre` use null; aggregation-only rules have no shell.
The pin includes per-stage and per-rule digests to identify which strings change.
Fixture and checkout roots become `<CAMPAIGN_ROOT>` and `<REPO>` before hashing; two independently located campaigns must give the same digest.
The normalized strings are also written to the fixture's `params_rendered.json` for inspection.
Script fingerprints and non-pre params that do not appear in a shell are outside this pin's scope.
Separate checks require `params` in both profiles' rerun triggers.

A pin change requires a campaign boundary on a fresh root.
Review the rendered strings, then regenerate and commit the pin with the intentional workflow change:

```bash
python -m pytest tests/workflow/test_params_pin.py --update-params-pin \
  -o addopts='' -q -p no:cacheprovider
```

## Failure modes and mutation probes

Each check has a mutation that must make it fail.
Apply mutations only to a disposable checkout, run the named test without `--update-params-pin`, and require assertion failures rather than collection/setup errors.

| Check | Mutation probes |
|---|---|
| `test_rule_set_matches_input_mode` | Invert `PERSISTS_PSF`; omit `star_cat_targets()`; rename `final_cat_merge` to `merge_final_cats`. |
| `test_clean_exposure_waits_on_persist_iff_psf` | Drop the persist edge; make it unconditional under fake PSFs; drop vignets consumers; remove the in-scope consumer filter. |
| `test_final_cat_merge_reads_every_ready_tile` | Drop one ready tile; append an out-of-scope tile. |
| `test_products_use_products_dir_and_run_name` | Rename either merged catalogue or the persist manifest; route products to scratch; derive `CAMPAIGN` from the products directory's basename. |
| `test_missing_run_fails_during_parse` | Remove `run` from `run_config.REQUIRED`; literal paths must still receive the required-key diagnostic, not a later `KeyError`. |
| `test_unit_pre_changes_at_campaign_boundary` | Append a line to `unit_pre`; change one rule's `params.pre`; change one shell; change a rendered thread count. |
| `test_params_pin_ignores_fixture_root` | Remove fixture-root normalization. |
| `test_params_is_a_rerun_trigger_under_both_profiles` | Remove `params` from candide or nibi's `rerun-triggers`. |
