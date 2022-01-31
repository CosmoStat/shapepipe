# Configuration

The pipeline requires a configuration file (by default called `conifg.ini`)
in order to be run. An example configuration file is provided in the
`example` directory.

The configuration parameters for the pipeline are:

## Default Options

1. `VERBOSE` : (`bool`) Set the verbosity level. Default value is `True`.
2. `RUN_NAME` : (`str`) The pipeline run name. Default value is `shapepipe_run`.
3. `RUN_DATETIME` : (`bool`) Option to add date and time to `RUN_NAME`. Default value is `True`.

## Execution Options

1. `MODULE` : (`str` or `list`) A valid module runner name (or a comma separated list of names).
2. `MODE` : (`str`) The pipeline execution mode. Options are `smp` or `mpi`. Default value is `smp`.

## File Options

1. `LOG_NAME` : (`str`) Current run log file name. Default value is `shapepipe`.
2. `RUN_LOG_NAME` : (`str`) Run history log file name. Default value is `shapepipe_runs`.
3. `INPUT_DIR` : (`str` or `list`) A valid directory containing input files for the first module or a comma separated list of directories. This parameter also recognizes the following special strings:
   - `last:MODULE` : This will point to the output directory of the last run of the specified module.
   - `PATTERN:MODULE` : This will point to the output directory of a specified module from a run matching the specified pattern.
4. `OUTPUT_DIR` : (`str`) A valid directory to write the pipeline output files.
5. `FILE_PATTERN` : (`str` or `list`) A list of string patterns to identify input files for the first module.
6. `FILE_EXT` : (`str` or `list`) A list of file extensions to identify input files for the first module.
7. `NUMBERING_SCHEME` : (`str`) A string indicating the expected numbering system for the input files (*e.g.* `000-0`). Single digits indicate integer values without limit, multiples of digits indicate integers with a maximum value. Standard characters can be placed around digits (*e.g.* `.`, `-`, `:`, *etc.*). Optionally a regular expression can also be passed if it is preceded by `RE:` (*e.g.* `RE:-\d{9}`).

## Job Options

1. `SMP_BATCH_SIZE` : (`int`) Number of SMP jobs to run in parallel. Default value is `1`.
2. `TIMEOUT` : (`int`) Timeout limit in seconds for a given job.

## Module Options

Additional module options can be added using the following structure:

```ini
[MODULE_NAME]
PARAMETER = PARAMETER VALUE
```

This mechanism can also be used to modify module decorator properties or append additional values to list properties as follows:

```ini
[MODULE_NAME]
ADD_PARAMETER = PARAMETER VALUE
```
