# Configuration

The pipeline requires a configuration file (by default called `conifg.ini`)
in order to be run. Example configuration files are provided in the
[example](https://github.com/CosmoStat/shapepipe/tree/develop/example)
directory.

The basic configuration parameters for the pipeline are described below.
Module-specific configuration options are provided in the API documentation.
For each option the expected data type is specified. Options labelled as
*optional* will default to the values mentioned if not specified manually and
options labelled as *required* must be provided in order for ShapePipe to
run.

## Default Options

The following options can be added to the `[DEFAULT]` section of the config
file

- `VERBOSE` : (`bool`, *optional*) Set the verbosity level. When set to `False`
  `shapepipe_run` will not print any output to the terminal. Default value is
  `True`.
- `RUN_NAME` : (`str`, *optional*) The pipeline run name. Default value is
  `shapepipe_run`.
- `RUN_DATETIME` : (`bool`, *optional*) Option to add the current date and time
  to `RUN_NAME`. Default value is `True`.

If you do not specify any options for this section ShapePipe will default to
the following settings.

```ini
[DEFAULT]
VERBOSE = True
RUN_NAME = shapepipe_run
RUN_DATETIME = True
```

## Execution Options

The following options can be added to the `[EXECUTION]` section of the config
file

- `MODULE` : (`str` or `list`, *required*) A valid module runner name (or a
  comma separated list of names).
- `MODE` : (`str`, *optional*) The pipeline execution mode. Options are `smp`
  (run with Joblib) or `mpi` (run with MPI). Default value is `smp`.

For example, if you provide the following options

```ini
[EXECUTION]
MODULE = sextractor_runner
```

the SExtractor module will be run using Joblib.

## File Options

The following options can be added to the `[FILE]` section of the config file

- `LOG_NAME` : (`str`, *optional*) Current run log file name. Default value is
  `shapepipe`.
- `RUN_LOG_NAME` : (`str`, *optional*) Run history log file name. Default value
  is `shapepipe_runs`.
- `INPUT_DIR` : (`str` or `list`, *required*) A valid directory containing
  input files for the first module or a comma separated list of directories.
  This parameter also recognises the following special strings:
   - `last:MODULE` : This will point to the output directory of the last run of
   the specified module.
   - `all:MODULE` : This will point to all the output directories in which the
   specified module was run.
   - `PATTERN:MODULE` : This will point to the output directory of a specified
   module from a run matching the specified pattern.
- `OUTPUT_DIR` : (`str`, *required*) A valid directory to write the pipeline
  output files.
- `FILE_PATTERN` : (`str` or `list`, *optional*) A list of string patterns to
  identify input files for the first module.
- `FILE_EXT` : (`str` or `list`, *optional*) A list of file extensions to
  identify input files for the first module.
- `NUMBERING_SCHEME` : (`str`, *optional*) A string indicating the expected
  numbering system for the input files (*e.g.* `000-0`). Single digits indicate
  integer values without limit, multiples of digits indicate integers with a
  maximum value. Standard characters can be placed around digits
  (*e.g.* `.`, `-`, `:`, *etc.*). *optional*ly a regular expression can also be
  passed if it is preceded by `RE:` (*e.g.* `RE:-\d{9}`).
- `NUMBER_LIST` : (`str` or `list`, *optional*) A list of number strings
  matching the numbering scheme or a file name.
- `CORRECT_FILE_PATTERN` : (`bool`, *optional*) Option to allow substring file
  patterns. Default value is `True`.

For example, if you provide the following options

```ini
[JOB]
INPUT_DIR = /home/username/my_input_dir
OUTPUT_DIR = /home/username/my_output_dir
FILE_PATTERN = galaxy
FILE_EXT = fits
NUMBERING_SCHEME = -00-0
```

ShapePipe will look for files of the form `galaxy-00-0.fits`,
`galaxy-00-1.fits`, etc. in the directory `/home/username/my_input_dir` and
save outputs to `/home/username/my_output_dir`. Note that the `FILE_PATTERN`
does not need to be complete, in other words files of the form
`mygalaxy-00-0.fits` would equally be found if `CORRECT_FILE_PATTERN = True`
were added to the options above.

Conversely, with the options

```ini
[JOB]
INPUT_DIR = last:psfex_runner
OUTPUT_DIR = /home/username/my_output_dir
FILE_PATTERN = psf
FILE_EXT = fits
NUMBERING_LIST = -001, -002
```

ShapePipe will look for the specific files `psf-001.fits` and `psf-002.fits`
in the output directory of the last run of the `psfex_runner`.

## Job Options

The following options can be added to the `[JOB]` section of the config file

- `SMP_BATCH_SIZE` : (`int`, *optional*) Number of SMP jobs to run in parallel.
  Note that this option is only valid for running in `smp` mode. Default value
  is `1`.
- `TIMEOUT` : (`str`, *optional*) Timeout limit in `HH:MM:SS` for a given job.
  If not specified no timeout limit is applied to the jobs.

For example, if you provide the following options

```ini
[JOB]
SMP_BATCH_SIZE = 4
TIMEOUT = 00:01:35
```

ShapePipe will distribute jobs across four cores and each job will time out
if not completed within 1min and 35s.

## Worker Options

The following options can be added to the `[WORKER]` section of the config file

- `PROCESS_PRINT_LIMIT` : (`int`, *optional*) The maximum number of characters
  allowed to print individual processes. Default value is `200`.

## Module Options

### Default Module Options

All ShapePipe modules accept the following options

- `INPUT_DIR` : (`str` or `list`, *optional*) Override the default input
  directories defined in `[FILE]`.
- `INPUT_MODULE` : (`str` or `list`, *optional*) Override the default input
  modules defined in the module runner.
- `FILE_PATTERN` : (`str` or `list`, *optional*) Override the default file
  pattern defined in the module runner.
- `FILE_EXT` : (`str` or `list`, *optional*) Override the default file
  extension defined in the module runner.
- `NUMBERING_SCHEME` : (`str`, *optional*) Override the default numbering
  scheme defined in the module runner.
- `DEPENDS` : (`str` or `list`, *optional*) Override the default Python
  dependencies defined in the module runner.
- `EXECUTES` : (`str` or `list`, *optional*) Override the default system
  executables defined in the module runner.

### Module-Specific Options

Additional module-specific options can be added using the following structure

```ini
[MODULE_NAME_RUNNER]
PARAMETER = PARAMETER VALUE
```

This mechanism can also be used to modify module decorator properties or append
additional values to list properties as follows

```ini
[MODULE_NAME_RUNNER]
ADD_PARAMETER = PARAMETER VALUE
```

### Multiple Module Runs

If a given module is run more than once, run specific parameter values can be
specified as follows

```ini
[MODULE_NAME_RUNNER_RUN_X]
PARAMETER = PARAMETER VALUE
```

where ``X`` is an integer greater than or equal to ``1``. This feature can be combined with the ``INPUT_DIR`` options. For example, for module *B* to access the first of two runs of module *A* you could something like set up below.

```ini

[EXECUTION]
MODULE = module_a_runner, module_b_runner, module_b_runner

[MODULE_A_RUNNER_RUN_1]
...

[MODULE_A_RUNNER_RUN_2]
...

[MODULE_B_RUNNER]
INPUT_DIR = last:module_a_runner_run_1

```
