# Module Development

This page provides details on how new modules can be implemented in ShapePipe.

## Module Package

Every ShapePipe module requires a *module package*, which is simply a directory
with the module name followed by `_package`, e.g. for a module called
`my_new_module` you would create a new directory called `my_new_module_package`
and put it in `shapepipe/modules`. Inside this directory you will need to
include a Python file (ideally named after your module, e.g.
`my_new_module.py`) and a `__init__.py` file with the following content.

```python
"""MY NEW MODULE PACKAGE.

This package contains the module for ``my_new_module``.

:Author:

:Parent module:

:Input:

:Output:

Description
===========


Module-specific config file entries
===================================

"""

__all__ = ['my_new_module.py']
```

You should provide a description of what your module does, what the config file
options are, the inputs and outputs, and any modules this module may depend on.

The Python file (e.g. `my_new_module.py`) will contain all of the code for
implementing the module operations.

## Module Runner

In addition to the module package, new modules will also require a
*module runner*.

The basic requirement for a new module runner is a single function decorated
with the `module_runner` wrapper that outputs the module `stdout` and
`stderr`.

```python
from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.module_name_package.module_name import ...


@module_runner(*args)
def example_module(*args):

    # DO SOMETHING

    return stdout, stderr
```

In the specific case of a module that executes an executable available on the
system, the module runner should also import the `execute` function.

```{note}
:class: margin
If no `stdout` or `stderr` are provided by the given module, the module
runner should simply return `None, None`.
```

```python
from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.module_name_package.module_name import ...
from shapepipe.pipeline.execute import execute


@module_runner(*args)
def example_module(*args):

    # DO SOMETHING
    command_line = ...
    # Execute command line
    stderr, stdout = execute(command_line)

    return stdout, stderr
```



The module runner decorator takes the following keyword arguments:

1. `version` : (`str`) The module version. Default value is `'0.0'`.
2. `input_module` :  (`str` or `list`) The name of a preceding module(s)
   whose output provide(s) the input to this module.
3. `file_pattern` : (`str` or `list`) The input file pattern(s) to look for.
4. `file_ext` : (`str` or `list`) The input file extensions(s) to look for.
5. `depends` : (`str` or `list`) The Python package(s) the module depends on.
6. `executes` : (`str` or `list`) The system executable(s) the module
   implements.
7. `numbering_scheme` : (`str`) The numbering scheme implemented by the module
   to find input files.
9. `run_method` : (`str`) The method by which the given module should be run.
   The options are `parallel` and `serial`. Default value is `parallel`.

The arguments passed to the module runner are the following:

1. `input_file_list` : The list of input files.
2. `run_dirs` : The run directories for the module output files.
3. `file_number_string` : The number pattern corresponding to the current
   process.
4. `config` : The config parser instance, which provides access to the
   configuration file parameter values. Module specific parameters can be
   passed using the following structure:

```python
 parameter_value = config.get(module_config_sec, 'PARAMETER')
```

5. `module_config_sec` : The name of the configuration file section for the
   current module.
6. `w_log` : The worker log instance, which can be used to record additional
   messages in the module output logs using the following structure:

```python
w_log.info('MESSAGE')
```
