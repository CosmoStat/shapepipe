# Module Development

New modules can be implemented in the pipeline by simply writing a *module runner*.

The basic requirement for a new module runner is a single function decorated with the `module_runner` wrapper that outputs the module `stdout` and `stderr`. *e.g.*:

```python
@module_runner()
def example_module(*args)

# DO SOMETHING

return stdout, stderr
```

The module runner decorator takes the following keyword arguments:

1. `version` : (`str`) The module version. Default value is `'0.0'`.
2. `input_module` :  (`str` or `list`) The name of a preceding module(s) whose output provide(s) the input to this module. Default value is `None`.
3. `file_pattern` : (`str` or `list`) The input file pattern(s) to look for. Default value is `''`.
4. `file_ext` : (`str` or `list`) The input file extensions(s) to look for. Default value is `''`.
5. `depends` : (`str` or `list`) The Python package(s) the module depends on. Default value is `[]`.
6. `executes` : (`str` or `list`) The system executable(s) the module implements. Default value is `[]`.
7. `numbering_scheme` : (`str`) The numbering scheme implemented by the module to find input files.
9. `run_method` : (`str`) The method by which the given module should be run. The options are `parallel` and `serial`. Default value is `parallel`.

The arguments passed to the module runner are the following:

1. `input_file_list` : The list of input files.
2. `run_dirs` : The run directories for the module output files.
3. `file_number_string` : The number pattern corresponding to the current process.
4. `config` : The config parser instance, which provides access to the configuration file parameter values. Module specific parameters can be passed using the following structure:

```python
 parameter_value = config.get(module_config_sec, 'PARAMETER')
```

5. `module_config_sec` : The name of the configuration file section for the current module.
6. `w_log` : The worker log instance, which can be used to record additional messages in the module output logs using the following structure:

```python
w_log.info('MESSAGE')
```
