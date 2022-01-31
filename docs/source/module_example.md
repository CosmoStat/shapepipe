# Module Examples

The following example module runners are provided in `shapepipe.modules`.

## Python Example

In this example a Python script using a `Dummy()` class is implemented. This module does not read inputs from any preceding module, but looks for files in the `INPUT_DIR` that match the file patterns `'numbers'` and `'letters'` with file extension `'.txt'`. This module depends on `numpy`.

As this module does not implement any system executable, it is not necessary to return a `stderr`. Instead any output content that should be recorded in the log can be returned, otherwise the module runner should simply
return `None, None`.

```python
@module_runner(
  version='1.0',
  file_pattern=['numbers', 'letters'],
  file_ext='.txt',
  depends='numpy',
  run_method='parallel',
)
def python_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
  """Define The Python Example Runner."""
  # Set output file name
  output_file_name = (
      f'{run_dirs["output"]}/pyex_output{file_number_string}.cat'
  )

  # Retrieve log message from config file
  message = config.get(module_config_sec, 'MESSAGE')

  # Create an instance of the Python example class
  py_ex_inst = python_example.PythonExample(0)

  # Read input files
  py_ex_inst.read_files(*input_file_list)

  # Write output files
  py_ex_inst.write_file(output_file_name, message)

  # Return file content and no stderr
  return py_ex_inst.content, None
```

## Executable Example

In this example the module runner call the system executable `head`. This module read input files from the `python_example` module output that match the file pattern `'process'` with file extension `'.cat'`.

```python
@module_runner(
    version='1.0',
    input_module='python_example_runner',
    file_pattern='pyex_output',
    file_ext='.cat',
    executes='head',
    run_method='parallel',
)
def execute_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Execute Example Runner."""
    command_line = f'head {input_file_list[0]}'
    output_file_name = (
        f'{run_dirs["output"]}/head_output{file_number_string}.txt'
    )

    stdout, stderr = execute(command_line)

    text_file = open(output_file_name, 'w')
    text_file.write(stdout)

    return stdout, stderr
```
