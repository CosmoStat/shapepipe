# Testing

After reading through the [basic execution](basic_execution.md) and
[configuration](configuration.md) instructions it may be helpful to run a few
tests.

ShapePipe includes some dummy modules for testing purposes

- `python_example_runner` : an example of a module written entirely in Python
- `serial_example_runner` : an example of a module that is run in serial mode
- `execute_example_runner` : an example of a module that calls a system
  executable

None of these modules do anything particularly interesting, but they can be run
to test that ShapePipe is up and running.

The [example](https://github.com/CosmoStat/shapepipe/tree/develop/example)
directory contains an example config file called `config.ini` that runs all of
these dummy modules. ShapePipe can be run with this config file as follows

```bash
shapepipe_run -c ./example/config.ini
```

The output of this run will be saved to `./example/output`. If you want a
more concrete understanding of what each of the configuration options does you
can modify `config.ini` and see what happens.
