# Understanding the API Documentation

```{note}
:class: margin
If you are already familiar with this type of documentation you can skip this
page.
```

This page aims to help ShapePipe users and developers understand the
Application Programming Interface (API) documentation.

## What Are API Docs?

The API documentation is designed to communicate in clear way what each class
and function in the package does. For example, what inputs they expect, what
outputs they provide, and what the various options do. These pages are
automatically generated from docstrings (i.e. the sections starting and ending
with three double quotes `"""`) written in the code.

## Standard API Docs

```{note}
:class: margin
If an optional argument does not explicitly specify the default parameter value
then the user should expect that this means the default will be `None`, `''`,
`[]`, etc. depending on the input data type.
```

All the classes/functions include a short description of what they do. This is
followed by a `Parameters` section containing a bullet point list of the
expected input arguments. For each parameter you will see in brackets the
expected input type (e.g. `int`, `float`, `list`, etc.) followed by a brief
description of what the argument is for. Parameters listed as *optional* in the
brackets do not need to be provided and will default to some predefined value.

For functions that return objects a `Returns` section will follow the
`Parameters` section. This will provide a brief description of what is provided
by this function. Following `Returns` you will always find `Return type`, which
specifies the data type of the returned object.

If the function raises an exception under certain conditions you will find a
`Raises` section containing a bullet point list of the type of exceptions
raised. Each point is followed by a short description of the conditions that
will lead to the exception being raised.

Some objects may contain a `Notes` section providing further details for the
class/function. Some objects may also contain a `See Also` section that links
to another related object.

Finally, all objects include `[source]` button that will allow you to view the
code implementation for the given class/function. You can switch back to the
API docs by clicking the corresponding `[docs]` button.

## ShapePipe Module Docs

In addition to the standard API docs, we provide a description of each ShapePipe
module defined in the `__init__.py` file for each module package. On these
pages you will find a header that specifies the package author, the
`Parent module` (i.e. a module that should be run before the module in
  question), and a brief description of the inputs and outputs. This is
  followed by a description of the module. Finally, there is a section that
  defines all of the module-specific [config file options](configuration.md).
  Each option is listed with the exact name that should be added to the config
  file, followed by the expected value type. Underneath this is a brief
  description of what this option does. For example, the following option for
  a module called `my_module`:

**MY_OPTION : *int***

would be included in the config file as follows.

```ini
[MY_MODULE_RUNNER]
MY_OPTION = 1
```

Similarly to the standard API docs, parameters listed as *optional* do not need
to be provided and will default to some predefined value.
